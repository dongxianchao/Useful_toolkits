#!/bin/bash
module load vasp/6.4.3-gpu

export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# Metropolis Monte Carlo sampling (restore to previous accepted step on reject)
set -euo pipefail

########## User parameters ##########
K=8.617333262e-5        # eV/K
T=300                  # K
VASP_CMD="srun -n 8 -c 32 --cpu-bind=cores --gpu-bind=none -G 8 vasp_std"   # VASP command
MAX_STEPS=0             # 0 => Infinity loop for 0, or others
# POSCAR format specifics (modify to your POSCAR)
HEAD_LINES=8            # POSCAR's header number
NGa=22                  # POSCAR's Ga number
NN=32                   # N atom number
NSc=10                  # Sc atom number
#####################################

# prepare directories
mkdir -p trajectory
mkdir -p dataset

# determine current indices (count subdirs)
trajectory_index=$(ls -1d trajectory/* 2>/dev/null | wc -l || true)
dataset_index=$(ls -1d dataset/* 2>/dev/null | wc -l || true)
trajectory_index=${trajectory_index:-0}
dataset_index=${dataset_index:-0}

# initial run (if OUTCAR exists, use it; otherwise run VASP once)
if [ -f OUTCAR ]; then
    E1=$(grep -i "TOTEN" OUTCAR | tail -n 1 | awk '{print $5}')
else
    echo "[INFO] No OUTCAR found    running initial VASP..."
    $VASP_CMD > vasp.out 2>&1
    E1=$(grep -i "TOTEN" OUTCAR | tail -n 1 | awk '{print $5}')
fi
echo "[INFO] Initial energy E1 = $E1 eV"

# backup current accepted POSCAR (this is the "previous accepted step" reference)
cp -f POSCAR POSCAR_backup

# main loop
step=0
while true; do
    if [ "$MAX_STEPS" -gt 0 ] && [ "$step" -ge "$MAX_STEPS" ]; then
        echo "[INFO] Reached MAX_STEPS=$MAX_STEPS. Exiting."
        break
    fi
    step=$((step+1))
    echo "================== Step $step =================="

    # save current POSCAR to allow restoration of this accepted state if reject
    cp -f POSCAR POSCAR_before_swap

    # randomly swap one Ga line and one Sc line in POSCAR
    # compute ranges
    ga_start=$((HEAD_LINES + 1))
    ga_end=$((HEAD_LINES + NGa))
    sc_start=$((HEAD_LINES + NGa + NN + 1))
    sc_end=$((HEAD_LINES + NGa + NN + NSc))

    # sanity checks
    if [ "$ga_start" -gt "$ga_end" ] || [ "$sc_start" -gt "$sc_end" ]; then
        echo "[ERROR] Wrong atom counts or HEAD_LINES. Check NGa/NN/NSc/HEAD_LINES."
        exit 1
    fi

    # pick random line numbers
    range_ga=$((ga_end - ga_start + 1))
    range_sc=$((sc_end - sc_start + 1))
    # Use $RANDOM; produce uniform selection
    line_ga=$((ga_start + RANDOM % range_ga))
    line_sc=$((sc_start + RANDOM % range_sc))
    echo "[INFO] swap lines: GA $line_ga <-> SC $line_sc"

    # get contents (use printf to preserve spaces)
    line_ga_content=$(sed -n "${line_ga}p" POSCAR)
    line_sc_content=$(sed -n "${line_sc}p" POSCAR)

    # create swapped POSCAR_tmp
    cp POSCAR POSCAR_tmp
    # replace lines (use awk to avoid sed delimiter issues)
    awk -v l1="$line_ga" -v t1="$line_sc_content" -v l2="$line_sc" -v t2="$line_ga_content" '
    NR==l1 { print t1; next }
    NR==l2 { print t2; next }
    { print }
    ' POSCAR_tmp > POSCAR_swapped
    mv POSCAR_swapped POSCAR

    # run VASP on new POSCAR
    echo "[INFO] Running VASP for trial structure..."
    $VASP_CMD > vasp.out 2>&1

    # read energy E2
    if ! grep -qi "TOTEN" OUTCAR; then
        echo "[ERROR] TOTEN not found in OUTCAR. Aborting."
        # restore POSCAR to backup then exit
        mv -f POSCAR_before_swap POSCAR
        exit 1
    fi
    E2=$(grep -i "TOTEN" OUTCAR | tail -n 1 | awk '{print $5}')
    echo "[INFO] Trial energy E2 = $E2 eV   (E1 = $E1 eV)"

    # compute p = exp(-(E2-E1)/(K*T))
    deltaE=$(awk -v e2="$E2" -v e1="$E1" 'BEGIN{print e2-e1}')
    exponent=$(awk -v dE="$deltaE" -v K="$K" -v T="$T" 'BEGIN{print -dE/(K*T)}')
    # compute p using awk's exp
    p=$(awk -v ex="$exponent" 'BEGIN{print exp(ex)}')
    # random x in [0,1)
    x=$(awk 'BEGIN{srand(); print rand()}')

    # increment dataset and save OUTCAR+CONTCAR
    dataset_index=$((dataset_index + 1))
    dst_dir="dataset/${dataset_index}"
    mkdir -p "$dst_dir"
    # copy if present
    [ -f OUTCAR ] && cp -f OUTCAR "$dst_dir/"
    [ -f CONTCAR ] && cp -f CONTCAR "$dst_dir/"

    echo "[INFO] E_diff=$deltaE , p=$p , x=$x"

    # acceptance check: accept if x < p
    acc="no"
    awk -v xx="$x" -v pp="$p" 'BEGIN{if (xx<pp) exit 0; else exit 1}' && acc="yes" || acc="no"

    if [ "$acc" = "yes" ]; then
        echo "[INFO] Move ACCEPTED."
        # create trajectory entry and copy outputs
        trajectory_index=$((trajectory_index + 1))
        tr_dir="trajectory/${trajectory_index}"
        mkdir -p "$tr_dir"
        [ -f OUTCAR ] && cp -f OUTCAR "$tr_dir/"
        [ -f CONTCAR ] && cp -f CONTCAR "$tr_dir/"

        # move CONTCAR -> POSCAR (update working POSCAR)
        if [ -f CONTCAR ]; then
            mv -f CONTCAR POSCAR
            # update backup to this accepted POSCAR
            cp -f POSCAR POSCAR_backup
        else
            echo "[WARN] CONTCAR not found on accept - POSCAR not updated."
        fi

        # update E1 to new accepted energy
        E1="$E2"
    else
        echo "[INFO] Move REJECTED. Restoring previous accepted POSCAR."
        # restore POSCAR from POSCAR_backup (previous accepted)
        if [ -f POSCAR_backup ]; then
            cp -f POSCAR_backup POSCAR
        else
            # fallback: restore pre-swap file
            mv -f POSCAR_before_swap POSCAR
        fi
    fi

    echo "-----------------------------------------------------"
done
