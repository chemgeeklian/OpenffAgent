#!/bin/bash -l
#PBS -l select=1:system=polaris
#PBS -l walltime=24:00:00
#PBS -q preemptable
#PBS -A FoundEpidem
#PBS -l filesystems=home:eagle
#PBS -N boltz_4x_parallel

module load conda 

conda activate openff312
cd /home/xlian/OpenffAgent/md_test

export CUDA_VISIBLE_DEVICES=0,1
python sim_test.py --tmp_dir ../tmp &> sim_test_ec03_1.log &

export CUDA_VISIBLE_DEVICES=2,3
python sim_test.py --tmp_dir ../wt_tmp &> sim_test_wt.log &

wait