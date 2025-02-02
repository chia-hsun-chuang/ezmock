#!/bin/bash
#SBATCH -o ../stdout/ezmock.%j.%N.out
#SBATCH -e ../stdout/ezmock.%j.%N.err
#SBATCH --get-user-env
### #SBATCH --clusters=mncp
###SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=all
#SBATCH --mail-user=r92222016@gmail.com
#SBATCH --export=NONE
#SBATCH --time=300:00:00
##SBATCH --cpus-per-task=1
##SBATCH -A serial
#SBATCH --cpus-per-task=16
#SBATCH -A 16cores
##
## force to run on brutus
##SBATCH  --mem=900000

source /opt/intel/composerxe/bin/compilervars.sh intel64

export OMP_NUM_THREADS=16
ulimit -c 10000
ulimit -s unlimited
time ./pro_EZmock_pk_cf < EZmock_v0.input.eBOSS_ELG

