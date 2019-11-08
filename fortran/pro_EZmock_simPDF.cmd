#!/bin/bash
##
## optional: energy policy tags
##
#@ job_type = MPICH
##@ class = general
##@ class = test
#@ class = fattest
#@ node = 1
###@ island_count= not needed for 
#@ total_tasks= 1
## other example
#@ wall_clock_limit = 2:00:00
##                     h  min  secs
#@ job_name = EZmock
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = /gpfs/scratch/h009z/di72yuc/zamock/fortran
#@ output = ./stdout/EZmock_d_$(jobid).out
#@ error = ./stdout/EZmock_d_$(jobid).err
#@ notification=always
#@ notify_user=r92222016@gmail.com
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
module unload mpi.ibm
module load mpi.intel
export OMP_NUM_THREADS=40
mpiexec -n 1 ./pro_EZmock_simPDF < EZmock_v0.input

