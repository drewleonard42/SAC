#!/bin/sh
#$ -l h_rt=01:00:00
#$ -cwd
#$ -pe openmpi-ib 16
#$ -l arch=intel*
#$ -l mem=8G
#$ -j y
#$ -N .out-2-.h5
#$ -P mhd
#$ -q mhd.q

cat $PE_HOSTFILE

/usr/local/mpi/pgi/openmpi/1.4.4/bin/mpirun sac_hdf5_mpi2


