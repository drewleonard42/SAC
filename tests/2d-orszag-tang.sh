set -e

# Get to correct directory and configure VAC
cd ~/stusac/sac/src/
make cleanall
./setvac -d=22 -phi=0 -z=0 -g=264,264 -p=mhd -u=nul -on=cd,mpi, -off=tvd

# Build
make vac;

# Change dirs again and run initiation script
cd ../
python inipar/orszagtang-test22.py

# Run VAC
#./vac < par/orszagtang-test22
cp par/orszagtang-test22 vac.par
mpirun -np 4 vac

# Gather results
./distribution data/orszagtang-test22_np0202.out data/orszagtang-test22.out
