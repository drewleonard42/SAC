set -e

# Get to correct directory and configure VAC
cd ~/stusac/sac/src/
make cleanall
./setvac -d=12 -phi=0 -z=0 -g=808 -p=mhd -u=viscosity -on=cd,mpi -off=tvd

# Build
make vac;

# Change dirs again and run initiation script
cd ../
python inipar/briowu-test12.py

# Run VAC
#./vac < par/briowu-test12
cp par/briowu-test12 vac.par
mpirun -np 1 vac

# Return to test directory
cd tests/
