set -e

# Get to correct directory and configure VAC
cd ~/stusac/sac/src/
make cleanall
./setvac -d=11 -phi=0 -z=0 -g=264 -p=mhd -u=viscosity -on=cd,mpi -off=tvd

# Build
make vacini;
make vac;

# Change dirs again and run initiation script
cd ../
python inipar/strongshock-test11.py

# Run VAC
#./vac < par/strongshock-test11
cp par/strongshock-test11 vac.par
mpirun -np 1 vac

# Return to test directory
cd tests/
