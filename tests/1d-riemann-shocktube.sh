set -e

# Get to correct directory and configure VAC
cd ~/stusac/sac/src/
make cleanall
./setvac -d=11 -phi=0 -z=0 -g=264 -p=mhd -u=viscosity -on=cd,mpi -off=resist

# Build
make vac;

# Change dirs again and run initiation script
cd ../
python inipar/riemann-test11.py

# Run VAC
cp par/riemann-test11 vac.par
mpirun -np 1 vac
#./vac < par/riemann-test11

#./distribution data/riemann-test11_np01.out data/riemann-test11.out 
