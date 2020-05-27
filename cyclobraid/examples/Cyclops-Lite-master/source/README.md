# Requirements

Tested to work on Ubuntu 20.04 LTS with:

python3
cython
openmpi (or similar)
mpi4py
numpy
scipy

# Example

python3 cyclobraid_driver.py --Nt=50 --tf=1.0 --FCF=0 --m=50 --maxiter=10 --maxlevels=2

Or, to use MPI:

mpirun -np 8 python3 cyclobraid_driver.py --Nt=50 --tf=1.0  --FCF=0 --m=50 --maxiter=10 --maxlevels=2
