# Requirements

Tested to work on Ubuntu 20.04 LTS with:

python3
cython
openmpi (or similar)
mpi4py
numpy
scipy

# Files
- cyclobraid.pyx        :  Cython interface file
- setup-cyclobraid.py   :  Install script for cyclobraid, 
                           Read comments at top for details 
                           on how to run and tweak for your 
                           platform
- cyclobraid_driver.py  :  Driver file for running test cases


# Example

You may have to first update $PYTHONPATH with the cyclobraid install location, e.g., 

    $ export PYTHONPATH="$HOME/.local/lib/python3.7:$PYTHONPATH"

or

    $ export PYTHONPATH="$HOME/.local/lib/python3.7/site-packages:$PYTHONPATH"

Two simple example runstrings are 

    $ python3 cyclobraid_driver.py --Nt=50 --tf=1.0 --FCF=0 --m=50 --maxiter=10 --maxlevels=2

Or, to use MPI:

    $ mpirun -np 8 python3 cyclobraid_driver.py --Nt=50 --tf=1.0  --FCF=0 --m=50 --maxiter=10 --maxlevels=2


# Reproducing data from Masters Thesis

Nick: Can you add runstrings here tah reproduce an entry from each Table in your thesis? 


