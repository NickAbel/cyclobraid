# Cyclobraid

# Requirements

Tested to work on Ubuntu 20.04 LTS with:

 - python3
 - cython
 - openmpi
 - mpi4py
 - numpy
 - scipy

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

# Usage Details

 - Nt        : Number of time points on the coarsest grid.
 - tf        : Final time (initial time set to 0.)
 - FCF       : Set to 1 for FCF-relaxation, 0 for F-relaxation (as in parareal,) 2 for FCFCF, etc.
 - m         : Coarsening factors. To specify for multilevel, can specify multiple coarsening factors in level order. Specifying only one --m argument will set a global coarsening factor for all levels.
 - maxiter   : Maximum number of Braid iterations. Pre-empted by halting tolerance.
 - maxlevels : Number of Braid levels. Set equal to 2, and FCF to 0, for a solve equivalent to parareal.

# More Example Runstrings


The example given above, shown here:

    $ python3 cyclobraid_driver.py --Nt=50 --tf=1.0 --FCF=0 --m=50 --maxiter=10 --maxlevels=2

Runs a two-level solve with F-relaxation. The time domain is [0,1] with 50 coarse time points, a coarsening factor of 50 hence 2500 fine grid time points, and 10 maximum iterations.

An example multilevel runstring is given by:

    $ python3 cyclobraid_driver.py --Nt=8 --tf=1. --FCF=0 --m=1 --m=4 --m=8 --maxlevels=3 --maxiter=8  

In particular, this gives a three-level time grid hierarchy of [256, 64, 8] time points.

# Reproducing data from Masters Thesis

Problem-specific parameters for the RSWE are set initially in cyclops_control.py and can be accessed in the control object.

The parameters set in cyclops_control.py are as follows:
  - Scale separation parameter epsilon is stored in control['epsilon']
  - Rossby radius of deformation is stored in control['deformation_radius']
  - Dissipation parameter mu is stored as control['mu']
  - Other parameters, such as start time (control['start_time']), spatial domain size and resolution (control['Lx'] and control['Nx'], respectively)
 
Other problem-specific parameters using the control object are set in cyclobraid.pyx:
  - Averaging window size eta and quadrature points M are stored in control['HMM_T0'] and control['HMM_M_bar']. They are set in cyclobraid.pyx and can be recomputed for multilevel or alternating window solves using the function recompute_T0_M() in cyclobraid.pyx
  - Final time t_f is stored in control['final_time'] and is set by command line argument in braid_init_py() in cyclobraid.pyx.
 
Note that in Cyclobraid, the halting tolerance is specified using the function call braid_SetAbsTol() inside braid_init_py() in cyclobraid.pyx. The tolerance in the control object for Cyclops is not used to determine halting in Cyclobraid.

The behavior of the solve is given in the my_step() function in cyclobraid.pyx. An if statement determines whether to use the fine solver or coarse asymptotic solver depending on the current level that my_step() is being called in. The recompute_T0_M() function can be called to recompute the averaging window parameters eta and M before calling the coarse solver. 

To perform reuse of prior values of asymptotic nonlinear products, the 'coarse_propagator_reuse' function can be called. The coarse_propagator_reuse function calls the strang_splitting_reuse function in RSWE_direct.py which accesses dicts containing nonlinear asymptotic terms computed in prior iterations to return a stale value. Currently, strang_splitting_reuse is set up to use the value from the previous iteration (control['iter']-1), but this can be changed to reuse any extant values in the dicts.

Code for wall timing calls to my_step() is available in my_step().
