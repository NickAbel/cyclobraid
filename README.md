# cyclobraid
This code repository consists of the XBraid parallel-in-time library (written in C),
the asymptotic parareal 1D rotating shallow water solver Cyclops-Lite (Python), and
an interface written in Cython that enables numerical routines in Cyclops-Lite to
be called by, and used with, the XBraid library. This interface is referred to as
Cyclobraid.


# Licenses
Licenses for XBraid are located in the cyclobraid/ directory. Licenses for
Cyclops-Lite are located in the directory 'cyclobraid/examples/Cyclops-Lite-master'
directory.


# Running and Compiling the Code
See the README in the directory 'cyclobraid/examples/Cyclops-Lite-master/source'


# Directory Structure

1. [here] -> cyclobraid : Main directory, holds an un-tarred version of Braid.
                          All the folders here correspond to the normal Braid
                          directory structure, e.g.,   
                       
    1. [here] -> cyclobraid -> braid, holds the source code for braid
                       
    and 
                       
    2. [here] -> cyclobraid -> examples, holds the standard braid examples

    3. The cyclobraid code is written as an "example" of Braid 


2. [here] -> cyclobraid -> examples -> Cyclops-Lite-Master : Is an un-tarring of Adam Peddles
                                                          Cyclops code, modified to work 
                                                          with Braid


3. [here] -> cyclobraid -> examples -> Cyclops-Lite-Master -> source : Holds the cyclobraid 
                                                                       Cython interface twixt
                                                                       Braid and Cyclops

    See the README in this 
    directory for compiling
    and running details

