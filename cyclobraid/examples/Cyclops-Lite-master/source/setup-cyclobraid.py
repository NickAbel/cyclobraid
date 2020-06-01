from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os
# import numpy

##
# Tested on High Sierra with Homebrew, and Ubuntu LTS.
#
# To Install,
#
# 1) Make sure that library_dirs and include_dirs point to the 
#    location of "braid"
#
# 2) Make sure to compile braid in the director
#    $ ./cyclobraid/braid
#    with
#    $ make clean; make
#
# 3) Type (using whatever install location you want)
#
#    $ python3 ex_01-setup.py install --prefix=$HOME/.local
#
#    Note that you may have to tweak the compilers and flags.
#    Some comments on this are below.
#   
#    For example, to wipe the previous compile and do a local install, 
#
#    $ rm -rf build
#    $ python3 setup-cyclobraid.py install --prefix=$HOME/.local
#
#    Sometimes, the PYTHONPATH environment variable needs to be wiped
#    of anything Python2 specific for the compile to work, e.g., 
#       $ export PYTHONPATH=""


# To Run, 
#
# 1) Make sure that the install directory and the location of MPI4PY 
#    is in your PYTHONPATH, e.g.,
# 
#     export PYTHONPATH="$HOME/.local/lib/python3.7"
#
# 2) Run cyclobraid_driver.py (see README)


##
# Other notes:
#  1) Some systems may need to find Numpy headers, which are located in 
#     include_dirs=["../braid", numpy.get_include()],
#  2) Some compilers may require "-fPIC" to be added to extra_compile_args
#
##





os.environ["CC"] = "mpicc"
os.environ["LDSHARED"] = "mpicc -shared"  # May need to comment out for High Sierra with Homebrew

cyclobraid_extension = Extension(
    name="cyclobraid",
    sources=["cyclobraid.pyx"],
    libraries=["braid"],
    library_dirs=["../../../braid"],
    include_dirs=["../../../braid"], ##, numpy.get_include()],
    extra_compile_args=["-Wno-incompatible-pointer-types", "-Wno-unused-function"]   
)
setup(
    name="cyclobraid",
    ext_modules=cythonize([cyclobraid_extension], language_level = "3")
)


