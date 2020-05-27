import cyclobraid          
import argparse            
from ctypes import *       

parser = argparse.ArgumentParser(description='Nt = Coarse grid time pts.\n tf: Final time (initial time set to 0.)\n FCF: Set to 1 for FCF, 0 for F-relaxation, 2 for FCFCF, etc. \n m: coarsening factors. To specify for multilevel, can specify multiple coarsening factors in level order. One --m argument will set the global coarsening factor.\n maxiter: Maximum number of Braid iterations. Pre-empted by halting tolerance. \n maxlevels: Number of Braid levels. Set equal to 2, and FCF to 0, for parareal. ')

parser.add_argument("--Nt",type=int)
parser.add_argument("--tf",type=float)
parser.add_argument("--FCF",type=int)
parser.add_argument("--m",action='append',type=int)
parser.add_argument("--maxiter",type=int)
parser.add_argument("--maxlevels",type=int)

args = parser.parse_args() 

cyclobraid.braid_init_py(args.Nt,args.tf,args.FCF,args.m,args.maxiter,args.maxlevels)
