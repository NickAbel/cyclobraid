import cyclobraid          
import argparse            
from ctypes import *       

parser = argparse.ArgumentParser(description='Nt = Coarse grid time pts.\n Nx = Points in space (Fourier modes.)\n tf = Final time\n eps = Scale separation parameter for RSW equations\n FCF = Set to 1 for FCF, 0 for F-relaxation\n m = coarsening factor\n')

parser.add_argument("--Nt",type=int)
parser.add_argument("--tf",type=float)
parser.add_argument("--FCF",type=int)
parser.add_argument("--m",action='append',type=int)
parser.add_argument("--maxiter",type=int)
parser.add_argument("--maxlevels",type=int)
parser.add_argument("--skip",type=int)
parser.add_argument("--avg_lvl_switch",type=int)
parser.add_argument("--l2_coeff_t0",type=float)

args = parser.parse_args() 

cyclobraid.braid_init_py(args.Nt,args.tf,args.FCF,args.m,args.maxiter,args.maxlevels,args.skip,args.avg_lvl_switch,args.l2_coeff_t0)
