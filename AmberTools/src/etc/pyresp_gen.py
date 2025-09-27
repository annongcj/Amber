#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import argparse as ap
from pyrespgen.prepin import dumpin1, dumpin2

def main_run(espdat,Istage,IIstage,\
				ptype,dtype,\
				nmol,charge,\
				QWT1,QWT2,PWT1,PWT2,\
				EXC12,EXC13,\
				DEPTH,\
				verbose,\
				strategy):
	if Istage and IIstage:
		dumpin1(espdat,Istage,ptype,dtype,nmol,charge,QWT1,PWT1,EXC12,EXC13,DEPTH,verbose,strategy)
		dumpin2(espdat,IIstage,ptype,dtype,nmol,charge,QWT2,PWT2,EXC12,EXC13,DEPTH,0,strategy)
	elif len(Istage) !=0 and len(IIstage) == 0:
		dumpin1(espdat,Istage,ptype,dtype,nmol,charge,QWT1,PWT1,EXC12,EXC13,DEPTH,verbose,strategy)
	elif len(Istage) ==0 and len(IIstage) != 0:
		dumpin2(espdat,IIstage,ptype,dtype,nmol,charge,QWT2,PWT2,EXC12,EXC13,DEPTH,verbose,strategy)
	else:
		print('Please specify one or two stage input file name with option --Istage or --IIstage')


parser = ap.ArgumentParser(description='PyRESP Generation')
parser.add_argument('--espdat','-i', help='Input file of esp data',required=True)
parser.add_argument('--Istage','-f1',help='Output file for 1st stage',default='pyrespgen.1st')
parser.add_argument('--IIstage','-f2',help='Output file for 2nd stage',default='pyrespgen.2nd')
parser.add_argument('--ptype','-p',help='Polarization type: chg | ind | perm ',default='perm')
parser.add_argument('--dtype','-d',help='Damping Function type: (1) applequist | (2) tinker | (3) exponential | (4) linear | (5) pgm',default='pgm')
parser.add_argument('--nmol','-n',type=int,help='Number of conformations',default=1)
parser.add_argument('--charge','-q',type=int,help='Total charge for this structure or conformer',default=0)
parser.add_argument('--QWT1','-qwt1',type=float,help='Charge restraint 1st stage',default=0.0005)
parser.add_argument('--QWT2','-qwt2',type=float,help='Charge restraint 2nd stage',default=0.001)
parser.add_argument('--PWT1','-pwt1',type=float,help='Permanent dipoles restraint 1st stage',default=0.0005)
parser.add_argument('--PWT2','-pwt2',type=float,help='Permanent dipoles restraint 2nd stage',default=0.001)
parser.add_argument('--EXC12','-exc12',type=int,help='include (0) or exclude (1) 1-2 interaction',default=0)
parser.add_argument('--EXC13','-exc13',type=int,help='include (0) or exclude (1) 1-3 interaction',default=0)
parser.add_argument('--DEPTH','-depth',type=int,help='Maximum depth for searching equilvalance atoms',default=3)
parser.add_argument('--verbose','-v',type=int,help='Print verbose information',default=0)
parser.add_argument('--strategy','-strategy',type=int,help='Choose Strategy for pGM-perm',default=2)


args = parser.parse_args()

if __name__ == '__main__':
	try:
		main_run(args.espdat,args.Istage,args.IIstage,\
						args.ptype,args.dtype,\
						args.nmol,args.charge,\
						args.QWT1,args.QWT2,args.PWT1,args.PWT2,\
						args.EXC12,args.EXC13,\
						args.DEPTH,\
						args.verbose,\
						args.strategy)
	except Exception as e:
		print(e)
		
