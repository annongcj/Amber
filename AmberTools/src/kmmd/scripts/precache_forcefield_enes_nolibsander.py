#!/usr/bin/env python3

'''
This script is intended to sweep a database of structures in .xyz format (fhi-aims quantum code output)
and evaluate classical forcefield energies based on parameters in an AMBER parmtop file.

The database of cached classical energies is then picked up by the SLUK code for quantum-assisted 
classical MD.
'''
import glob
import sys, getopt, subprocess
import numpy as np
import os
import tempfile
import math

EV_TO_KCAL_MOL=23.06035

sander_inpstr = \
""" sander gas phase
&cntrl
ntx             = 1,        !5:Read coordinates and velocities
cut             = 99,       !Cut off (direct) interactions after xxA.
nstlim          = 0,        !Number of MD steps
igb             = 6,        !0 or 6 : periodic/infinite vacuum.
/
"""

def xyzToCrd( xyzFile, crdFile ):
    f     = open( xyzFile, "r" )
    lines = f.readlines()
    f.close()
    N_ats = int(lines[0])
    
    ##format in theory cannot have more than 99999 ats. Real code probably forgives having more.
    column_count = 0
    g = open( crdFile, "w" )
    if len(lines[1]) > 20:
       g.write(lines[1][:20]+"\n") ##amber restart file format says comment can only be 20 chars.
    else:
       g.write(lines[1])   
       
    g.write("%i\n" % N_ats)  
    for line in lines[2:]:
       L     = line.split()
       if len(L) == 4:
          atcrds = (float(L[1]), float(L[2]), float(L[3]))
          for crd in atcrds:
             g.write("%12.7f" % crd)
             column_count += 1
             if column_count == 6:
                g.write("\n")
                column_count = 0
    if column_count != 0:
       g.write("\n")
    g.close()
          
     

def eval_ene_sander_syscall( xyzFile, topFile, sander_exe_path="sander", override_temp = None ):

    '''
    Make a system call to sander (in a tmp directory), get the energy, and clean up.
    '''
    
    workdir = tempfile.TemporaryDirectory()
    wd_name = workdir.name
    
    
    if override_temp is not None:
        wd_name = override_temp
    
    in_name       = os.path.join(wd_name,"md")
    sander_in     = open(in_name,"w")
    sander_in.write(sander_inpstr)
    sander_in.close()
    
    ##reformat aims output for amber input.
    xyzToCrd( xyzFile, "%s.crd" % in_name )
    
    
    cmd = [sander_exe_path, "-O",\
                 "-i %s" % in_name,\
                 "-o %s.out" % in_name,\
                 "-p %s" % topFile,\
                 "-r %s.rst" % in_name,\
                 "-inf %s.mdinfo" % in_name,\
                 "-c %s.crd" % in_name]
    
    try:
       returned = subprocess.check_output(" ".join(cmd),\
                                           stderr=subprocess.STDOUT, shell=True)
    except Exception as e:
       print("could not execute sander command in temp directory: _%s_" % cmd)
       print( e )
    
    f = open("%s.out" % in_name, "r")
    for line in f:
       if "EPtot      =" in line:
          E_kcal = float(line.split()[-1])
          f.close()
          return E_kcal 
    f.close()  
    return
    
def getAimsOutputStrucs( DB_dir ):
    '''
    Scan a directory for aims outputs and harvest coordinates for amber calc.
    '''
    xyzfiles = glob.glob(DB_dir+"/**/*.xyz", recursive=True)
    if len( xyzfiles ) < 1:
       search = DB_dir+"/*.xyz"
       print("globbing %s" % search )
       xyzfiles = glob.glob(search)
       print( xyzfiles )
    
    
    fnames_checked  = []
    enes_ev         = []
    for fname in xyzfiles:
    
       f = open(fname, 'r')
       at_count = 0
       U = None
       for line in f:
           if line[:5] == "#U_ev":
               U = float(line.split()[-1])
           else:
               L = line.split()
               if len(L) == 4:
                   at_count += 1
               if len(L) == 1:
                   nAts  = int(line)
       f.close()
       if at_count == nAts and U is not None:
           fnames_checked.append(fname)
           enes_ev.append(U)
       else:
           print("filname %s : problems to parse atoms and energy in ev" % fname)
    
    print("got a quantum calculation DB of %i files in %s" %\
             (len(fnames_checked), DB_dir))
    
    return fnames_checked, enes_ev               
                   
    
    
    
    
if __name__ == "__main__":
    """Simple test harness for sander energy eval."""
    
    
    
    if len(sys.argv) < 2:
        raise ValueError(\
        """Require argument: a directory with some quantum calculations.
           Individual file format is:
           
           int_num_atoms
           #U_ev: float_energy_electron_volts
           char_element float_x float_y float_z
           char_element float_x float_y float_z
           char_element float_x float_y float_z
           char_element float_x float_y float_z
           .
           .
           
           Files should end with suffix ".xyz"
           
        """)
        
    print("caching enes, arguments are:", sys.argv)
    db_dir = sys.argv[1]
    if len( sys.argv ) == 3:
        fix_weight = float(sys.argv[2])
    else:
        fix_weight = 1.0
    
    
    fnames, enes_ev = getAimsOutputStrucs( db_dir )
    
    db_subdir_path  = os.path.dirname( fnames[0] )
    cache_fname     = os.path.join( db_subdir_path, "cache_ffenes.txt" )
    f_cache         = open( cache_fname, "w" )
    print(cache_fname)
    
    top = glob.glob(db_subdir_path+"/*.top")
    if len(top) != 1:
        raise ValueError("Expecting exactly one parmtop in the DB_subdir, got:", top)
    else:
        top = top[0]
    
    enes_class_kcal = []
    for f, Uev in zip(fnames, enes_ev):
        usetmp = "."
       # usetmp = None 
        E_kcal =\
           eval_ene_sander_syscall( xyzFile=f, topFile=top, sander_exe_path="sander", override_temp=usetmp )
        if E_kcal is None:
           raise ValueError("Could not process file based on %s" % f)
        enes_class_kcal.append( E_kcal ) 
        print(f, Uev, E_kcal)
        Qene_kcal =  Uev * EV_TO_KCAL_MOL
        f_fname   =  os.path.basename( f )
        f_cache.write("%s %.16f %.5f %.5f\n" % (os.path.join(db_subdir_path,f_fname), Qene_kcal, E_kcal, fix_weight) )
    
       
       
       
                                           
    
