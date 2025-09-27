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
       try:
           for line in f:
               if line[:5] == "#U_ev":
                  U = float(line.split()[-1])
               else:
                  L = line.split()
                  if len(L) == 4:
                     at_count += 1
                  if len(L) == 1:
                     nAts  = int(line)
       except Exception as e:
           print("problem parsing file: %s\nline:\n%s" % (fname, line))
           raise e
           
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
    helptext = \
    """Simple test harness for sander energy eval.

       This version is intended for solvation energy calculations, so not 
       soooo simple: need three directories full of xyz files with 1:1:whatever correspondance
       directories must be passed as arguments, in order:

       ./precache_forcefield_enes_nolibsander_subtractWaterShells.py  SOLUTE_SOLVENT SOLVENT SOLUTE lambda

       The fourth argument is a lambda value on 0..1 which weights the solvated [0] versus unsolvated [1]
       versions of the structures.

       the script should loop over the directories, creating triples of structures, and then 
       save a database with only two sets of files:

       file1_1.xyz  Qene Eclass 1-lambda
       file2_1.xyz  Qene Eclass 1-lambda 
       . 
       .
       file1_2.xyz  Qene Eclass lambda
       file2_2.xyz  Qene Eclass lambda
       . 
       .

       Energies in the first set of files are the values (SOLUTE_SOLVENT - SOLVENT) for the given struc.

       Energies in the second set of files are the values SOLUTE only, so as to characterise the system
       genuinely in vacuum, where it will typically be much stiffer compared to the solvated
       structure.   
    """
    
    
    
    if len(sys.argv) < 5:
        raise ValueError(\
           """
           %s
          
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
           
        """ % helptext)
        
    print("caching enes, arguments are:", sys.argv)
    solute_solvent_dir = sys.argv[1]
    shell_dir          = sys.argv[2]
    gasphase_dir       = sys.argv[3]
    lambda_value       = float(sys.argv[4])


    fnames_lists = []
    qene_lists   = []
    cene_lists   = []    
    for db_dir in [solute_solvent_dir, shell_dir, gasphase_dir]:
        fnames, enes_ev = getAimsOutputStrucs( db_dir )
        fnames_lists.append( fnames )
        qene_lists.append( enes_ev )

    if len(qene_lists[0]) != len(qene_lists[1]):
        raise ValueError( "different numbers of quentum enes, %i versus %i for %s vs %s" %\
                      (len(qene_lists[0]),len(qene_lists[1]), solute_solvent_dir, shell_dir)  )


    ##Do (and cache) classical energy evaluations
    sys_datasets = []
    for ii, db_dir in enumerate([solute_solvent_dir, shell_dir, gasphase_dir]):
        fnames = fnames_lists[ii]
        qenes  = qene_lists[ii]
        db_subdir_path  = os.path.dirname( fnames[0] )
        cache_fname     = os.path.join( db_subdir_path, "cache_ffenes.txt" )
        f, q, c         = [], [], []

        try:
            f_cache = open( cache_fname, "r" )
            for line in f_cache:
                print("read cache line: %s" % line[:len(line)-1])
                file_fname, Qene_kcal, E_kcal = line.split()
                f.append(file_fname)
                q.append(float(Qene_kcal))
                c.append(float(E_kcal))
            f_cache.close()
            sys_datasets.append( [f,q,c]  )
            continue
        except:
            print("could not find cache %s for reading, so generating energies" % cache_fname)


        f_cache         = open( cache_fname, "w" )
        print(cache_fname)
    
        top = glob.glob(db_subdir_path+"/*.top")
        if len(top) != 1:
                 raise ValueError("Expecting exactly one parmtop in the DB_subdir, got:", top)
        else:
                 top = top[0]
    
        enes_class_kcal = []
        for fname, Uev in zip(fnames, qenes):
            usetmp = "."
       # usetmp = None 
            E_kcal =\
                   eval_ene_sander_syscall( xyzFile=fname, topFile=top, sander_exe_path="sander", override_temp=usetmp )
            if E_kcal is None:
                   raise ValueError("Could not process file based on %s" % fname)
            Qene_kcal =  Uev * EV_TO_KCAL_MOL
            f_fname   =  os.path.basename( fname )
            f_cache.write("%s %.16f %.16f\n" % (os.path.join(db_subdir_path,f_fname), Qene_kcal, E_kcal) )
            f.append( os.path.join(db_subdir_path,f_fname) )
            q.append( Qene_kcal )
            c.append( E_kcal )

        f_cache.close()
        sys_datasets.append( [f,q,c]  )     
     
    DB_outfile     = open("DB_lambda_%s.txt" % sys.argv[4], "w")
    n_pts_solvated = len( sys_datasets[0][0] )

    for i in range(n_pts_solvated):
         DB_outfile.write("%s %.18f %.18f %.4f -1\n" % \
              (sys_datasets[0][0][i], \
               sys_datasets[0][1][i] - sys_datasets[1][1][i],\
               sys_datasets[0][2][i] - sys_datasets[1][2][i],\
               lambda_value))



    print("wrote %i solvated datapoints" % n_pts_solvated)
    n_pts_dry = len( sys_datasets[2][0] )
    for i in range(n_pts_dry):
         DB_outfile.write("%s %.18f %.18f %.4f 1\n" % \
              (sys_datasets[2][0][i],\
               sys_datasets[2][1][i],\
               sys_datasets[2][2][i],\
               1.0-lambda_value))
    print("wrote %i dry datapoints" % n_pts_dry)
    DB_outfile.close()





                                       
    
