#!/usr/bin/env python3

from edgembar import DiscoverEdges
import os
from pathlib import Path

#
# The output directory (where the edge xml input files are to be written
#
odir = Path(".")

#
# The format string describing the directory structure.
# The {edge} {env} {stage} {trial} {traj} {ene} placeholders are used
# to extract substrings from the path; only the {edge} {traj} and {ene}
# are absolutely required.  The {env} placeholder must be either 
# 'target' or 'reference', or you must supply the directory string
# with the target and reference optional arguments.
# If the {env} placeholder is missing, then
# 'complex' is assumed.
#
# Full example:
#    s = r"dats/{trial}/free_energy/{edge}_ambest/{env}/{stage}/efep_{traj}_{ene}.dat"
# Minimal example:
#    s = r"dats/{edge}/efep_{traj}_{ene}.dat"

s = r"dats/{edge}/{env}/{stage}/{trial}/efep_{traj}_{ene}.dat"

exclusions=None
# exclusions=["trial1"]
edges = DiscoverEdges(s,
                      target="complex",
                      reference="solvated",
                      exclude_trials=exclusions)

#
# In some instances, one may have computed a stage with lambda values
# going in reverse order relative to the thermodynamic path that leads
# from the reactants to the products. We can reverse the order of the
# files to effectively negate the free energy of each state (essentially
# treating the lambda 0 state as the lambda 1 state).
#
#for edge in edges:
#    for trial in edge.GetAllTrials():
#        if trial.stage.name == "STAGE":
#            trial.reverse()


#
# There may be an analytic free energy correction that is applied
# to some transformations, such as when Boresch restraints are used.
# We can add a constant onto a trial free energy.
#
# shifts = {}
# shifts["1a~1b"] = 1.0 # kcal/mol
# for edge in edges:
#     if not edge.name in shifts:
#         continue
#     for e in edge.GetEnvs():
#         if e.name != "target":
#             continue
#         for stage in e.stages:
#             if stage.name != "STAGE":
#                 continue
#             #   One could have different shifts for each trial
#             #for trial in stage.trials:
#             #    # if trial.name == "t1":
#             #    trial.SetShift( 1.0 )
#             #   Or one could apply a shift to all trials
#             stage.SetShift( shifts[edge.name] )
#
#

if not odir.is_dir():
    os.makedirs(odir)
    
for edge in edges:
    fname = odir / (edge.name + ".xml")
    edge.WriteXml( fname )

    
