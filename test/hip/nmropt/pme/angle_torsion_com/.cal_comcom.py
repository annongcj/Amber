import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
import pytraj as pt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import os
#/mnt/raidc2/haichit/anaconda3/bin/python

def cal_dist(top):
  traj =  pt.load(top)
  traj.autoimage()
  com1_com2 = pt.distance(traj, '(:10-40,60-90)&(@CA) (:109-139,159-189)&(@CA)')
  Yarray = com1_com2
  print ("COM mon - mon distance:")
  print (Yarray)

def cal_dist1(top):
  traj =  pt.load(top)
  traj.autoimage()
  com1_com2 = pt.distance(traj, ':31,75,85,89@N,CA,C,O :130,174,184,188@N,CA,C,O')
  Yarray = com1_com2
  print ("COM res 317585 89 distance:")
  print (Yarray)


def cal_hand(top):
  traj =  pt.load(top)
  traj.autoimage()
  #handed_dih = pt.dihedral(traj, ':50,51, :87-92 :186-191 :151,152')
  handed_dih = pt.dihedral(traj, ':48,49,52,53@CA,N,C :87-92@CA,N,C :186-191@CA,N,C :147,148,151,152@CA,N,C')
  Yarray = handed_dih
  print ("handedness dihedral:")
  print (Yarray)

def cal_ang(top):
  traj =  pt.load(top)
  traj.autoimage()
  com_Ang = pt.angle(traj, ':32,57,58,75,76@CA,N,C :2,3,96,97,101,102,195,196@CA,N,C :131,156,157,174,175@CA,N,C')

  Yarray = com_Ang
  print ("COM Ang:")
  print (Yarray)

# openHIVPR is an open structure from GBneck2 14sbonlysc simualtion 

list_pdb = ['mypdb', 'leap_1hhp', 'openHIVPR']
for i in list_pdb:
  print (i)
  cal_dist('./%s.pdb' %(i)) 
  cal_dist1('./%s.pdb' %(i)) 
  cal_hand('./%s.pdb' %(i)) 
  cal_ang('./%s.pdb' %(i))

