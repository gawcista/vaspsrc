import re
import os
import subprocess
import math
from numpy import array
from subprocess import Popen, PIPE, call

global NKPOINTS,Kpoints,TITLE

def getnum(s):
   return re.findall(r"[-+]?[0-9]*\.?[0-9]+",s)

def initialization():
   global NKPOINTS,Kpoints
   Kpoints = []
   NKPOINTS = 0

def read_IBZKPT():
   fin = open('scf/IBZKPT')
   lines = fin.readlines()
   fin.close

   global NKPOINTS,Kpoints,TITLE
   TITLE = lines[0]
   nk = int(getnum(lines[1])[0])
   NKPOINTS += nk
   for i in range(nk):
      kx = float(getnum(lines[3+i])[0])
      ky = float(getnum(lines[3+i])[1])
      kz = float(getnum(lines[3+i])[2])
      weight = float(getnum(lines[3+i])[3])
      Kpoints.append([kx,ky,kz,weight])

def read_kpoint_band():
   fin = open('KPOINTS.band')
   lines = fin.readlines()
   fin.close()

   global NKPOINTS,Kpoints
   num_points = int(getnum(lines[1])[0])
   num_lines = math.ceil((len(lines)-4)/3)
   NKPOINTS += num_points*num_lines
   for i in range(num_lines):
      kx0 = float(getnum(lines[4+3*i])[0])
      ky0 = float(getnum(lines[4+3*i])[1])
      kz0 = float(getnum(lines[4+3*i])[2])
      kx1 = float(getnum(lines[5+3*i])[0])
      ky1 = float(getnum(lines[5+3*i])[1])
      kz1 = float(getnum(lines[5+3*i])[2])
      k0 = array([kx0,ky0,kz0,0])
      k1 = array([kx1,ky1,kz1,0])
      for j in range(num_points):
         Kpoints.append(k0+(k1-k0)*j/(num_points-1))

def print_to_file(filename='KPOINTS.hse'):
   fout = open(filename,'w')
   fout.write(TITLE)
   fout.write(str(NKPOINTS)+'\n')
   fout.write('Reciprocal lattice\n')
   for value in Kpoints:
      fout.write(str(value[0])+' '+str(value[1])+' '+str(value[2])+' '+str(int(value[3]))+' '+'\n')
   fout.close()
initialization()
read_IBZKPT()
read_kpoint_band()
print_to_file()
