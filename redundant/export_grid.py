import csv
import numpy as np
import matplotlib.pyplot as plt
from math import *
from collections import defaultdict
import matplotlib.cm as cm
import sys
import os
import matplotlib.animation as animation
from scipy import interpolate
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'



with open("./p17000_end.csv", 'r') as file:
  csvreader = csv.reader(file)
  i=0
  for line in csvreader:
    i=i+1
    

import csv
import numpy
reader = csv.reader(open("./p17000_end.csv", "r"), delimiter=",")
x = list(reader)
result = numpy.array(x).astype("float")
result = np.transpose(result) 




x = np.arange(0, 334, 1)
y = np.arange(0, 20,  1)
xx, yy = np.meshgrid(x, y)
z = result
f = interpolate.interp2d(x, y, z, kind='cubic')



xnew = np.arange(0, 334,  334/5000)
ynew = np.arange(0, 20,  1/3)
znew = f(xnew, ynew)


print(np.shape(znew)) 


np.savetxt("p17000_1st_cycle_10ms.csv", np.transpose(znew), delimiter=",")


with open("./rho17000_end.csv", 'r') as file:
  csvreader = csv.reader(file)
  i=0
  for line in csvreader:
    i=i+1
    

import csv
import numpy
reader = csv.reader(open("./rho17000_end.csv", "r"), delimiter=",")
x = list(reader)
result = numpy.array(x).astype("float")
result = np.transpose(result) 




x = np.arange(0, 334, 1)
y = np.arange(0, 20,  1)
xx, yy = np.meshgrid(x, y)
z = result
f = interpolate.interp2d(x, y, z, kind='cubic')



xnew = np.arange(0, 334, 334/5000)
ynew = np.arange(0, 20,  1/3)
znew = f(xnew, ynew)


print(np.shape(znew)) 


np.savetxt("rho17000_1st_cycle_10ms.csv", np.transpose(znew), delimiter=",")


with open("./u17000_end.csv", 'r') as file:
  csvreader = csv.reader(file)
  i=0
  for line in csvreader:
    i=i+1
    

import csv
import numpy
reader = csv.reader(open("./u17000_end.csv", "r"), delimiter=",")
x = list(reader)
result = numpy.array(x).astype("float")
result = np.transpose(result) 




x = np.arange(0, 334, 1)
y = np.arange(0, 20,  1)
xx, yy = np.meshgrid(x, y)
z = result
f = interpolate.interp2d(x, y, z, kind='cubic')



xnew = np.arange(0, 334, 334/5000)
ynew = np.arange(0, 20,  1/3)
znew = f(xnew, ynew)


print(np.shape(znew)) 


np.savetxt("u17000_1st_cycle_10ms.csv", np.transpose(znew), delimiter=",")


with open("./v17000_end.csv", 'r') as file:
  csvreader = csv.reader(file)
  i=0
  for line in csvreader:
    i=i+1
    

import csv
import numpy
reader = csv.reader(open("./v17000_end.csv", "r"), delimiter=",")
x = list(reader)
result = numpy.array(x).astype("float")
result = np.transpose(result) 




x = np.arange(0, 334, 1)
y = np.arange(0, 20,  1)
xx, yy = np.meshgrid(x, y)
z = result
f = interpolate.interp2d(x, y, z, kind='cubic')



xnew = np.arange(0, 334, 334/5000)
ynew = np.arange(0, 20,  1/3)
znew = f(xnew, ynew)


print(np.shape(znew)) 


np.savetxt("v17000_1st_cycle_10ms.csv", np.transpose(znew), delimiter=",")





