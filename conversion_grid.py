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

Name=84400

with open("./p"+str(Name)+".csv", 'r') as file:
  csvreader = csv.reader(file)
  i=0
  for line in csvreader:
    i=i+1
    

import csv
import numpy
reader = csv.reader(open("./p"+str(Name)+".csv", "r"), delimiter=",")
x = list(reader)
result = numpy.array(x).astype("float")
result = np.transpose(result) 

print(np.shape(result)[0])


x = np.arange(0, np.shape(result)[1], 1)
y = np.arange(0, np.shape(result)[0],  1)
xx, yy = np.meshgrid(x, y)
z = result
f = interpolate.interp2d(x, y, z, kind='cubic')



xnew = np.arange(0, np.shape(result)[1], 5)
ynew = np.arange(0, np.shape(result)[0],  1)
znew = f(xnew, ynew)


print(np.shape(znew)) 


np.savetxt("p_new.csv", np.transpose(znew), delimiter=",")


with open("./rho"+str(Name)+".csv", 'r') as file:
  csvreader = csv.reader(file)
  i=0
  for line in csvreader:
    i=i+1
    

import csv
import numpy
reader = csv.reader(open("./rho"+str(Name)+".csv", "r"), delimiter=",")
x = list(reader)
result = numpy.array(x).astype("float")
result = np.transpose(result) 




x = np.arange(0, np.shape(result)[1], 1)
y = np.arange(0, np.shape(result)[0],  1)
xx, yy = np.meshgrid(x, y)
z = result
f = interpolate.interp2d(x, y, z, kind='cubic')



xnew = np.arange(0, np.shape(result)[1], 5)
ynew = np.arange(0, np.shape(result)[0],  1)
znew = f(xnew, ynew)


print(np.shape(znew)) 


np.savetxt("rho_new.csv", np.transpose(znew), delimiter=",")


with open("./u"+str(Name)+".csv", 'r') as file:
  csvreader = csv.reader(file)
  i=0
  for line in csvreader:
    i=i+1
    

import csv
import numpy
reader = csv.reader(open("./u"+str(Name)+".csv", "r"), delimiter=",")
x = list(reader)
result = numpy.array(x).astype("float")
result = np.transpose(result) 




x = np.arange(0, np.shape(result)[1], 1)
y = np.arange(0, np.shape(result)[0],  1)
xx, yy = np.meshgrid(x, y)
z = result
f = interpolate.interp2d(x, y, z, kind='cubic')



xnew = np.arange(0, np.shape(result)[1], 5)
ynew = np.arange(0, np.shape(result)[0],  1)
znew = f(xnew, ynew)


print(np.shape(znew)) 


np.savetxt("u_new.csv", np.transpose(znew), delimiter=",")


with open("./v"+str(Name)+".csv", 'r') as file:
  csvreader = csv.reader(file)
  i=0
  for line in csvreader:
    i=i+1
    

import csv
import numpy
reader = csv.reader(open("./v"+str(Name)+".csv", "r"), delimiter=",")
x = list(reader)
result = numpy.array(x).astype("float")
result = np.transpose(result) 




x = np.arange(0, np.shape(result)[1], 1)
y = np.arange(0, np.shape(result)[0],  1)
xx, yy = np.meshgrid(x, y)
z = result
f = interpolate.interp2d(x, y, z, kind='cubic')




xnew = np.arange(0, np.shape(result)[1], 5)
ynew = np.arange(0, np.shape(result)[0],  1)
znew = f(xnew, ynew)


print(np.shape(znew)) 


np.savetxt("v_new.csv", np.transpose(znew), delimiter=",")





