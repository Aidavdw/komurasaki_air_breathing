import numpy as np
import matplotlib.pyplot as plt
from math import *
from collections import defaultdict
import matplotlib.cm as cm
import sys
import os
import matplotlib.animation as animation
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'





# PARAMETER DICTIONARIES
N_valve=4
p_wall = []
time=[]
Po=1013e2
os.mkdir('Post_treatment_data')


### Pressure at wall figure 


filename='OUTPUT/p_wall.dat'
file = open(filename,"r")
lineIndex =0 


i=0
for line in file.readlines():
	i=i+1
	if i>1:
		A=line.split(" ")
		time.append(float(A[0]))
		p_wall.append((float(A[1])/Po))



fig=plt.figure()
fig.set_size_inches(12.5, 10.5)

plt.rc('text',usetex=True)
plt.plot(time, p_wall)

plt.xlabel('time in $s$',fontSize=20)
plt.ylabel('Pressure at wall in $atm$', fontSize=20)


fig.savefig('Post_treatment_data/P_wall_history.png', dpi=150) 



### Mean flow rate 
MFR = []
time=[]
filename='OUTPUT/mfr_intake.dat'
file = open(filename,"r")
lineIndex =0 

i=0
for line in file.readlines():
	i=i+1
	if i>1:
		A=line.split(" ")
		time.append(float(A[0]))
		MFR.append(float(A[1]))



fig=plt.figure()
fig.set_size_inches(12.5, 10.5)

plt.rc('text',usetex=True)
plt.plot(time, MFR)
plt.xlabel('time in $s$',fontSize=20)
plt.ylabel('Mean flow rate $m^{-3}s^{-1}$', fontSize=20)


fig.savefig('Post_treatment_data/MFR_history.png', dpi=150) 


### Mean density
Rho_mean = []
time_d=[]
filename='OUTPUT/RHO_tube_mean.dat'
file = open(filename,"r")
lineIndex =0 

i=0
for line in file.readlines():
	i=i+1
	if i>1:
		A=line.split(" ")
		time_d.append(float(A[0]))
		Rho_mean.append(float(A[1]))



fig=plt.figure()
fig.set_size_inches(12.5, 10.5)

plt.rc('text',usetex=True)
plt.plot(time_d, Rho_mean)
plt.xlabel('time in $s$',fontSize=20)
plt.ylabel('Mean tube density $kg.m^{-3}$', fontSize=20)


fig.savefig('Post_treatment_data/Rho_tube_mean_history.png', dpi=150) 


### Mean pressure
p_mean = []
time_p=[]
filename='OUTPUT/P_tube_mean.dat'
file = open(filename,"r")
lineIndex =0 

i=0
for line in file.readlines():
	i=i+1
	if i>1:
		A=line.split(" ")
		time_p.append(float(A[0]))
		p_mean.append(float(A[1]))



fig=plt.figure()
fig.set_size_inches(12.5, 10.5)

plt.rc('text',usetex=True)
plt.plot(time_p, p_mean)
plt.xlabel('time in $s$',fontSize=20)
plt.ylabel('Mean tube pressure $(Pa)$', fontSize=20)


### Mean temperature
T_mean = []
time_t=[]
filename='OUTPUT/T_tube_mean.dat'
file = open(filename,"r")
lineIndex =0 

i=0
for line in file.readlines():
	i=i+1
	if i>1:
		A=line.split(" ")
		time_t.append(float(A[0]))
		T_mean.append(float(A[1]))



fig=plt.figure()
fig.set_size_inches(12.5, 10.5)

plt.rc('text',usetex=True)
plt.plot(time_t, T_mean)
plt.xlabel('time in $s$',fontSize=20)
plt.ylabel('Mean tube temperature $(K)$', fontSize=20)

fig.savefig('Post_treatment_data/T_tube_mean_history.png', dpi=150) 



### Pressure / density comparison 

Rho_mean_r=[]
p_mean_r=[]

fig=plt.figure()
fig.set_size_inches(12.5, 10.5)

plt.rc('text',usetex=True)

for k in range(len(Rho_mean)): 
	Rho_mean_r.append(Rho_mean[k]/Rho_mean[0])

plt.plot(time_d, Rho_mean_r, label=r'$\frac{\rho}{\rho_0}$')

for k in range(len(p_mean)): 
	p_mean_r.append(p_mean[k]/p_mean[0])



plt.plot(time_p, p_mean_r, label=r'$\frac{p}{p_0}$')


plt.legend(fontsize=30)
plt.xlabel('time in $s$',fontSize=20)




fig.savefig('Post_treatment_data/Rho_p_tube_mean_ratio_history.png', dpi=150) 







### PFR 

PFR = []
time=[]
filename='OUTPUT/pfr.dat'
file = open(filename,"r")
lineIndex =0 

i=0
for line in file.readlines():
	i=i+1
	if i>1:
		A=line.split(" ")
		time.append(float(A[0]))
		PFR.append(float(A[1]))



fig=plt.figure()
fig.set_size_inches(12.5, 10.5)

plt.rc('text',usetex=True)
plt.plot(time, PFR)
plt.xlabel('time in $s$',fontSize=20)
plt.ylabel('Partial filing rate (no unit)', fontSize=20)


fig.savefig('Post_treatment_data/PFR_history.png', dpi=150) 



### Pressure ratio for each valve 

Pratio = []
time=[]


filename='OUTPUT/p_ratio.dat'
file = open(filename,"r")
lineIndex =0 

i=0
for line in file.readlines():
	i=i+1
	if i>1:
		A=line.split(" ")
		time.append(float(A[0]))
		values=[];
		for k in range(N_valve): 
			values.append(float(A[k+1]))
		Pratio.append(values)



Pratio=np.transpose(Pratio)
print(Pratio)
fig=plt.figure()
fig.set_size_inches(12.5, 10.5)

plt.rc('text',usetex=True)
for k in range(N_valve): 
	plt.plot(time, Pratio[k][:], label='valve number'+str(k))




plt.legend()
plt.xlabel('time in $s$',fontSize=20)
plt.ylabel('Pressure ratio at between reed valve', fontSize=20)


fig.savefig('Post_treatment_data/Pressure_at_valves_history.png', dpi=150) 


### Displacement for each valve 

Disp = []
time=[]


filename='OUTPUT/y_tip.dat'
file = open(filename,"r")
lineIndex =0 

i=0
for line in file.readlines():
	i=i+1
	if i>1:
		A=line.split(" ")
		time.append(float(A[0]))
		values=[];
		for k in range(N_valve): 
			values.append(float(A[k+1]))
		Disp.append(values)



Disp=np.transpose(Disp)
print(Disp)
fig=plt.figure()
fig.set_size_inches(12.5, 10.5)

plt.rc('text',usetex=True)
for k in range(N_valve): 
	plt.plot(time, Disp[k][:], label='valve number'+str(k))

print(max(Disp[0][:]))


plt.legend()
plt.xlabel('time in $s$',fontSize=20)
plt.ylabel('Displacement of the extreme part of each valves (in $m$)', fontSize=20)


fig.savefig('Post_treatment_data/Displacement_at_valves_history.png', dpi=150) 










