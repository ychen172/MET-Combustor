from DropEvapFun import EvapCalc
import matplotlib.pyplot as plt
import numpy as np
#Parameters
YOxi = "O2:0.2314,N2:0.7622,H2O:0.0064"
YFue = "C2H5OH:1"
ReacMech = 'CRECK2019-NOx-NoBin.yaml'
Pref = 485501 #Pa
Tinf = 487.233 #K
Tdrop = 298 #K Droplet Temperature
YOxiInf = 1.0 #Farfield Composition
rhol = 789 #kg/m3
cpl = 2570 #J/kg/K
LHVapor = 918187.9354880721 #J/kg
LHVheat = 27728322.542190198 #J/kg
FAst = 0.11107613244759687 #Stoichiometric Fuel Air Ratio
TsatRef = 298 #K Fuel saturation temperature
PsatRef = 0.008e6 #Pa 0.008 MPa 
MachLiner = 0.35
Gamma = 1.4
fracMassEvap = 0.1
#Test Conditions
rdIniList = [50e-6,30e-6,10e-6] #m radius of droplet
#Solve
Result = []
for i in range(len(rdIniList)):
    Result.append(EvapCalc(YOxi,YFue,ReacMech,Pref,Tinf,Tdrop,YOxiInf,rdIniList[i],rhol,cpl,LHVapor,LHVheat,FAst,TsatRef,PsatRef,MachLiner,Gamma,fracMassEvap))
"""
rdLst = Result[0]
time = Result[1]
DelTMLst = Result[2]
NuLst = Result[3]
ReLst = Result[4]
PrLst = Result[5]
velRelLst = Result[6]
mdotFLst = Result[7]
rFLst = Result[8]
TFLst = Result[9]
TSLst = Result[10]
YSLst = Result[11]
"""
#Plot Radius Evolution
fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for i in range(len(rdIniList)):
    plt.plot(Result[i][1]*1e3, (Result[i][0]*1e6)**2)
plt.ylabel('Square of Droplet Radius [um2]')
plt.xlabel('Time [ms]')
plt.title('Relative Mach Number '+str(MachLiner))
plt.grid(True)
plt.savefig("rdSquare.jpg")

#Plot T Flame
fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for i in range(len(rdIniList)):
    plt.plot(Result[i][1][:-1]*1e3, Result[i][9][:-1], label="rd: "+str(rdIniList[i]*1e6)+" um")
plt.ylabel('Flame Temperature [K]')
plt.ylim([np.mean(Result[i][9][:-1])-5,np.mean(Result[i][9][:-1])+5])
plt.xlabel('Time [ms]')
plt.title('Relative Mach Number '+str(MachLiner))
plt.grid(True)
plt.savefig("TFlame.jpg")

#Plot Droplet LifeTime
fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for i in range(len(rdIniList)):
    plt.plot(rdIniList[i]*1e6,Result[i][1][-1]*1e3,'+')
plt.ylabel('Droplet Lifetime [ms]')
plt.xlabel('Initial Droplet Radius [um]')
plt.title('Relative Mach Number '+str(MachLiner))
plt.grid(True)
plt.savefig("Droplife.jpg")