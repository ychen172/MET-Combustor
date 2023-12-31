from DropEvapFun import EvapCalc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#Parameters
YOxi = "O2:0.2314,N2:0.7622,H2O:0.0064"
YFue = "NC12H26:1"#"C2H5OH:1"
ReacMech = 'CRECK2019-NOx-NoBin.yaml'
Tdrop = 298 #K Droplet Temperature
YOxiInf = 1.0 #Farfield Composition
rhol = 749 #kg/m3 nDodecane: 749 kg/m3  Ethanol: 789 kg/m3
cpl = 2211.764 #J/kg/K nDodecane: 2211.764 J/kg/K Ethanol: 2570 J/kg/K
LHVapor = 256000 #J/kg nDodecane: 256000 #J/kg  Ethanol: 918187.9354880721 J/kg
LHVheat = 44462059.86766551 #J/kg nDodecane: 44462059.86766551 Ethanol: 27728322.542190198 J/kg
FAst = 0.06660059875330497 #Stoichiometric Fuel Air Ratio nDodecane: 0.06660059875330497  Ethanol: 0.11107613244759687
TsatRef = 450 #K Fuel saturation temperature nDodecane: 450K  Ethanol: 298K
PsatRef = 35750 #Pa nDodecane: 35750 Pa  Ethanol: 0.008e6 Pa
MachLiner = 0.35
Gamma = 1.4
fracMassEvap = 0.1
#Read Cycle Data
CycleName = "ACSRQLDerived_DualFine_JetAJetAOutPutTot.csv"
Cycle = pd.read_csv(CycleName)
TCycleLst = Cycle['Tt3[K]'].values[:]
PCycleLst = Cycle['Pt3[Pa]'].values[:]
PrefLst = PCycleLst[-1:] #485501 #Pa
TinfLst = TCycleLst[-1:] #487.233 #K
#Test Conditions
rdIniList = [30e-6,25e-6,20e-6,15e-6,10e-6,5e-6] #m radius of droplet
#Solve
Result = []
for i in range(len(rdIniList)):
    Result.append([])
    for j in range(len(TinfLst)):
        Result[i].append(EvapCalc(YOxi,YFue,ReacMech,PrefLst[j],TinfLst[j],Tdrop,YOxiInf,rdIniList[i],rhol,cpl,LHVapor,LHVheat,FAst,TsatRef,PsatRef,MachLiner,Gamma,fracMassEvap))
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
    plt.plot(Result[i][0][1]*1e3, (Result[i][0][0]*1e6)**2)
plt.ylabel('Square of Droplet Radius [um2]')
plt.xlabel('Time [ms]')
plt.title('Relative Mach Number '+str(MachLiner)+' Tinf: '+str(np.round(TinfLst[0],3))+' K Pinf: '+str(np.round(PrefLst[0],3))+' Pa')
plt.grid(True)
plt.savefig("rdSquare.jpg")

#Plot T Flame
fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for i in range(len(rdIniList)):
    plt.plot(Result[i][0][1][:-1]*1e3, Result[i][0][9][:-1], label="rd: "+str(rdIniList[i]*1e6)+" um")
plt.ylabel('Flame Temperature [K]')
plt.xlabel('Time [ms]')
plt.legend()
plt.title('Relative Mach Number '+str(MachLiner)+' Tinf: '+str(np.round(TinfLst[0],3))+' K Pinf: '+str(np.round(PrefLst[0],3))+' Pa')
plt.grid(True)
plt.savefig("TFlame.jpg")

#Single Droplet Lifetime
fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for i in range(len(rdIniList)):
    plt.plot(rdIniList[i]*1e6, Result[i][0][1][-1]*1e3, 'k.')
plt.ylabel('Droplet Lifetime [ms]')
plt.xlabel('Droplet Radius [um]')
plt.title('Relative Mach Number '+str(MachLiner)+' Tinf: '+str(np.round(TinfLst[0],3))+' K Pinf: '+str(np.round(PrefLst[0],3))+' Pa')
plt.grid(True)
plt.savefig("SingleDroplife.jpg")

#Combined Droplet Lifetime
DropLifeTime = []
for i in range(len(rdIniList)):
    DropLifeTime.append([])
    for j in range(len(TinfLst)):
        DropLifeTime[i].append(Result[i][j][1][-1])
fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for i in range(len(rdIniList)):
    plt.plot(TinfLst, np.array(DropLifeTime[i])*1e3,label="rd: "+str(rdIniList[i]*1e6)+" um")
plt.ylabel('Droplet Lifetime [ms]')
plt.xlabel('Combustor Inlet Temperature [K]')
plt.legend()
plt.title('Relative Mach Number '+str(MachLiner))
plt.grid(True)
plt.savefig("Droplife.jpg")