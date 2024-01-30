import numpy as np
#Usr Input
mdotMeas = np.array([0.015999, 0.04995421, 0.277])# measured flowrates kg/s[Pz Sz Dz]
DiaCur = np.array([6.528,7.000,3.394])# diamters of holes mm
numHoleCur = np.array([12,24,300])# number of holes on each sec
incSZHole = True
numHoleNew = np.array([28,300]) #Sz Dz
#Design Parameter
mdotTarg = np.array([0.0231, 0.0849, 0.276])# target flowrates at three planes kg/s[Pz Sz Dz]
#Get the individual mass injection
mdotMeas = np.concatenate([[mdotMeas[0]],np.diff(mdotMeas)])
mdotTarg = np.concatenate([[mdotTarg[0]],np.diff(mdotTarg)])
#Get Current Total Area
AreaCur = np.pi*((DiaCur/2)**2)*numHoleCur #Curr Total Area [mm2]
#Keep Pz port unchanged
fracSzDzMea = (mdotMeas[1]+mdotMeas[2])/sum(mdotMeas)
fracSzDzTar = (mdotTarg[1]+mdotTarg[2])/sum(mdotTarg)
fracChangeSzDz = (fracSzDzTar-fracSzDzMea)/fracSzDzMea # %Change of SZ DZ combined area required
AreaModSZDZ = AreaCur*1.0
AreaModSZDZ[1] = AreaModSZDZ[1]*(1+fracChangeSzDz)
AreaModSZDZ[2] = AreaModSZDZ[2]*(1+fracChangeSzDz)
#Get SZ DZ mass distribution correct
fracSZMeas = mdotMeas[1]/(mdotMeas[1]+mdotMeas[2])
fracSZTarg = mdotTarg[1]/(mdotTarg[1]+mdotTarg[2])
AreaChangeOneSide = np.abs((fracSZTarg-fracSZMeas)/2)*(AreaModSZDZ[1]+AreaModSZDZ[2])
AreaModFinal = AreaModSZDZ*1.0
if incSZHole:
    AreaModFinal[1] += AreaChangeOneSide
    AreaModFinal[2] -= AreaChangeOneSide
else:
    AreaModFinal[2] += AreaChangeOneSide
    AreaModFinal[1] -= AreaChangeOneSide 
#Compute new hole distribution
DiaNew = DiaCur*1
DiaNew[1] = np.sqrt((AreaModFinal[1]/numHoleNew[0])/np.pi)*2
DiaNew[2] = np.sqrt((AreaModFinal[2]/numHoleNew[1])/np.pi)*2
print("New SZ Dia: "+str(DiaNew[1])+" mm with "+str(numHoleNew[0])+" holes")
print("New DZ Dia: "+str(DiaNew[2])+" mm with "+str(numHoleNew[1])+" holes")






