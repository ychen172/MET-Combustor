import numpy as np
import cantera as ct
from scipy.optimize import fsolve
from scipy.optimize import least_squares
# Function list
def DM_Fun(gamma,Ma):
    return np.sqrt(gamma)*Ma*(1+0.5*(gamma-1)*(Ma**2))**((1+gamma)/(2*(1-gamma)))
def Pres_Calc(P,Pt,Gam,R,phi,YFuel,YOxid,Tini,Aliner,mdot):
    Mach = np.sqrt((2/(Gam-1))*((Pt/P)**((Gam-1)/Gam) - 1))
    gas = ct.Solution('gri30.yaml')
    gas.set_equivalence_ratio(phi,fuel=YFuel,oxidizer=YOxid,basis='mass')
    gas.TP = Tini,P
    gas.equilibrate('HP')
    TBurnt = gas.T
    RhoBurnt = gas.density_mass
    SOS = np.sqrt(Gam*R*TBurnt)
    VelRHS = SOS*Mach
    VelLHS = mdot/(RhoBurnt*Aliner)
    outputLst = np.array([TBurnt,RhoBurnt,SOS])
    return VelRHS,VelLHS,outputLst
def DM_machcalc(vars,extArgs):
    mdot     = extArgs[0] #kg/s
    R        = extArgs[1] #J/(kg*K)
    Tt       = extArgs[2] #K
    Pt       = extArgs[3] #Pa
    area     = extArgs[4] #m2
    gamma    = extArgs[5] 
    Ma       = vars[0]
    residual = (mdot*np.sqrt(R*Tt))/(Pt*area) - DM_Fun(gamma,Ma)
    return [residual]
def Pres_Balance(vars,extArgs):
    Pt     = extArgs[0] #kg/s
    Gam    = extArgs[1] #J/(kg*K)
    R      = extArgs[2] #K
    phi    = extArgs[3] #Pa
    YFuel  = extArgs[4] #m2
    YOxid  = extArgs[5] 
    Tini   = extArgs[6]
    Aliner = extArgs[7]
    mdot   = extArgs[8]
    P      = vars[0]
    VelRHS,VelLHS,dummy = Pres_Calc(P,Pt,Gam,R,phi,YFuel,YOxid,Tini,Aliner,mdot)
    residual = VelRHS - VelLHS
    return [residual]
def holePattCalc(numHTot,fracH,nRow,charaLen,dj,Dl=[9e9,10e9],isFront=False):
    numH = numHTot * fracH
    numHRow = np.round(numH/nRow)
    numHRow += np.mod(numHRow,2) #Make only even number
    AGap = charaLen/(nRow+1) #Axial Gap[m]
    if isFront:
        Dli=Dl[0]
        Dlo=Dl[1]
        RadH = np.linspace(0.5*Dli+AGap,0.5*Dlo-AGap,nRow) #[m][Array]Radius array of the frontal holes
    else:
        RadH = Dl[0]/2
    TGap = np.sin(((2*np.pi)/numHRow)*0.5)*RadH*2 #Tangential Gap[m][Array]
    AGapRatio = (AGap-dj)/dj
    TGapRatio = (TGap-dj)/dj
    return numHRow,RadH,AGap,AGapRatio,TGapRatio

### User Input
# Design switch
EmbeddedVaporizer = False #If fuel vaporizer is embedded within liner, then inner liner diameter will be corrected
# Combustor inlet gas property
mdot3 = 0.268   #Inlet mass flow [kg/s]
Tt3   = 346.758 #Total temperature [K]
Pt3   = 1.70e5  #Total pressure [Pa]
R3    = 288.512 #Gas constant [J/(kg*K)]
gam3  = 1.397   #Gamma
# Combustor outlet gas property
mdot4 = 0.2708   #Outlet mass flow [kg/s] 
Tt4   = 802.588 #Total temperature [K]
Pt4   = 1.63e5  #Total pressure [Pa]
R4    = 288.469 #Gas constant [J/(kg*K)]
gam4  = 1.344   #Gamma
# Combustor geometry
Dco   = 0.145   #Outer casing diameter [m]
Dci   = 0.070   #Inner casing diameter [m]
Dinleto = 0.14 #Outer inlet diameter [m]
Dinleti = 0.1288 #Inner inlet diameter [m]
Ma4   = 0.1     #Exit Mach number
YFuel = {'C3H8':1}
YOxid = {'O2':0.230,'N2':0.769,'H2O':(1-0.230-0.769)} #Oxidizer must contain only O2 N2 and H2O
phiPz = 1.8 #Primary zone
phiSz = 0.6 #Secondary zone
fracLossPz = 0.40 #Primary zone
fracLossSz = 0.40 #Secondary zone
TauTot = 4e-3 #Total residence time [s]
LszLpzR = 2.2 #length of secondary zone over primary zone 
PHRjpz = 0.112 #Fraction of the liner height for penetration
PHRsz = 0.5 #Fraction of the liner height for penetration
PHRdz = 0.10 #Fraction of the liner height for penetration
velRatJetFuel = 7.4 #Ufuel/UJetAirPz
fracPzFuelPipe = 0.194 #Frac of primary zone air use for vaporization
numPzFuelInj = 10 #Number of primary zone fuel injectors
thiPipeFuel = 0.508/1000 #[m] Thickness of copper and stainless steel pipeline
thiPipeAir = 0.508/1000 #[m] Thickness of copper and stainless steel pipeline
# Reference parameters
Patm  = 101325  #Reference atmosphere pressure [Pa]
Tfuel = 298 #Temperature of fuel before injection K
# Hole Pattern
fracPzFront = 0.60 
fracPzOut = 0.40
nRowPzFront = 4
nRowPzOut = 3
nRowPzIn = 2
fracSzOut = 1.0
nRowSzOut = 3
nRowSzIn = 3
fracHLenSz = 1.0 #fraction of frontal length over which holes are punched (Only for quenchin holes)
fracDzOut = 0.6
nRowDzOut = 14
nRowDzIn = 14
### User Input End
print("Combustor Design Python Code: \n\n\n\n\n\n")
#Print out operating conditions
print("####Print out operating conditions####")
print("======================================")
mdotf = mdot4-mdot3 #Total fuel flow rate [kg/s]
print("Inlet Air Flow: "+str(np.round(mdot3,7))+" kg/s")
print("Inlet Fuel Flow: "+str(np.round(mdotf,7))+" kg/s")
print("Outlet Exhaust Flow: "+str(np.round(mdot4,7))+" kg/s")
print("Inlet Air Total Temperature: "+str(np.round(Tt3,3))+" K")
print("Inlet Fuel Total Temperature: "+str(np.round(Tfuel,3))+" K")
print("Outlet Exhaust Total Temperature: "+str(np.round(Tt4,3))+" K")
print("Inlet Air Total Pressure: "+str(np.round(Pt3,3))+" Pa")
print("Outlet Exhaust Total Pressure: "+str(np.round(Pt4,3))+" Pa")
print("Inlet Air Mass Fraction: "+str(YOxid))
print("Inlet Fuel Mass Fraction: "+str(YFuel))
print(" ")

#Reference Parameters
print("####Compute reference parameters####")
print("====================================")
Aref = (np.pi/4)*(Dco**2 - Dci**2) #Reference area [m2]
Ainlet = (np.pi/4)*(Dinleto**2 - Dinleti**2) #Inlet area [m2]
#Mach Numbers
Maref = fsolve(DM_machcalc,[0.4],args = [mdot3,R3,Tt3,Pt3,Aref,gam3])[0] #Reference Mach number & try initial guess between [0 and 1]
Mainlet = fsolve(DM_machcalc,[0.4],args = [mdot3,R3,Tt3,Pt3,Ainlet,gam3])[0] #Reference Mach number & try initial guess between [0 and 1]
#Reference Dynamic Pressure
rhot3 = Pt3/(R3*Tt3) #Total density [kg/m3]
rho3 = rhot3*(1+0.5*(gam3-1)*(Mainlet**2))**(1/(1-gam3))
Uref = mdot3/(rho3*Aref) #Reference velocity [m/s]
qref = (rho3*(Uref**2))/2 #Reference dynamic pressure [Pa]
dPt34 = Pt3-Pt4 #Stagnation pressure loss across combustor [Pa]
ArefVer2 = np.sqrt((R3/2)*((mdot3*np.sqrt(Tt3)/Pt3)**2)*((dPt34/qref)/(dPt34/Pt3))) #Reference area version 2 [m2]
print("Casing Front Mach Number: "+str(np.round(Maref,3)))
print("Combustor Inlet Mach Number: "+str(np.round(Mainlet,3)))
print("Difference Using Two Area Calculation Method: "+str(np.round(np.abs(Aref-ArefVer2)*100/Aref,3))+" %")
print(" ")

#Liner Area
print("####Compute linear crossectional area####")
print("=========================================")
Aft = (mdot4*np.sqrt(R4*Tt4))/(Pt4*DM_Fun(gam4,Ma4)) #Liner area [m2]
AftVer2 = (1.621e-2)*((mdotf*np.sqrt(Tt3))/(Pt3/Patm))*np.sqrt(Pt3/dPt34) #Liner area version 2 [m2] 
factAft = Aft/Aref #Area factor
factAftVer2 = AftVer2/Aref #Area factor version 2
Dlo = np.sqrt(0.5*((4*Aft)/np.pi + Dco**2 + Dci**2)) #Liner outer diameter [m]
Dli = np.sqrt(Dco**2 + Dci**2 - Dlo**2) #Liner inner diameter [m]
Hl = 0.5*(Dlo-Dli) #Liner height [m]
Hc = 0.5*(Dco-Dci) #Casing height [m]
print("Calculated Area Ratio between Liner and Casing: "+str(np.round(factAft,3)))
print("Good Reference Area Ratio between Liner and Casing: "+str(np.round(factAftVer2,3)))
print("Liner Height: "+str(np.round(Hl*1e3,3))+" mm")
print("Casing Height: "+str(np.round(Hc*1e3,3))+" mm")
print("Liner Inner Diameter: "+str(np.round(Dli*1e3))+" mm")
print("Liner Outer Diameter: "+str(np.round(Dlo*1e3))+" mm")
print("Casing Inner Diameter: "+str(np.round(Dci*1e3))+" mm")
print("Casing Outer Diameter: "+str(np.round(Dco*1e3))+" mm")
print("Combustor Inlet Inner Diameter: "+str(np.round(Dinleti*1e3))+" mm")
print("Combustor Inlet Outer Diameter: "+str(np.round(Dinleto*1e3))+" mm")
print(" ")

#Compute dump pressure loss
print("####Compute pressure loss in dump####")
print("=====================================")
ARexpand = Ainlet/Aref
MNexpand = 1+0.5*(gam3-1)*(Mainlet**2)
PtLossFracDump = ((1-ARexpand)**2)*(1-MNexpand**(gam3/(1-gam3))) #dPt/Ptinlet
print("Dump Expansion Area Ratio: "+str(np.round(1/ARexpand,3)))
print("Before Expansion Mach Number: "+str(np.round(Mainlet,3)))
print("Dump Expansion Percentage Pressure Loss: "+str(np.round(PtLossFracDump*100,3))+" %")
print(" ")

#Inside liner air distribution
print("####Compute air distribution####")
print("================================")
gas = ct.Solution('gri30.yaml')
gas.set_equivalence_ratio(1.0,fuel=YFuel,oxidizer=YOxid,basis='mass')
YOxid_ST = gas['O2'].Y[0]+gas['N2'].Y[0]+gas['H2O'].Y[0]
YFuel_ST = 1-YOxid_ST
FAR_ST = YFuel_ST/YOxid_ST #Stoichiometric fuel air ratio 
mdotaPz = mdotf/(phiPz*FAR_ST) #Primary zone air flow [kg/sec]
mdotaSz = mdotf/(phiSz*FAR_ST) - mdotaPz #Additional air injected throught secondary zone [kg/sec]
mdotaDz = mdot3 - mdotaPz - mdotaSz #Dilution zone air injected [kg/sec]
phi4 = (mdotf/mdot3)/FAR_ST #Total equivalence ratio
fracAirPz = mdotaPz/mdot3 #Fractional primary zone air
fracAirSz = mdotaSz/mdot3 #Fractional secondary zone air
fracAirDz = mdotaDz/mdot3 #Fractional Dilution zone air
print("Stoichiometric Fuel Air Ratio: "+str(np.round(FAR_ST,4)))
print("Percentage Primary Zone Air: "+str(np.round(fracAirPz*100,3))+" %")
print("Percentage Secondary Zone Air: "+str(np.round(fracAirSz*100,3))+" %")
print("Percentage Dilution Zone Air: "+str(np.round(fracAirDz*100,3))+" %")
print("Primary Zone Equivlalence Ratio: "+str(np.round(phiPz,3)))
print("Secondary Zone Equivlalence Ratio: "+str(np.round(phiSz,3)))
print("Overall Combustor Equivalence Ratio: "+str(np.round(phi4,3)))
print(" ")

#Calculate density and dynamic pressure
print("####Compute Dynamic Pressure within each zone####")
print("=================================================")
Ptlo = Pt3 - Pt3*(PtLossFracDump) #[Pa]
dPtliner = Ptlo-Pt4 #[Pa]
Ptpz = Ptlo - dPtliner*(fracLossPz) #[Pa]
Ptsz = Ptpz - dPtliner*(fracLossSz) #[Pa]
Ptdz = Pt4 #[Pa]
Ppz = least_squares(Pres_Balance,[Ptpz*0.99],bounds=(Patm,Ptpz),args = [[Ptpz,gam4,R4,phiPz,YFuel,YOxid,Tt3,Aft,(mdotaPz+mdotf)]]).x[0]
dummy,Upz,outTemp = Pres_Calc(Ppz,Ptpz,gam4,R4,phiPz,YFuel,YOxid,Tt3,Aft,(mdotaPz+mdotf))
print('Primary Zone Pressure Solution Error: '+str(abs(Upz-dummy)))
Tpz = outTemp[0]
rhopz = outTemp[1]
SOSpz = outTemp[2]
Psz = least_squares(Pres_Balance,[Ptsz*0.99],bounds=(Patm,0.5*(Ptpz+Ptsz)),args = [[0.5*(Ptpz+Ptsz),gam4,R4,0.5*(phiPz+phiSz),YFuel,YOxid,Tt3,Aft,(mdotaPz+mdotf+0.5*mdotaSz)]]).x[0]
dummy,Usz,outTemp = Pres_Calc(Psz,0.5*(Ptpz+Ptsz),gam4,R4,0.5*(phiPz+phiSz),YFuel,YOxid,Tt3,Aft,(mdotaPz+mdotf+0.5*mdotaSz))
print('Secondary Zone Pressure Solution Error: '+str(abs(Usz-dummy)))
Tsz = outTemp[0]
rhosz = outTemp[1]
SOSsz = outTemp[2]
Pdz = least_squares(Pres_Balance,[Ptdz*0.99],bounds=(Patm,0.5*(Ptsz+Ptdz)),args = [[0.5*(Ptsz+Ptdz),gam4,R4,0.5*(phiSz+phi4),YFuel,YOxid,Tt3,Aft,(mdotaPz+mdotf+mdotaSz+0.5*mdotaDz)]]).x[0]
dummy,Udz,outTemp = Pres_Calc(Pdz,0.5*(Ptsz+Ptdz),gam4,R4,0.5*(phiSz+phi4),YFuel,YOxid,Tt3,Aft,(mdotaPz+mdotf+mdotaSz+0.5*mdotaDz))
print('Dilution Zone Pressure Solution Error: '+str(abs(Udz-dummy)))
Tdz = outTemp[0]
rhodz = outTemp[1]
SOSdz = outTemp[2]
qpz = 0.5*(rhopz*Upz**2) #Primary zone dynamic pressure [Pa]
qsz = 0.5*(rhosz*Usz**2) #Secondary zone dynamic pressure [Pa]
qdz = 0.5*(rhodz*Udz**2) #Dilution zone dynamic pressure [Pa]
Mapz = Upz/SOSpz #Primary zone mean Mach number
Masz = Usz/SOSsz #Secondary zone mean Mach number
Madz = Udz/SOSdz #Dilution zone mean Mach number
print("Primary zone mean Mach number: "+str(np.round(Mapz,3)))
print("Secondary zone mean Mach number: "+str(np.round(Masz,3)))
print("Dilution zone mean Mach number: "+str(np.round(Madz,3)))
print("Primary zone static temperature: "+str(np.round(Tpz,3))+" K")
print("Secondary zone static temperature: "+str(np.round(Tsz,3))+" K")
print("Dilution zone static temperature: "+str(np.round(Tdz,3))+" K")
print("Outside Liner Stagnation Pressure: "+str(np.round(Ptlo,3))+" Pa")
print("Primary Zone Stagnation Pressure: "+str(np.round(Ptpz,3))+" Pa")
print("Secondary Zone Stagnation Pressure: "+str(np.round(Ptsz,3))+" Pa")
print(" ")

#Calculate jet velocity
print("####Compute Mean Jet Velocity Within Each Zone####")
print("==================================================")
Ujpz = np.sqrt(2*(Ptlo-Ppz)/rho3) #Primary zone jet velocity [m/s]
Ujsz = np.sqrt(2*(Ptlo-Psz)/rho3) #Secondary zone jet velocity [m/s]
Ujdz = np.sqrt(2*(Ptlo-Pdz)/rho3) #Dilution zone jet velocity [m/s]
Majpz = Ujpz/SOSpz #Primary zone cooling jet Mach number
Majsz = Ujsz/SOSsz #Secondary zone cooling jet Mach number
Majdz = Ujdz/SOSdz #Dilution zone cooling jet Mach number
print("Primary Zone Cooling Jet Mach Number: "+str(np.round(Majpz,3)))
print("Secondary Zone Cooling Jet Mach Number: "+str(np.round(Majsz,3)))
print("Dilution Zone Cooling Jet Mach Number: "+str(np.round(Majdz,3)))
print(" ")

#Calculate total effect jet area
print("####Compute the effective jet area based on reference####")
print("=========================================================")
AjPzTEff = mdotaPz/(rho3*Ujpz) #Total effective jet area for primary zone [m2]
AjSzTEff = mdotaSz/(rho3*Ujsz) #Total effective jet area for secondary zone [m2]
AjDzTEff = mdotaDz/(rho3*Ujdz) #Total effective jet area for primary zone [m2]
AjTEffEst = Aref/np.sqrt(dPtliner/qref) #Reference Estimated total jet hole area version [m2]
print('Primary Zone Total Effective Jet Area: '+str(np.round(AjPzTEff,6))+' m2')
print('Secondary Zone Total Effective Jet Area: '+str(np.round(AjSzTEff,6))+' m2')
print('Dilution Zone Total Effective Jet Area: '+str(np.round(AjDzTEff,6))+' m2')
print('Computed All Zone Total Effective Jet Area: '+ str(np.round(AjPzTEff+AjSzTEff+AjDzTEff,6))+' m2')
print('Reference All Zone Total Effective Jet Area: '+ str(np.round(AjTEffEst,6))+' m2')
print(" ")

#Correct total effect jet area(scaled the calculated characteristic jet total effective area such that the sum of area match with estimated total area)
print("####Correct the effective jet area based on reference####")
print("=========================================================")
AeffCorrFact = AjTEffEst/(AjPzTEff+AjSzTEff+AjDzTEff) #Scaling factor Area Estimated over Area Calculated
AjPzTEff = AjPzTEff*AeffCorrFact
AjSzTEff = AjSzTEff*AeffCorrFact
AjDzTEff = AjDzTEff*AeffCorrFact
print("Corrected Total Effective Jet Area: "+str(np.round(AjPzTEff+AjSzTEff+AjDzTEff,6))+" m2")
print(" ")

#Compute Fuel Pipe Area
print("####Compute Fuel Pipe Area####")
print("==============================")
gas.TPY = Tfuel,Ppz,"C3H8:1"
rhoFuel = gas.density_mass #[kg/m3]
UjFuel = Ujpz*velRatJetFuel #[m/s]
PtFuelUp = 0.5*(UjFuel**2)*rhoFuel + Ppz #Upstream pressure needed for JetFuel [Pa] 
ApPzFuelTEff = mdotf/(rhoFuel*UjFuel) #Effective total fuel area for primary zone vaporizer [m2]
ApPzFuelEff = ApPzFuelTEff/numPzFuelInj #Effective fuel area for primary zone vaporizer [m2]
print("Fuel Density: "+str(np.round(rhoFuel,3))+" kg/m3")
print("Fuel Velocity: "+str(np.round(UjFuel,3))+" m/s")
print("Fuel Tank Pressure: "+str(np.round(PtFuelUp,3))+" Pa")
print(" ")

#Compute Air Pipe Area
print("####Compute Air Pipe Area####")
print("=============================")
mdotaPzPipe = mdotaPz*fracPzFuelPipe #primary zone air mass flow used for vaporizer [kg/s]
mdotaPzJet = mdotaPz*(1-fracPzFuelPipe) #primary zone air mass flow used for liner cooling [kg/s]
ApPzAirTEff = AjPzTEff*fracPzFuelPipe #Effective total air area for primary zone vaporizer [m2]
AjPzAirTEff = AjPzTEff*(1-fracPzFuelPipe) #Effective total air area for primary zone cooling [m2]
ApPzAirEff = ApPzAirTEff/numPzFuelInj #Effective air area for primary zone vaporizer [m2]
print(" ")

#Calculate effective jet dimater
print("####Compute Effective Jet Diameter####")
print("======================================")
JGRpz = (0.5*rho3*Ujpz**2)/(qpz) #Dynamic pressure ratio between jet and gas
JGRsz = (0.5*rho3*Ujsz**2)/(qsz) #Dynamic pressure ratio between jet and gas
JGRdz = (0.5*rho3*Ujdz**2)/(qdz) #Dynamic pressure ratio between jet and gas
MRjPz = (mdotaPz+mdotf)/((mdotaPz+mdotf)+mdotaPzJet) #Mass flow ratio for primary zone cooling (mgas/(mgas+mjet))
MRSz = (mdotaPz+mdotf)/((mdotaPz+mdotf)+mdotaSz) #Mass flow ratio for secondary zone (mgas/(mgas+mjet))
MRDz = (mdotaPz+mdotf+mdotaSz)/((mdotaPz+mdotf+mdotaSz)+mdotaDz) #Mass flow ratio (mgas/(mgas+mjet))
djeffjPz = (PHRjpz*Hl)/(1.25*JGRpz**0.5*MRjPz) #Effective air diameter for primary zone cooling [m]
djeffSz = (PHRsz*Hl)/(1.25*JGRsz**0.5*MRSz) #Effective air diameter for secondary zone cooling [m]
djeffDz = (PHRdz*Hl)/(1.25*JGRdz**0.5*MRDz) #Effective air diameter for dilution zone cooling [m]
print("Effective Primary Zone Hole diameter: "+str(np.round(djeffjPz*1e3,3))+" mm")
print("Effective Secondary Zone Hole diameter: "+str(np.round(djeffSz*1e3,3))+" mm")
print("Effective Dilution Zone Hole diameter: "+str(np.round(djeffDz*1e3,3))+" mm")
print(" ")

#Calculate number of holes
print("####Compute Number of holes####")
print("===============================")
numHjpz = AjPzAirTEff/(np.pi*(djeffjPz**2)/4) #Number of holes for primary zone air cooling (can be zero!!!)
numHsz = AjSzTEff/(np.pi*(djeffSz**2)/4) #Number of holes for secondary zone air cooling
numHdz = AjDzTEff/(np.pi*(djeffDz**2)/4) #Number of holes for dilution zone air cooling
print("Number of primary zone holes(Exact&Fractional): "+str(numHjpz))
print("Number of secondary zone holes(Exact&Fractional): "+str(numHsz))
print("Number of dilution zone holes(Exact&Fractional): "+str(numHdz))
print(" ")

#Calculate penetration depth of primary zone fuel air jet
print("####Estimated primary zone fuel-air jet penetration(Pure Air Stream Assumed)####")
print("================================================================================")
MRpPz = (mdotaPz+mdotf)/((mdotaPz+mdotf)+mdotaPzPipe) #Mass flow ratio for primary zone vaporizer (mgas/(mgas+mjet))
dpPzAirEff = np.sqrt(ApPzAirEff*4/np.pi) #Effective air diameter for primary zone vaporizer [m]
PHRppz = (dpPzAirEff*1.25*JGRpz**0.5*MRpPz)/Hl  #Fractional Penetration of primary zone vaporizer stream [fraction of liner height]
print("Primary zone vaporizer jet penetration depth: "+str(np.round(PHRppz*100,3))+" % of liner height")
print(" ")

#Compute discharge coefficients
print("####Compute discharge coefficients####")
print("======================================")
JARpz = (0.5*rho3*Ujpz**2)/qref # Jet annulus dynamic pressure ratio for primary zone
JARsz = (0.5*rho3*Ujsz**2)/qref # Jet annulus dynamic pressure ratio for primary zone
JARdz = (0.5*rho3*Ujdz**2)/qref # Jet annulus dynamic pressure ratio for primary zone
alpha = 0.25 #impericial factor
CdJetpz = (1.25*(JARpz-1))/(4*JARpz**2 - JARpz*(2-alpha)**2)**0.5 # Primary zone discharge coefficient
CdPipepz = (1.65*(JARpz-1))/(4*JARpz**2 - JARpz*(2-alpha)**2)**0.5 # Primary zone discharge coefficient
Cdsz = (1.25*(JARsz-1))/(4*JARsz**2 - JARsz*(2-alpha)**2)**0.5 # Secondary zone discharge coefficient
Cddz = (1.25*(JARdz-1))/(4*JARdz**2 - JARdz*(2-alpha)**2)**0.5 # Dilution zone discharge coefficient
print("Primary zone cooling holes Discharge Coeff: "+str(np.round(CdJetpz,3)))
print("Primary zone vaporizer Discharge Coeff: "+str(np.round(CdPipepz,3)))
print("Secondary zone cooling holes Discharge Coeff: "+str(np.round(Cdsz,3)))
print("Dilution zone cooling holes Discharge Coeff: "+str(np.round(Cddz,3)))
print(" ")

#Compute actual hole size for cooling
print("####Compute actual diameters for cooling####")
print("============================================")
djPz = djeffjPz/(CdJetpz**0.5) # actual air diameter for primary zone cooling [m]
djSz = djeffSz/(Cdsz**0.5) # actual air diameter for secondary zone cooling [m]
djDz = djeffDz/(Cddz**0.5) # actual air diameter for dilution zone cooling [m]
print("Actual Primary Zone Hole diameter: "+str(np.round(djPz*1e3,3))+" mm")
print("Actual Secondary Zone Hole diameter: "+str(np.round(djSz*1e3,3))+" mm")
print("Actual Dilution Zone Hole diameter: "+str(np.round(djDz*1e3,3))+" mm")
print(" ")

#Compute actual hole size for vaporizer
print("####Compute actual diameters for vaporizer####")
print("==============================================")
dpPzFuelEff = np.sqrt(ApPzFuelEff*4/np.pi) #Effective fuel diameter for primary zone vaporizer [m]
djFuelIn = dpPzFuelEff/(CdPipepz**0.5) # actual fuel inner diameter for primary zone vaporizer [m]
djFuelOut = djFuelIn + 2*thiPipeFuel # actual fuel outer diameter for primary zone vaporizer [m]
ApPzMixEff = ApPzAirEff + np.pi*(djFuelOut**2)/4 #Effective mixed area for primary zone vaporizer [m2]
dpPzMixEff = np.sqrt(ApPzMixEff*4/np.pi) #Effevtive mixed diameter for the primary zone vaporizer [m] 
dpPzMixIn = dpPzMixEff/(CdPipepz**0.5) # actual mixed inner diameter for primary zone vaporizer [m]
dpPzMixOut = dpPzMixIn + 2*thiPipeAir # actual mixed outer diameter for primary zone vaporizer [m]
print("Effective Diameter of vaporizer fuel stream: "+str(np.round(dpPzFuelEff*1e3,3))+" mm")
print("Effective Diameter of vaporizer mixture stream: "+str(np.round(dpPzMixEff*1e3,3))+" mm")
print("Actual Inner Diameter of vaporizer fuel pipe: "+str(np.round(djFuelIn*1e3,3))+" mm")
print("Actual Outer Diameter of vaporizer fuel pipe: "+str(np.round(djFuelOut*1e3,3))+" mm")
print("Actual Inner Diameter of vaporizer mixture pipe: "+str(np.round(dpPzMixIn*1e3,3))+" mm")
print("Actual Outer Diameter of vaporizer mixture pipe: "+str(np.round(dpPzMixOut*1e3,3))+" mm")
print("Fuel Pipe Wall Thickness: "+str(np.round(thiPipeFuel*1e3,3))+" mm")
print("Air Pipe Wall Thickness: "+str(np.round(thiPipeAir*1e3,3))+" mm")
print("Number of vaporizers: "+str(numPzFuelInj))
print(" ")

#Calculate Liner Length
print("####Compute liner length####")
print("============================")
LenPz = Hl #Perfect recirculation primary zone length height ratio to calculate parimary zone length [m]
TauPz = LenPz/Upz #Residence time primary zone [s]
LenSz = LenPz*LszLpzR #Scaling ratio to calculate secondary zone length [m]
TauSz = LenSz/Usz #Residence time secondary zone [s]
TauDz = TauTot-TauPz-TauSz #Leftover dilution zone residence time [s]
LenDz = TauDz*Udz #Length of dilution zone length [m]
LenTot = LenPz+LenSz+LenDz #Total length is the sum [m]
print("Total Residence Time: "+str(np.round(TauTot*1e3,3))+" ms")
print("Primary zone length: "+str(np.round(LenPz*1e3,3))+" mm")
print("Secondary zone length: "+str(np.round(LenSz*1e3,3))+" mm")
print("Dilution zone length: "+str(np.round(LenDz*1e3,3))+" mm")
print("Primary zone range: "+str(0)+" mm to "+str(np.round(LenPz*1e3,3))+" mm")
print("Secondary zone range: "+str(np.round(LenPz*1e3,3))+" mm to "+str(np.round((LenPz+LenSz)*1e3,3))+" mm")
print("Dilution zone range: "+str(np.round((LenPz+LenSz)*1e3,3))+" mm to "+str(np.round(LenTot*1e3,3))+" mm")
print(" ")

#Compensate for the additional area from Pipe Injection
if EmbeddedVaporizer:
    print("####Correct for embedded fuel vaporizer####")
    print("===========================================")
    DliNew = np.sqrt(Dlo**2 - numPzFuelInj*(dpPzMixOut**2) - Aft*4/np.pi) #actual inner diameter of the liner[m]
    HlNew = 0.5*(Dlo-DliNew) #Liner height [m]
    print("Corrected Liner Height: "+str(np.round(HlNew*1e3,3))+" mm")
    print("Corrected Liner Inner Diameter: "+str(np.round(DliNew*1e3,3))+" mm")
    print(" ")
else:
    DliNew = Dli
    HlNew = Hl

#Calculate the hole pattern
print("####Compute Hole Pattern####")
print("============================")
if numHjpz>0.9:
    if fracPzFront>0:
        numHpzFrontRow,RadHpzFront,AGapHpzFront,AGapHpzFrontRatio,TGapHpzFrontRatio = holePattCalc(numHjpz,fracPzFront,nRowPzFront,HlNew,djPz,Dl=[DliNew,Dlo],isFront=True)
        print("Primary Zone Front Face Has "+str(int(numHpzFrontRow))+" holes per row &&&& "+str(nRowPzFront)+" rows with diameter "+str(np.round(djPz*1e3,3))+" mm")
        print("Primary Zone Front Rows Radii:"+str(np.round(RadHpzFront*1e3,3))+" mm")
        print("Primary Zone Front Gap Ratio (Radial):"+str(np.round(AGapHpzFrontRatio,3))+" &&&& (Tangential):"+str(np.round(TGapHpzFrontRatio,3)))
    if fracPzOut>0:
        numHpzOutRow,RadHpzOut,AGapHpzOut,AGapHpzOutRatio,TGapHpzOutRatio = holePattCalc(numHjpz,fracPzOut,nRowPzOut,LenPz,djPz,Dl=[Dlo],isFront=False)
        print("Primary Zone Outer Face Has "+str(int(numHpzOutRow))+" holes per row &&&& "+str(nRowPzOut)+" rows with diameter "+str(np.round(djPz*1e3,3))+" mm")
        print("Primary Zone Outer Rows Axial Locations: "+str(np.round(AGapHpzOut*1e3,3))+" mm")
        print("Primary Zone Outer Gap Ratio (Axial):"+str(np.round(AGapHpzOutRatio,3))+" &&&& (Tangential):"+str(np.round(TGapHpzOutRatio,3)))
    if (1-fracPzOut-fracPzFront)>0:
        numHpzInRow,RadHpzIn,AGapHpzIn,AGapHpzInRatio,TGapHpzInRatio = holePattCalc(numHjpz,(1-fracPzOut-fracPzFront),nRowPzIn,LenPz,djPz,Dl=[DliNew],isFront=False)
        print("Primary Zone Inner Face has "+str(int(numHpzInRow))+" holes per row &&&& "+str(nRowPzIn)+" rows with diameter "+str(np.round(djPz*1e3,3))+" mm")
        print("Primary Zone Inner Rows Axial Locations: "+str(np.round(AGapHpzIn*1e3,3))+" mm")
        print("Primary Zone Inner Gap Ratio (Axial):"+str(np.round(AGapHpzInRatio,3))+" &&&& (Tangential):"+str(np.round(TGapHpzInRatio,3)))
    print(" ")
if numHsz>0.9:
    if fracSzOut>0:        
        numHszOutRow,RadHszOut,AGapHszOut,AGapHszOutRatio,TGapHszOutRatio = holePattCalc(numHsz,fracSzOut,nRowSzOut,LenSz*fracHLenSz,djSz,Dl=[Dlo],isFront=False)
        AGapHszOut += LenPz #Add offset [m]
        print("Secondary Zone Outer Face Has "+str(int(numHszOutRow))+" holes per row &&&& "+str(nRowSzOut)+" rows with diameter "+str(np.round(djSz*1e3,3))+" mm")
        print("Secondary Zone Outer Rows Axial Locations: "+str(np.round(AGapHszOut*1e3,3))+" mm")
        print("Secondary Zone Outer Gap Ratio (Axial):"+str(np.round(AGapHszOutRatio,3))+" &&&& (Tangential):"+str(np.round(TGapHszOutRatio,3)))
    if (1-fracSzOut)>0:
        numHszInRow,RadHszIn,AGapHszIn,AGapHszInRatio,TGapHszInRatio = holePattCalc(numHsz,(1-fracSzOut),nRowSzIn,LenSz*fracHLenSz,djSz,Dl=[DliNew],isFront=False)
        AGapHszIn += LenPz #Add offset [m]
        print("Secondary Zone Inner Face has "+str(int(numHszInRow))+" holes per row &&&& "+str(nRowSzIn)+" rows with diameter "+str(np.round(djSz*1e3,3))+" mm")
        print("Secondary Zone Inner Rows Axial Locations: "+str(np.round(AGapHszIn*1e3,3))+" mm")
        print("Secondary Zone Inner Gap Ratio (Axial):"+str(np.round(AGapHszInRatio,3))+" &&&& (Tangential):"+str(np.round(TGapHszInRatio,3)))
    print(" ")
if numHdz>0.9:
    if fracDzOut>0:
        numHdzOutRow,RadHdzOut,AGapHdzOut,AGapHdzOutRatio,TGapHdzOutRatio = holePattCalc(numHdz,fracDzOut,nRowDzOut,LenDz,djDz,Dl=[Dlo],isFront=False)
        AGapHdzOut += LenPz + LenSz #Add offset [m]
        print("Dilution Zone Outer Face Has "+str(int(numHdzOutRow))+" holes per row &&&& "+str(nRowDzOut)+" rows with diameter "+str(np.round(djDz*1e3,3))+" mm")
        print("Dilution Zone Outer Rows Axial Locations: "+str(np.round(AGapHdzOut*1e3,3))+" mm")
        print("Dilution Zone Outer Gap Ratio (Axial):"+str(np.round(AGapHdzOutRatio,3))+" &&&& (Tangential):"+str(np.round(TGapHdzOutRatio,3)))
    if (1-fracDzOut)>0:
        numHdzInRow,RadHdzIn,AGapHdzIn,AGapHdzInRatio,TGapHdzInRatio = holePattCalc(numHdz,(1-fracDzOut),nRowDzIn,LenDz,djDz,Dl=[DliNew],isFront=False)
        AGapHdzIn += LenPz + LenSz #Add offset [m]
        print("Dilution Zone Inner Face has "+str(int(numHdzInRow))+" holes per row &&&& "+str(nRowDzIn)+" rows with diameter "+str(np.round(djDz*1e3,3))+" mm")
        print("Dilution Zone Inner Rows Axial Locations: "+str(np.round(AGapHdzIn*1e3,3))+" mm")
        print("Dilution Zone Inner Gap Ratio (Axial):"+str(np.round(AGapHdzInRatio,3))+" &&&& (Tangential):"+str(np.round(TGapHdzInRatio,3)))
    print(" ")
print('end')
