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
PHRpz = 0.112 #Fraction of the liner height for penetration
PHRsz = 0.5 #Fraction of the liner height for penetration
PHRdz = 0.10 #Fraction of the liner height for penetration
Ent2JetRat = 1.0 #Entrained peneration to Jet Peneration Ratio 
velRatJetFuel = 7.4 #Ufuel/UJetAirPz
fracPzFuelPipe = 0.194 #Frac of primary zone air use for vaporization
numPzFuelInj = 10 #Number of primary zone fuel injectors
thiPipeFuel = 0.508/1000 #[m] Thickness of copper and stainless steel pipeline
thiPipeAir = 0.508/1000 #[m] Thickness of copper and stainless steel pipeline
# Reference parameters
Patm  = 101325  #Reference atmosphere pressure [Pa]
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

#Reference Area
Aref = (np.pi/4)*(Dco**2 - Dci**2) #Reference area [m2]
Ainlet = (np.pi/4)*(Dinleto**2 - Dinleti**2) #Inlet area [m2]
rhot3 = Pt3/(R3*Tt3) #Total density [kg/m3]
Uref = mdot3/(rhot3*Aref) #Reference velocity [m/s]
qref = (rhot3*(Uref**2))/2 #Reference dynamic pressure [Pa]
dPt34 = Pt3-Pt4 #Stagnation pressure loss across combustor [Pa]
ArefVer2 = np.sqrt((R3/2)*((mdot3*np.sqrt(Tt3)/Pt3)**2)*((dPt34/qref)/(dPt34/Pt3))) #Reference area version 2 [m2]
print('Aref: '+str(Aref)+' m2 while ArefVer2: '+str(ArefVer2)+' m2')
print('Uref: '+str(Uref)+' m/s and qref: '+str(qref)+' Pa')
#Liner Area
Maref = fsolve(DM_machcalc,[0.4],args = [mdot3,R3,Tt3,Pt3,Aref,gam3])[0] #Reference Mach number & try initial guess between [0 and 1]
print('Maref: '+str(Maref))
Mainlet = fsolve(DM_machcalc,[0.4],args = [mdot3,R3,Tt3,Pt3,Ainlet,gam3])[0] #Reference Mach number & try initial guess between [0 and 1]
print('Mainlet: '+str(Mainlet))
Aft = (mdot4*np.sqrt(R4*Tt4))/(Pt4*DM_Fun(gam4,Ma4)) #Liner area [m2]
mdotf = mdot4-mdot3 #Total fuel flow rate [kg/s]
AftVer2 = (1.621e-2)*((mdotf*np.sqrt(Tt3))/(Pt3/Patm))*np.sqrt(Pt3/dPt34) #Liner area version 2 [m2] 
factAft = Aft/Aref #Area factor
factAftVer2 = AftVer2/Aref #Area factor version 2 
print('Aft: '+str(Aft)+' m2 while AftVer2: '+str(AftVer2)+' m2')
print('factAft: '+str(factAft)+' while factAftVer2: '+str(factAftVer2)+', But good factor around 0.66')
Dlo = np.sqrt(0.5*((4*Aft)/np.pi + Dco**2 + Dci**2)) #Liner outer diameter [m]
Dli = np.sqrt(Dco**2 + Dci**2 - Dlo**2) #Liner inner diameter [m]
print('Dco: '+str(Dco)+' m Dci: '+str(Dci)+' m')
print('Dlo: '+str(Dlo)+' m Dli: '+str(Dli)+' m')
Hl = 0.5*(Dlo-Dli) #Liner height [m]
Hc = 0.5*(Dco-Dci) #Casing height [m]
print('Hc: '+str(Hc)+' m Hl: '+str(Hl)+' m')

#Compute dump pressure loss
ARexpand = Ainlet/Aref
MNexpand = 1+0.5*(gam3-1)*(Mainlet**2)
PtLossFracDump = ((1-ARexpand)**2)*(1-MNexpand**(gam3/(1-gam3))) #dPt/Ptinlet
print("ARexpand: "+str(ARexpand) + " PtLossFracDump: "+str(PtLossFracDump))

#Inside liner air distribution
gas = ct.Solution('gri30.yaml')
gas.set_equivalence_ratio(1.0,fuel=YFuel,oxidizer=YOxid,basis='mass')
YOxid_ST = gas['O2'].Y[0]+gas['N2'].Y[0]+gas['H2O'].Y[0]
YFuel_ST = 1-YOxid_ST
FAR_ST = YFuel_ST/YOxid_ST #Stoichiometric fuel air ratio 
mdotaPz = mdotf/(phiPz*FAR_ST) #Primary zone air flow [kg/sec]
mdotaSz = mdotf/(phiSz*FAR_ST) - mdotaPz #Additional air injected throught secondary zone [kg/sec]
mdotaDz = mdot3 - mdotaPz - mdotaSz #Dilution zone air injected [kg/sec]
phi4 = (mdotf/mdot3)/FAR_ST #Total equivalence ratio
print('mdotaPz: '+str(mdotaPz)+' kg/s mdotaSz: '+str(mdotaSz)+' kg/s mdotaDz: '+str(mdotaDz)+' kg/s')
print('phi4: '+str(phi4))

#Compute air distribution
print('perAirPz: '+str(mdotaPz*100/mdot3) + ' % perAirSz: '+str(mdotaSz*100/mdot3) +' % perAirDz: '+str(mdotaDz*100/mdot3))

#Calculate density and dynamic pressure
Ptlo = Pt3 - Pt3*(PtLossFracDump) #[Pa]
dPtliner = Ptlo-Pt4 #[Pa]
Ptpz = Ptlo - dPtliner*(fracLossPz) #[Pa]
Ptsz = Ptpz - dPtliner*(fracLossSz) #[Pa]
Ptdz = Pt4 #[Pa]
Ppz = least_squares(Pres_Balance,[Ptpz*0.99],bounds=(Patm,Ptpz),args = [[Ptpz,gam3,R3,phiPz,YFuel,YOxid,Tt3,Aft,(mdotaPz+mdotf)]]).x[0]
dummy,Upz,outTemp = Pres_Calc(Ppz,Ptpz,gam3,R3,phiPz,YFuel,YOxid,Tt3,Aft,(mdotaPz+mdotf))
print('Residual: '+str(abs(Upz-dummy)))
Tpz = outTemp[0]
rhopz = outTemp[1]
SOSpz = outTemp[2]
print('Upz: '+str(Upz)+' m/s Tpz: '+str(Tpz)+' K rhopz: '+str(rhopz)+' kg/m3 SOSpz: '+str(SOSpz)+' m/s')
Psz = least_squares(Pres_Balance,[Ptsz*0.99],bounds=(Patm,Ptsz),args = [[Ptsz,0.5*(gam3+gam4),0.5*(R3+R4),phiSz,YFuel,YOxid,Tt3,Aft,(mdotaPz+mdotf+mdotaSz)]]).x[0]
dummy,Usz,outTemp = Pres_Calc(Psz,Ptsz,0.5*(gam3+gam4),0.5*(R3+R4),phiSz,YFuel,YOxid,Tt3,Aft,(mdotaPz+mdotf+mdotaSz))
print('Residual: '+str(abs(Usz-dummy)))
Tsz = outTemp[0]
rhosz = outTemp[1]
SOSsz = outTemp[2]
print('Usz: '+str(Usz)+' m/s Tsz: '+str(Tsz)+' K rhosz: '+str(rhosz)+' kg/m3 SOSsz: '+str(SOSsz)+' m/s')
Pdz = least_squares(Pres_Balance,[Ptdz*0.99],bounds=(Patm,Ptdz),args = [[Ptdz,gam4,R4,phi4,YFuel,YOxid,Tt3,Aft,(mdotaPz+mdotf+mdotaSz+mdotaDz)]]).x[0]
dummy,Udz,outTemp = Pres_Calc(Pdz,Ptdz,gam4,R4,phi4,YFuel,YOxid,Tt3,Aft,(mdotaPz+mdotf+mdotaSz+mdotaDz))
print('Residual: '+str(abs(Udz-dummy)))
Tdz = outTemp[0]
rhodz = outTemp[1]
SOSdz = outTemp[2]
print('Udz: '+str(Udz)+' m/s Tdz: '+str(Tdz)+' K rhodz: '+str(rhodz)+' kg/m3 SOSdz: '+str(SOSdz)+' m/s')
qpz = 0.5*(rhopz*Upz**2) #Primary zone dynamic pressure [Pa]
qsz = 0.5*(rhopz*Upz**2) #Secondary zone dynamic pressure [Pa]
qdz = 0.5*(rhosz*Usz**2) #Dilution zone dynamic pressure [Pa]

#Calculate jet velocity
Ujpz = np.sqrt(2*(Ptlo-Ppz)/rhot3) #Primary zone jet velocity [m/s]
Ujsz = np.sqrt(2*(Ptlo-Ppz)/rhot3) #Secondary zone jet velocity [m/s]
Ujdz = np.sqrt(2*(Ptlo-Psz)/rhot3) #Dilution zone jet velocity [m/s]
print('Ujpz: '+str(Ujpz)+' m/s Ujsz: '+str(Ujsz)+' m/s Ujdz: '+str(Ujdz)+' m/s')

#Calculate total effect jet area
AjPzTEff = mdotaPz/(rhot3*Ujpz) #Total effective jet area for primary zone [m2]
AjSzTEff = mdotaSz/(rhot3*Ujsz) #Total effective jet area for secondary zone [m2]
AjDzTEff = mdotaDz/(rhot3*Ujdz) #Total effective jet area for primary zone [m2]
AjPzTEffVer2 = Aref/np.sqrt((Ptlo-Ppz)/qref) #Second version [m2]
AjSzTEffVer2 = Aref/np.sqrt((Ptlo-Ppz)/qref) #Second version [m2]
AjDzTEffVer2 = Aref/np.sqrt((Ptlo-Psz)/qref) #Second version [m2] (Maybe this formula is only applicable to Dz)
print('AjPzTEff: '+str(AjPzTEff)+' m2 AjPzTEffVer2: '+str(AjPzTEffVer2)+' m2')
print('AjSzTEff: '+str(AjSzTEff)+' m2 AjSzTEffVer2: '+str(AjSzTEffVer2)+' m2')
print('AjDzTEff: '+str(AjDzTEff)+' m2 AjDzTEffVer2: '+str(AjDzTEffVer2)+' m2')

#Compute Fuel Pipe Area
gas.TPY = Tpz,Ppz,"C3H8:1"
rhoFuel = gas.density_mass #[kg/m3]
UjFuel = Ujpz*velRatJetFuel #[m/s]
PtFuelUp = 0.5*(UjFuel**2)*rhoFuel + Ppz #[Pa] Upstream pressure needed for JetFuel
AjFuelTEff = mdotf/(rhoFuel*UjFuel) #[m2] Total area for the JetFuel
print("rhoFuel: "+str(rhoFuel)+" kg/m3 UjFuel: "+str(UjFuel)+" m/s PtFuelUp: "+str(PtFuelUp)+" Pa AjFuelTEff: "+str(AjFuelTEff)+" m2")

#Compute Air Pipe Area
mdotaPzPipe = mdotaPz*fracPzFuelPipe #Used for vaporizer
mdotaPzJet = mdotaPz*(1-fracPzFuelPipe) #Used for liner cooling
APipeAirTEff = mdotaPzPipe/(rhot3*Ujpz) #Effective total pipe air area
APipeMixTEff = AjFuelTEff+APipeAirTEff #Effective total pipe mixture area
AJetPzAirTEff = AjPzTEff-APipeAirTEff #Effective total jet pz air area
AjFuelEff = AjFuelTEff/numPzFuelInj 
APipeMixEff = APipeMixTEff/numPzFuelInj
djFuelEff = np.sqrt(AjFuelEff*4/np.pi)
dPipeMixEff = np.sqrt(APipeMixEff*4/np.pi)
print("djFuelEff: "+str(djFuelEff)+" m dPipeMixEff: "+str(dPipeMixEff)+" m")

#Calculate effective jet dimater
JGRpz = (0.5*rhot3*Ujpz**2)/(qpz) #Dynamic pressure ratio between jet and gas
JGRsz = (0.5*rhot3*Ujsz**2)/(qsz) #Dynamic pressure ratio between jet and gas
JGRdz = (0.5*rhot3*Ujdz**2)/(qdz) #Dynamic pressure ratio between jet and gas
MRPz = (mdotaPzPipe+mdotf)/(mdotaPzPipe+mdotf+mdotaPzJet)
MRSz = (mdotaPzPipe+mdotf+mdotaPzJet)/(mdotaPzPipe+mdotf+mdotaPzJet+mdotaSz)
MRDz = (mdotaPzPipe+mdotf+mdotaPzJet+mdotaSz)/(mdotaPzPipe+mdotf+mdotaPzJet+mdotaSz+mdotaDz)
djeffPz = (PHRpz*Hl/Ent2JetRat)/(1.25*JGRpz**0.5*MRPz) #Primary zone jet diameter [m]
djeffSz = (PHRsz*Hl/Ent2JetRat)/(1.25*JGRsz**0.5*MRSz) #Primary zone jet diameter [m]
djeffDz = (PHRdz*Hl/Ent2JetRat)/(1.25*JGRdz**0.5*MRDz) #Primary zone jet diameter [m]
print('djeffPz: '+str(djeffPz)+' m djeffSz: '+str(djeffSz)+' m djeffDz: '+str(djeffDz)+' m')

#Calculate number of holes
numHpz = AJetPzAirTEff/(0.25*np.pi*djeffPz**2)
numHsz = AjSzTEff/(0.25*np.pi*djeffSz**2)
numHdz = AjDzTEff/(0.25*np.pi*djeffDz**2)
print('numHpz: '+str(numHpz)+' numHsz: '+str(numHsz)+' numHdz: '+str(numHdz))

#Compute actual hole size
JARpz = (0.5*rhot3*Ujpz**2)/qref # Jet annulus dynamic pressure ratio for primary zone
JARsz = (0.5*rhot3*Ujsz**2)/qref # Jet annulus dynamic pressure ratio for primary zone
JARdz = (0.5*rhot3*Ujdz**2)/qref # Jet annulus dynamic pressure ratio for primary zone
alpha = 0.25 #impericial factor
CdJetpz = (1.25*(JARpz-1))/(4*JARpz**2 - JARpz*(2-alpha)**2)**0.5 # Primary zone discharge coefficient
CdPipepz = (1.65*(JARpz-1))/(4*JARpz**2 - JARpz*(2-alpha)**2)**0.5 # Primary zone discharge coefficient
Cdsz = (1.25*(JARsz-1))/(4*JARsz**2 - JARsz*(2-alpha)**2)**0.5 # Secondary zone discharge coefficient
Cddz = (1.25*(JARdz-1))/(4*JARdz**2 - JARdz*(2-alpha)**2)**0.5 # Dilution zone discharge coefficient
print('CdJetpz: '+str(CdJetpz)+' CdPipepz: '+str(CdPipepz)+' Cdsz: '+str(Cdsz)+' Cddz: '+str(Cddz))

#Compute actual hole size
djPz = djeffPz/(CdJetpz**0.5) # actual hole diameter [m]
djSz = djeffSz/(Cdsz**0.5) # actual hole diameter [m]
djDz = djeffDz/(Cddz**0.5) # actual hole diameter [m]
dPipeMix = dPipeMixEff/(CdPipepz**0.5) # actual pipe diameter [m]
djFuel = djFuelEff/(CdPipepz**0.5) # actual fuel pipe diameter [m]
dPipeMix = dPipeMix + 2*thiPipeFuel # Take into account the pipe thickness of fuel injector
print('djPz: '+str(djPz)+' m djSz: '+str(djSz)+' m djDz: '+str(djDz)+' m dPipeMix: '+str(dPipeMix) + ' m djFuel' + str(djFuel) + ' m')
dPipeMixOut = dPipeMix + 2*thiPipeAir # Outer diameter [m] Take into account the pipe thickess of the vaporizer
print('dPipeMixOut: '+str(dPipeMixOut)+" m")

#Calculate Liner Length
LenPz = Hl #Due to primary recirculation zone [m]
TauPz = LenPz/Upz #Residence time primary zone [s]
LenSz = LenPz*LszLpzR #Secondary zone [m]
TauSz = LenSz/Usz #Secondary zone [s]
TauDz = TauTot-TauPz-TauSz #Dilution zone [s]
LenDz = TauDz*Udz #Dilution zone [m]
LenTot = LenPz+LenSz+LenDz #Total length [m]
print('TauPz: '+str(TauPz)+' s TauSz: '+str(TauSz)+' s TauDz: '+str(TauDz)+' s TauTot: '+str(TauTot)+' s')
print('LenPz: '+str(LenPz)+' m LenSz: '+str(LenSz)+' m LenDz: '+str(LenDz)+' m LenTot: '+str(LenTot)+' m')

#Compensate for the additional area from Pipe Injection
DliNew = np.sqrt(Dlo**2 - numPzFuelInj*(dPipeMixOut**2) - Aft*4/np.pi) #actual inner diameter of the liner[m]
print("Dli: "+str(Dli)+ " m DliNew: "+str(DliNew)+ " m")

#Calculate the hole pattern
print("\n\n\n#####\n#####\nCalculate Hole Pattern")
numHpzFrontRow,RadHpzFront,AGapHpzFront,AGapHpzFrontRatio,TGapHpzFrontRatio = holePattCalc(numHpz,fracPzFront,nRowPzFront,Hl,djPz,Dl=[Dli,Dlo],isFront=True)
numHpzOutRow,RadHpzOut,AGapHpzOut,AGapHpzOutRatio,TGapHpzOutRatio = holePattCalc(numHpz,fracPzOut,nRowPzOut,LenPz,djPz,Dl=[Dlo],isFront=False)
numHpzInRow,RadHpzIn,AGapHpzIn,AGapHpzInRatio,TGapHpzInRatio = holePattCalc(numHpz,(1-fracPzOut-fracPzFront),nRowPzIn,LenPz,djPz,Dl=[Dli],isFront=False)
print("Front pz "+str(numHpzFrontRow)+" holes times "+str(nRowPzFront)+" rows with diameter "+str(djPz)+" m")
print("Out pz "+str(numHpzOutRow)+" holes times "+str(nRowPzOut)+" rows with diameter "+str(djPz)+" m")
print("In pz "+str(numHpzInRow)+" holes times "+str(nRowPzIn)+" rows with diameter "+str(djPz)+" m")
print("Radius Front pz:"+str(RadHpzFront)+" m")
print("Axial gap Out pz: "+str(AGapHpzOut)+" m In pz: "+str(AGapHpzIn)+" m")
print("Gap Ratio Front pz (A/T):"+str(AGapHpzFrontRatio)+" & "+str(TGapHpzFrontRatio))
print("Gap Ratio Out pz (A/T):"+str(AGapHpzOutRatio)+" & "+str(TGapHpzOutRatio))
print("Gap Ratio In pz (A/T):"+str(AGapHpzInRatio)+" & "+str(TGapHpzInRatio))
print("\n")
numHszOutRow,RadHszOut,AGapHszOut,AGapHszOutRatio,TGapHszOutRatio = holePattCalc(numHsz,fracSzOut,nRowSzOut,LenSz*fracHLenSz,djSz,Dl=[Dlo],isFront=False)
numHszInRow,RadHszIn,AGapHszIn,AGapHszInRatio,TGapHszInRatio = holePattCalc(numHsz,(1-fracSzOut),nRowSzIn,LenSz*fracHLenSz,djSz,Dl=[Dli],isFront=False)
print("Out sz "+str(numHszOutRow)+" holes times "+str(nRowSzOut)+" rows with diameter "+str(djSz)+" m")
print("In sz "+str(numHszInRow)+" holes times "+str(nRowSzIn)+" rows with diameter "+str(djSz)+" m")
print("Axial gap Out sz: "+str(AGapHszOut)+" m In sz: "+str(AGapHszIn)+" m")
print("Gap Ratio Out sz (A/T):"+str(AGapHszOutRatio)+" & "+str(TGapHszOutRatio))
print("Gap Ratio In sz (A/T):"+str(AGapHszInRatio)+" & "+str(TGapHszInRatio)) #-0.5axial ratio will be the target if only one row
print("\n")
numHdzOutRow,RadHdzOut,AGapHdzOut,AGapHdzOutRatio,TGapHdzOutRatio = holePattCalc(numHdz,fracDzOut,nRowDzOut,LenDz,djDz,Dl=[Dlo],isFront=False)
numHdzInRow,RadHdzIn,AGapHdzIn,AGapHdzInRatio,TGapHdzInRatio = holePattCalc(numHdz,(1-fracDzOut),nRowDzIn,LenDz,djDz,Dl=[Dli],isFront=False)
print("Out dz "+str(numHdzOutRow)+" holes times "+str(nRowDzOut)+" rows with diameter "+str(djDz)+" m")
print("In dz "+str(numHdzInRow)+" holes times "+str(nRowDzIn)+" rows with diameter "+str(djDz)+" m")
print("Axial gap Out dz: "+str(AGapHdzOut)+" m In dz: "+str(AGapHdzIn)+" m")
print("Gap Ratio Out dz (A/T):"+str(AGapHdzOutRatio)+" & "+str(TGapHdzOutRatio))
print("Gap Ratio In dz (A/T):"+str(AGapHdzInRatio)+" & "+str(TGapHdzInRatio))
print('end')
