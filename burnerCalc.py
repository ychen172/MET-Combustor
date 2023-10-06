import numpy as np
from scipy.optimize import fsolve
# This is a microturbojet engine combustor design
mdot_at =  0.268 #kg/s
mdot_exit = 0.270 #kg/s
Tt3 = 346.758 #K
Tt4 = 802.588 #K
Pt3 = 1.70e5 #Pa
Pt4 = 1.63e5 #Pa
R_3 = 288.512 #J/(kg*K) Universal gas constant for specific gas
R_4 = 288.469 #J/(kg*K)
Dia_case = 0.178 #m
Dia_shaft = 0.130 #0.025 #m
gamma3 = 1.397
gamma4 = 1.344
Ma_linerDes = 0.1
tau_Des = 4.e-3 #sec
fracPZLength = 0.231 #0.231 to 0.276
fracSZLength = 0.174 #0.174 to 0.207
FA_st = 0.0638 #Stoichiometric Fuel Air Ratio by mass
phi_PZ = 1.11
phi_SZ = 0.59
fracPZAirAxi = 0.6 #60% enter from axially direction
alpha1 = 1.25 #discharge factor for surface hole
YmaxFrac = 0.7 #Fraction of liner height for jet penetration
ProductToAirRatio = 0.5 #1 for full product and 0 for pure air 
Cd_overwrite = True
Cd_overwriteVal = 0.6
#### Reference Area
A_case = np.pi*(Dia_case/2)**2 #m2
A_shaft = np.pi*(Dia_shaft/2)**2 #m2
A_ref = A_case - A_shaft #m2
## Compute mach number from area
def DM_Fun(gamma,Ma):
    return np.sqrt(gamma)*Ma*(1+0.5*(gamma-1)*(Ma**2))**((1+gamma)/(2*(1-gamma)))

def DM_machcalc(vars,extArgs):
    mdot = extArgs[0] #kg/s
    R = extArgs[1] #J/(kg*K)
    Tt = extArgs[2] #K
    Pt = extArgs[3] #Pa
    area = extArgs[4] #m2
    gamma = extArgs[5] 
    Ma = vars[0]
    residual = (mdot*np.sqrt(R*Tt))/(Pt*area) - DM_Fun(gamma,Ma)
    return [residual]
Ma_case = fsolve(DM_machcalc,[0.4],args = [mdot_at,R_3,Tt3,Pt3,A_ref,gamma3])[0] #If crash then try different initial mach guess betwee [0 and 1]
# print('mach number residual is: '+str((mdot_at*np.sqrt(R_3*Tt3))/(Pt3*A_ref) - np.sqrt(gamma3)*Ma_case*(1+0.5*(gamma3-1)*(Ma_case**2))**((1+gamma3)/(2*(1-gamma3)))))
print("Casing Mach :" + str(np.round(Ma_case,3))+" Liner Mach:"+str(np.round(Ma_linerDes,3)))
## continue
rhot3 = Pt3/(R_3*Tt3) #kg/m3
rho3 = rhot3*(1+0.5*(gamma3-1)*(Ma_case**2))**(1/(1-gamma3)) # kg/m3
U_ref = mdot_at/(rho3*A_ref) # m/s
## Verify Mach Number
T3 = Tt3*((1+0.5*(gamma3-1)*(Ma_case**2))**(-1)) #K
SOS3 = np.sqrt(gamma3*R_3*T3) #m/s
Ma_caseV2 = U_ref/SOS3
## Reference Dynamic Head
q_ref = 0.5*rho3*(U_ref**2) #Pa

#### Calculate Liner Area
A_liner = (mdot_exit*np.sqrt(R_4*Tt4))/(Pt4*DM_Fun(gamma4,Ma_linerDes)) #m2
areaLinerFactor = A_liner/A_ref 
print('Area factor is: '+str(np.round(areaLinerFactor,3))+' Recommended Factor is 0.666')
Dia_ave = 0.5*(Dia_case + Dia_shaft)
Dia_linIn = Dia_ave - (A_liner/(Dia_ave*np.pi)) #m
Dia_linOu = 2*Dia_ave - Dia_linIn #m
# print("Verify Area: "+str(np.pi*(Dia_linOu/2)**2 - np.pi*(Dia_linIn/2)**2)+'&'+str(A_liner))
Hei_case = (Dia_case-Dia_shaft)/2 #m
Hei_liner = (Dia_linOu-Dia_linIn)/2 #m
print('Casing: Height:'+str(np.round(Hei_case*1000,3))+'mm OD:'+str(np.round(Dia_case*1000,3))+'mm ID:'+str(np.round(Dia_shaft*1000,3))+'mm')
print('Liner: Height:'+str(np.round(Hei_liner*1000,3))+'mm OD:'+str(np.round(Dia_linOu*1000,3))+'mm ID:'+str(np.round(Dia_linIn*1000,3))+'mm')
A_annOuter = A_ref - A_liner #m2 Outler annulus area
A_annInner = np.pi*(Dia_linIn/2)**2 - A_shaft #m2 Inner annulus area

#### Calculate Liner Length
T4 = Tt4*((1+0.5*(gamma4-1)*(Ma_linerDes**2))**(-1)) #K
SOS4 = np.sqrt(gamma4*R_4*T4) #m/s speed of sound
Vel4 = SOS4*Ma_linerDes #m/s
Len_liner = Vel4*tau_Des #m
print('Liner: Length:'+str(np.round(Len_liner*1000,3))+'mm & designed tau:'+str(np.round(tau_Des*1000,3))+'ms')
Len_PZ = Len_liner*fracPZLength
Len_SZ = Len_liner*fracSZLength
Len_DZ = Len_liner-Len_PZ-Len_SZ
print('PZ_End:'+str(np.round(Len_PZ*1000,3))+'mm')
print('SZ_End:'+str(np.round((Len_PZ+Len_SZ)*1000,3))+'mm')


#### Compute fuel flow rate
mdot_ft = mdot_exit - mdot_at #kg/sec
mdot_apz = mdot_ft/(phi_PZ*FA_st) #kg/sec
mdot_apz_axi = mdot_apz*fracPZAirAxi #kg/sec
mdot_apz_ann = mdot_apz-mdot_apz_axi #kg/sec
mdot_asz = mdot_ft/(phi_SZ*FA_st) - mdot_apz #kg/sec secondary jet flow rate
mdot_adz = mdot_at - mdot_apz - mdot_asz #kg/sec
print('Air Distribution (kg/sec): PZ_axi '+str(np.round(mdot_apz_axi,3))+' PZ_ann '+str(np.round(mdot_apz_ann,3)))
print('Air Distribution (kg/sec): SZ '+str(np.round(mdot_asz,3))+' DZ '+str(np.round(mdot_adz,3)))

#### Compute the hole size
P3 = Pt3*(1+0.5*(gamma3-1)*(Ma_case**2))**(gamma3/(1-gamma3)) #Pa
P_liner = Pt3*(1+0.5*(gamma3-1)*(Ma_linerDes**2))**(gamma3/(1-gamma3)) #Pa
dP_liner = P3-P_liner #Pa Liner pressure drop
rho_liner = rhot3*(1+0.5*(gamma3-1)*(Ma_linerDes**2))**(1/(1-gamma3)) # kg/m3
q_ann = 0.5*rho_liner*(Vel4**2) #Pa annulus dyanmic head 
K_drop = 1+(dP_liner/q_ann) #Pressure drop coefficient
##Compute axial PZ holes
JMF2AMF_axiPZ = mdot_apz_axi/mdot_apz_ann #ratio of jet mass flow to annulus mass flow
JMF2AMF_annPZ = mdot_apz_ann/mdot_apz_axi
JMF2AMF_SZ    = mdot_asz/mdot_apz
JMF2AMF_DZ    = mdot_adz/(mdot_apz+mdot_asz)
if Cd_overwrite:
    Cd_axiPZ = Cd_overwriteVal
    Cd_annPZ = Cd_overwriteVal
    Cd_SZ = Cd_overwriteVal
    Cd_DZ = Cd_overwriteVal
else:
    Cd_axiPZ = (alpha1*(K_drop-1))/(np.sqrt(4*(K_drop**2)-(K_drop*(2-JMF2AMF_axiPZ)**2)))
    Cd_annPZ = (alpha1*(K_drop-1))/(np.sqrt(4*(K_drop**2)-(K_drop*(2-JMF2AMF_annPZ)**2)))
    Cd_SZ = (alpha1*(K_drop-1))/(np.sqrt(4*(K_drop**2)-(K_drop*(2-JMF2AMF_SZ)**2)))
    Cd_DZ = (alpha1*(K_drop-1))/(np.sqrt(4*(K_drop**2)-(K_drop*(2-JMF2AMF_DZ)**2)))
Vel_jet = np.sqrt((2*dP_liner)/rho3) #m/sec
J_mom_axiPZ = (rho3*Vel_jet)/(mdot_apz_ann/A_liner) #Momentum flux ratio
J_mom_annPZ = (rho3*Vel_jet)/(mdot_apz_axi/A_liner)
J_mom_SZ = (rho3*Vel_jet)/(mdot_apz/A_liner)
J_mom_DZ = (rho3*Vel_jet)/((mdot_apz+mdot_asz)/A_liner)
Y_max = Hei_liner*YmaxFrac #m penetration distance of jet
dj_axiPZ = Y_max/(1.25*np.sqrt(J_mom_axiPZ)*ProductToAirRatio) #m
dj_annPZ = Y_max/(1.25*np.sqrt(J_mom_annPZ)*ProductToAirRatio) #m
dj_SZ = Y_max/(1.25*np.sqrt(J_mom_SZ)*ProductToAirRatio) #m
dj_DZ = Y_max/(1.25*np.sqrt(J_mom_DZ)*ProductToAirRatio) #m
N_axiPZ = mdot_apz_axi/(0.25*np.pi*(dj_axiPZ**2)*rho3*Vel_jet)
N_annPZ = mdot_apz_ann/(0.25*np.pi*(dj_annPZ**2)*rho3*Vel_jet)
N_SZ = mdot_asz/(0.25*np.pi*(dj_SZ**2)*rho3*Vel_jet)
N_DZ = mdot_adz/(0.25*np.pi*(dj_DZ**2)*rho3*Vel_jet)
dh_axiPZ = dj_axiPZ/np.sqrt(Cd_axiPZ)
dh_annPZ = dj_annPZ/np.sqrt(Cd_annPZ)
dh_SZ = dj_SZ/np.sqrt(Cd_SZ)
dh_DZ = dj_DZ/np.sqrt(Cd_DZ)

print('end')
