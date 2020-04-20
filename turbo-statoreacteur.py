import numpy as np
import matplotlib.pyplot as plt
import csv

#variables
#Tt4=1600
T0=217.0
gamma=1.4
cp=1004 #cp massique
PCI=42.8*(10**6)
Patm=1.0
g=9.81
r=287.0
R=8.314
Tb=2000.0
T_H2=20.0

#flux echangeur de chaleur
flux_1_5=50000.0 #J


#conditions nominales
Dn=180.0
Tn=288.0
Pn=1.0
PIn=6.67

#definition des listes
P=[]
T=[]
Pt=[]
Tt=[]

#CRS-CFR
cr = csv.reader(open("CRS-CFR.csv","rb"))
CFR=[]
CRS=[]
for row in cr:
    CFR.append(float(row[0][0:15]))
    CRS.append(float(row[0][18:32]))


#RPR-CFR tous les valeurs qui sont des pourcentages sont multiplies par 100
cr = csv.reader(open("RPR-CFR.csv","rb"))
RPR=[]
CFR=[]
for row in cr:
    CFR.append(float(row[0][0:15]))
    RPR.append(float(row[0][18:32]))


def cp_air(T) :
    if T<1000 :
        a_N2=[0.02926640E+02,0.14879768E-02,-0.05684760E-05,0.10097038E-09,-0.06753351E-13,-0.09227977E+04,0.05980528E+02]
        a_O2=[3.28253784E+00,1.48308754E-03,-7.57966669E-07,2.09470555E-10,-2.16717794E-14,-1.08845772E+03,5.45323129E+00]
        Cp_N2 = R*(a_N2[0] + a_N2[1]*T + a_N2[2]*(T**2) + a_N2[3]*(T**3) + a_N2[4]*(T**4))
        Cp_O2 = R*(a_O2[0] + a_O2[1]*T + a_O2[2]*(T**2) + a_O2[3]*(T**3) + a_O2[4]*(T**4))
        Cp_air=0.78*Cp_N2*14+0.21*Cp_O2*32
    else :
        a_N2=[0.03298677E+02,0.14082404E-02,-0.03963222E-04,0.05641515E-07,-0.02444854E-10,-0.10208999E+04,0.03950372E+02]
        a_O2=[3.78245636E+00,-2.99673416E-03,9.84730201E-06,-9.68129509E-09,3.24372837E-12,-1.06394356E+03,3.65767573E+00]
        Cp_N2 = R*(a_N2[0] + a_N2[1]*T + a_N2[2]*(T**2) + a_N2[3]*(T**3) + a_N2[4]*(T**4))
        Cp_O2 = R*(a_O2[0] + a_O2[1]*T + a_O2[2]*(T**2) + a_O2[3]*(T**3) + a_O2[4]*(T**4))
        Cp_air=0.78*Cp_N2*14+0.21*Cp_O2*32
    return(Cp_air)

def cp_H2(T) :
    if T<1000 :
        a_H2=[3.33727920E+00,-4.94024731E-05,4.99456778E-07,-1.79566394E-10,2.00255376E-14,-9.50158922E+02,-3.20502331E+00]
        Cp_H2 = R*(a_H2[0] + a_H2[1]*T + a_H2[2]*(T**2) + a_H2[3]*(T**3) + a_H2[4]*(T**4))
    else :
        a_H2=[2.34433112E+00,7.98052075E-03,-1.94781510E-05,2.01572094E-08,-7.37611761E-12,-9.17935173E+02,6.83010238E-01]
        Cp_H2 = R*(a_H2[0] + a_H2[1]*T + a_H2[2]*(T**2) + a_H2[3]*(T**3) + a_H2[4]*(T**4))
    return(Cp_H2)


def Tt0 (M0) : 
    return(T0*(1+(gamma-1)/2*M0*M0))

def P0 (z) :
    return(Patm*np.exp(-7*g*z/(2*cp*T0)))

def Pt0 (M0,z) :
    return(P0(z)*((1+(gamma-1)/2*M0*M0)**(gamma/(gamma-1))))

def Pt1 (M0,z) :
    if M0<1 :
        Pt1=Pt0(M0,z)*0.95
    else :
        Pt1=Pt0(M0,z)*(0.95-0.075*((M0-1)**1.35))
    return (Pt1)

def Tt1 (M0,z) :
    return (Tt0(M0)*((Pt0(M0,z)/Pt1(M0,z))**((1-gamma)/gamma)))

def T1(M0,z):
    return (Tt1(M0,z)/(1+(gamma-1)/2*(M0**2)))

def T2 (M0,z) :
    flux_sortie=cp_air(T1(M0,z))*float(T1(M0,z))+float(flux_1_5)
    
    a=200.0
    b=3500.0
    while abs(b-a)>0.01 :
        m=(a+b)/2
        
        if (cp_air(m)*m)>flux_sortie :
            b=m
        else :
            a=m
    return((a+b)/2)

def Tt2(M0,z) :
    return (T2(M0,z)*(1+(gamma-1)/2*(M0**2)))


def Pt2 (M0,z) :
    return (Pt1(M0,z)*(Tt2(M0,z)/Tt1(M0,z))**(gamma/(gamma-1)))


def CRS_to_CFR (CRS_0) :
    k=1
    while CRS[k]<CRS_0 :
        k+=1
    t=(CRS_0-CRS[k-1])/(CRS[k]-CRS[k-1])
    CFR_0=CFR[k]*t+CFR[k-1]*(1-t)
    return (CFR_0)
    
def CFR_to_RPR (CFR_0) :
    k=1
    while CFR[k]<CFR_0 :
        k+=1
    t=(CFR_0-CFR[k-1])/(CFR[k]-CFR[k-1])
    RPR_0=RPR[k]*t+RPR[k-1]*(1-t)
    return (RPR_0)

def theta2(M0,z):
    return(Tt2(M0,z)/Tn)

def delta2 (M0,z) :
    return (Pt2(M0,z)/Patm)

def CRS_0 (M0,z) :
    return(100/((theta2(M0,z))**0.5))

def RPR_0 (M0,z) :
    RPR_0=CFR_to_RPR(CRS_to_CFR(CRS_0(M0,z)))
    return(RPR_0)

def CFR_0(M0,z):
    CFR_0=CRS_to_CFR (CRS_0(M0,z))
    return (CFR_0)

def PI (M0,z) :
    return (PIn*RPR_0(M0,z))

def Pt3(M0,z) :
    return (Pt2(M0,z)*PI(M0,z))

def Tt3(M0,z) :
    return (Tt2(M0,z)*(PI(M0,z)**((gamma-1)/gamma)))

def T3(M0,z):
    return (Tt3(M0,z))

def D_air (M0,z) :
    return (Dn*CFR_0(M0,z)*delta2 (M0,z)/(theta2(M0,z)**(0.5)))

def flux_massique_3_4(M0,D_H2,z): 
    return ((((D_H2*cp_H2(Tb)+D_air(M0,z)*cp_air(Tb))*Tb)-(D_H2*cp_H2(T_H2)*T_H2)-(D_air(M0,z)*cp_air(T3(M0,z))*T3(M0,z)))/(D_H2+D_air(M0,z)))

def Pt4(M0,z) :
    return (Pt3(M0,z)*((Tb/T3(M0,z))**(gamma/(gamma-1))))

def Tt4(M0,z):
    return (Tt3(M0,z)*((Pt4(M0,z)/Pt3(M0,z))**((gamma-1)/gamma)))


def T5(M0,D_H2,z) :
    flux_sortie=((cp_air(Tb)*D_air(M0,z)+cp_H2(Tb)*D_H2)*Tb/(D_H2+D_air(M0,z)))+flux_1_5-(flux_massique_3_4(M0,D_H2,z))
    a=0.0
    b=3500.0
    while abs(b-a)>0.01 :
        
        m=(a+b)/2
        if cp_air(m)*m>flux_sortie :
            b=m
        else :
            a=m
    return((a+b)/2)

def Pt5(M0,D_H2,z) :
    
    
    return (Pt4(M0,z)*((T5(M0,D_H2,z)/Tb)**(gamma/(gamma-1))))

def M9 (M0,D_H2,z) :
    
    res=Pt5(M0,D_H2,z)/P0(z)
    res=res**((gamma-1)/gamma)
    
    res=2*(res-1)/(gamma-1)
    
    res=res**0.5
    return(res)



def a0(M0) :
    return((gamma*r*T0)**0.5)

def Fsp (M0,D_H2,z) :
    res=M9(M0,D_H2,z)*(((Tt4(M0,z)/T0)/(Tt3(M0,z)/Tt2(M0,z))/(Tt0(M0)/T0))**0.5)
    res=a0(M0)*(res-M0)
    return(res)

def f (M0,z) :
    res=Tt4(M0,z)/T0-Tt3(M0,z)/Tt2(M0,z)*Tt0(M0,z)/T0
    res=res*T0
    res=res*gamma*r/(gamma-1)
    res=res/PCI
    return(res)
    
def v9 (M0,D_H2,z) :
    return((Fsp(M0,D_H2,z)/a0(M0)+M0)*a0(M0))

def rendement_p(M0,D_H2,z):
    return((2*M0)/(v9(M0,D_H2,z)/a0(M0)+M0))

def rendement_th(M0,z) :
    return (1-(1/(Tt3(M0,z)/Tt2(M0,z)*Tt0(M0,z)/T0)))

def rendement_global(M0,D_H2,z) :
    return (rendement_p(M0,D_H2,z)*rendement_th(M0,z))



