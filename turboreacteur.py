import numpy as np
import matplotlib.pyplot as plt
import csv

#variables
Tt4=1600
T0=217
gamma=1.4
cp=1004 #cp massique
PCI=42.8*(10**6)
Patm=1
g=9.81
r=287

#conditions nominales
Dn=180
Tn=288
Pn=1
PIn=6.67

#definition des listes
P=[]
T=[]
Pt=[]
Tt=[]

###CRS-CFR
##cr = csv.reader(open("CRS-CFR.csv","rb"))
##CFR=[]
##CRS=[]
##for row in cr:
##    CFR.append(float(row[0][0:15]))
##    CRS.append(float(row[0][18:32]))
##
##
###RPR-CFR tous les valeurs qui sont des pourcentages sont multiplies par 100
##cr = csv.reader(open("RPR-CFR.csv","rb"))
##RPR=[]
##CFR=[]
##for row in cr:
##    CFR.append(float(row[0][0:15]))
##    RPR.append(float(row[0][18:32]))




def Tt0 (M0) : 
    return(T0*(1+(gamma-1)/2*M0*M0))

def P0 (z=3000) :
    return(Patm*np.exp(-7*g*z/(2*cp*T0)))

def Pt0 (M0) :
    return(P0()*((1+(gamma-1)/2*M0*M0)**(gamma/(gamma-1))))

def Pt1 (M0) :
    if M0<1 :
        Pt1=Pt0()*0.95
    else :
        Pt1=Pt0(M0)*(0.95-0.075*((M0-1)**1.35))
    return (Pt1)

def Tt1 (M0) :
    return (Tt0(M0)*((Pt0(M0)/Pt1(M0))**((1-gamma)/gamma)))

def Tt2 (M0) :
    return (Tt1(M0))

def Pt2 (M0) :
    return (Pt1(M0))


##def CRS_to_CFR (CRS_0) :
##    k=1
##    while CRS[k]<CRS_0 :
##        k+=1
##    t=(CRS_0-CRS[k-1])/(CRS[k]-CRS[k-1])
##    CFR_0=CFR[k]*t+CFR[k-1]*(1-t)
##    return (CFR_0)
##    
##def CFR_to_RPR (CFR_0) :
##    k=1
##    while CFR[k]<CFR_0 :
##        k+=1
##    t=(CFR_0-CFR[k-1])/(CFR[k]-CFR[k-1])
##    RPR_0=RPR[k]*t+RPR[k-1]*(1-t)
##    return (RPR_0)
##
##def theta2(M0):
##    return(Tt2(M0)/Tn)
##
##def CRS_0 (M0) :
##    return(100/((theta2(M0))**0.5))
##
##def RPR_0 (M0) :
##    RPR_0=CFR_to_RPR(CRS_to_CFR(CRS_0(M0)))
##    return(RPR_0)
##
##def PI (M0,PIn) :
##    return (PIn*RPR_0(M0))

def Pt3(M0,PIc) :
    return (Pt2(M0)*PIc)

def Tt3(M0,PIc) :
    return (Tt2(M0)*(PIc**((gamma-1)/gamma)))

def a0(M0) :
    return(M0*((gamma*r*T0)**0.5))

def Fsp (M0,PIc) :
    res = Tt3(M0,PIc)/Tt2(M0)
    res = res -1
    res = Tt0(M0)/Tt4*res
    res = 1-res
    res = res*Tt3(M0,PIc)/Tt2(M0)*Tt0(M0)/T0
    res = res -1
    res = res*Tt4/Tt3(M0,PIc)
    res = res*2/(gamma-1)
    res =res**0.5
    res =res - M0
    res = res *a0(M0)
    return(res)

def f (M0,PIc) :
    res=Tt4/T0-Tt3(M0,PIc)/Tt2(M0)*Tt0(M0)/T0
    res=res*T0
    res=res*gamma*r/(gamma-1)
    res=res/PCI
    return(res)
    
def v9 (M0,PIc) :
    return((Fsp(M0,PIc)/a0(M0)+M0)*a0(M0))

def rendement_p(M0,PIc):
    return((2*M0)/(v9(M0,PIc)/a0(M0)+M0))

def rendement_th(M0,PIc) :
    return (1-(1/(Tt3(M0,PIc)/Tt2(M0)*Tt0(M0)/T0)))

    
