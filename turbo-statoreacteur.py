import numpy as np
import matplotlib.pyplot as plt
import csv
import random

#variables
#Tt4=1600

gamma=1.4
#cp=1004 #cp massique
#PCI=42.8*(10**6)
Patm=1.0
g=9.81
r=287.0
R=8.314
Tb=2000.0
T_H2=20.0


#Valeurs numériques
M_H2=2*0.001
M_O2=32*0.001
M_N2=28.0*0.001

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
CFR1=[]
CRS=[]
for row in cr:
    CFR1.append(float(row[0][0:15]))
    CRS.append(float(row[0][18:32]))


#RPR-CFR tous les valeurs qui sont des pourcentages sont multiplies par 100
cr = csv.reader(open("RPR-CFR.csv","rb"))
RPR=[]
CFR2=[]
for row in cr:
    CFR2.append(float(row[0][0:15]))
    RPR.append(float(row[0][18:32]))


def cp_air(T) :
    if T<1000.0 :
        a_N2=[0.02926640E+02,0.14879768E-02,-0.05684760E-05,0.10097038E-09,-0.06753351E-13,-0.09227977E+04,0.05980528E+02]
        a_O2=[3.28253784E+00,1.48308754E-03,-7.57966669E-07,2.09470555E-10,-2.16717794E-14,-1.08845772E+03,5.45323129E+00]
        Cp_N2 = R*(a_N2[0] + a_N2[1]*T + a_N2[2]*(T**2.0) + a_N2[3]*(T**3.0) + a_N2[4]*(T**4.0))
        Cp_O2 = R*(a_O2[0] + a_O2[1]*T + a_O2[2]*(T**2.0) + a_O2[3]*(T**3.0) + a_O2[4]*(T**4.0))
        Cp_air=0.78*Cp_N2/M_N2+0.21*Cp_O2/M_O2
    else :
##        a_N2=[0.03298677E+02,0.14082404E-02,-0.03963222E-04,0.05641515E-07,-0.02444854E-10,-0.10208999E+04,0.03950372E+02]
##        a_O2=[3.78245636E+00,-2.99673416E-03,9.84730201E-06,-9.68129509E-09,3.24372837E-12,-1.06394356E+03,3.65767573E+00]
##        Cp_N2 = R*(a_N2[0] + a_N2[1]*T + a_N2[2]*(T**2.0) + a_N2[3]*(T**3.0) + a_N2[4]*(T**4.0))
##        Cp_O2 = R*(a_O2[0] + a_O2[1]*T + a_O2[2]*(T**2.0) + a_O2[3]*(T**3.0) + a_O2[4]*(T**4.0))
##        Cp_air=0.78*Cp_N2/M_N2+0.21*Cp_O2/M_O2
        return(cp_air(999))
    return(Cp_air)

def cp_H2(T) :
    if T<1000.0 :
        a_H2=[3.33727920E+00,-4.94024731E-05,4.99456778E-07,-1.79566394E-10,2.00255376E-14,-9.50158922E+02,-3.20502331E+00]
        Cp_H2 = R*(a_H2[0] + a_H2[1]*T + a_H2[2]*(T**2.0) + a_H2[3]*(T**3.0) + a_H2[4]*(T**4.0))
    else :
##      a_H2=[2.34433112E+00,7.98052075E-03,-1.94781510E-05,2.01572094E-08,-7.37611761E-12,-9.17935173E+02,6.83010238E-01]
##      Cp_H2 = R*(a_H2[0] + a_H2[1]*T + a_H2[2]*(T**2.0) + a_H2[3]*(T**3.0) + a_H2[4]*(T**4.0))
        return(cp_H2(999))
    return(Cp_H2/M_H2)

def T0(z):
##    a=0.0065
##    Ti=15.0
##    res=Ti-(a*float(z))
##    res=273.15+res
##    if (res<3.0) :
##        res=3.0
##    
##    return(res)
    return(217.0)


def Tt0 (M0,z) : 
    return(T0(z)*(1.0+(gamma-1.0)/2.0*M0*M0))

def P0 (z) :
    return(Patm*np.exp(-7.0*g*z/(2.0*cp_air(T0(z))*T0(z))))

def Pt0 (M0,z) :
    return(P0(z)*((1.0+(gamma-1.0)/2.0*M0*M0)**(gamma/(gamma-1.0))))

def Pt1 (M0,z) :
    if M0<1.0 :
        Pt1=Pt0(M0,z)*0.95
    else :
        Pt1=Pt0(M0,z)*(0.95-0.075*((M0-1)**1.35))
    return (Pt1)

def Tt1 (M0,z) :
    return (Tt0(M0,z)*((Pt0(M0,z)/Pt1(M0,z))**((1.0-gamma)/gamma)))

def T1(M0,z):
    return (Tt1(M0,z)/(1.0+(gamma-1.0)/2.0*(M0**2.0)))

def T2 (M0,z) :
    flux_sortie=cp_air(T1(M0,z))*float(T1(M0,z))+float(flux_1_5)
    
    a=200.0
    b=3500.0
    while abs(b-a)>0.01 :
        m=(a+b)/2.0
        
        if (cp_air(m)*m)>flux_sortie :
            b=m
        else :
            a=m
    return((a+b)/2.0)

def Tt2(M0,z) :
    return (T2(M0,z)*(1.0+(gamma-1.0)/2.0*(M0**2.0)))


def Pt2 (M0,z) :
    return (Pt1(M0,z)*((Tt2(M0,z)/Tt1(M0,z))**(gamma/(gamma-1.0))))


def CRS_to_CFR (CRS_0) :
    k=0
    while (k<len(CRS)-1) and (CRS[k]<CRS_0) :
            k+=1
    if k>=len(CRS)-1 :
        return (CFR1[len(CRS)-1])
    else :
        t=(CRS_0-CRS[k-1])/(CRS[k]-CRS[k-1])
        CFR_0=CFR1[k]*t+CFR1[k-1]*(1-t)
    return (CFR_0)
    
def CFR_to_RPR (CFR_0) :
    k=0
    while (k<len(RPR)-1) and (CFR2[k]<CFR_0)  :
        k+=1
    if k>=len(RPR)-1 :
        return (RPR[len(RPR)-1])
    else :
        t=(CFR_0-CFR2[k-1])/(CFR2[k]-CFR2[k-1])
        RPR_0=RPR[k]*t+RPR[k-1]*(1-t)
        return (RPR_0)

def theta2(M0,z):
    
    return(Tt2(M0,z)/Tn)

def delta2 (M0,z) :
    return (Pt2(M0,z)/Patm)

def CRS_0 (M0,z) :
    return(100.0/((theta2(M0,z))**0.5))

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
    return (Tt2(M0,z)*(PI(M0,z)**((gamma-1.0)/gamma)))

def T3(M0,z):
    return (Tt3(M0,z))

def D_air (M0,z) :
    return (M0/5.0*Dn*CFR_0(M0,z)*delta2 (M0,z)/(theta2(M0,z)**(0.5)))

def flux_massique_3_4(M0,D_H2,z): 
    return ((((D_H2*cp_H2(Tb)+D_air(M0,z)*cp_air(Tb))*Tb)-(D_H2*cp_H2(T_H2)*T_H2)-(D_air(M0,z)*cp_air(T3(M0,z))*T3(M0,z)))/(D_H2+D_air(M0,z)))

def Pt4(M0,z) :
    return (Pt3(M0,z)*((Tb/T3(M0,z))**(gamma/(gamma-1.0))))

def Tt4(M0,z):
    return (Tt3(M0,z)*((Pt4(M0,z)/Pt3(M0,z))**((gamma-1.0)/gamma)))


def T5(M0,D_H2,z) :
    flux_sortie=((cp_air(Tb)*D_air(M0,z)+cp_H2(Tb)*D_H2)*Tb/(D_H2+D_air(M0,z)))+flux_1_5-(flux_massique_3_4(M0,D_H2,z))
    a=0.0
    b=3500.0
    while abs(b-a)>0.01 :
        
        m=(a+b)/2.0
        if cp_air(m)*m>flux_sortie :
            b=m
        else :
            a=m
    return((a+b)/2.0)

def Pt5(M0,D_H2,z) :
    
    
    return (Pt4(M0,z)*((T5(M0,D_H2,z)/Tb)**(gamma/(gamma-1.0))))

def M9 (M0,D_H2,z) :
    
    res=Pt5(M0,D_H2,z)/P0(z)
    res=res**((gamma-1.0)/gamma)
    
    res=2.0*(res-1.0)/(gamma-1.0)
    if res<0 : 
        res=res**0.5
    return(res)



def a0(z) :
    
    return((gamma*r*T0(z))**0.5)

def Fsp (M0,D_H2,z) :
    res=M9(M0,D_H2,z)*(((Tt4(M0,z)/T0(z))/(Tt3(M0,z)/Tt2(M0,z))/(Tt0(M0,z)/T0(z)))**0.5)
    res=a0(z)*(res-M0)
    return(res)

def f (M0,z) :
    res=Tt4(M0,z)/T0(z)-Tt3(M0,z)/Tt2(M0,z)*Tt0(M0,z)/T0(z)
    res=res*T0(z)
    res=res*gamma*r/(gamma-1)
    res=res/PCI
    return(res)
    
def v9 (M0,D_H2,z) :
    return((Fsp(M0,D_H2,z)/a0(z)+M0)*a0(z))

def rendement_p(M0,D_H2,z):
    return((2*M0)/(v9(M0,D_H2,z)/a0(z)+M0))

def rendement_th(M0,z) :
    return (1.0-(1.0/(Tt3(M0,z)/Tt2(M0,z)*Tt0(M0,z)/T0(z))))

def rendement_global(M0,D_H2,z) :
    return (rendement_p(M0,D_H2,z)*rendement_th(M0,z))



dt=1.0
T=600.0
angle_gamma=5.0/360.0*2.0*np.pi
angle_alpha=5.0/360.0*2.0*np.pi
f=10.0
Cz=2.0*np.pi*angle_alpha
Cx=Cz/f
Masse_mol=28.976E-03
Rmax=7000.0
m_engin=200*1000

def masse_vol(z) :
    return (P0(z)*Masse_mol/(R*T0(z)))

def f_D_H2_liste() :
    res=[0.01,0.01]
    for i in range (int(T/dt)) :
        res.append(float(i)*0.005+0.01)
    return (res)


def f_x_z_liste(D_H2_liste):
    x=[0.0,0.0]
    z=[0.0,0.0]
    M=[0.0,0.0]
    t=[0.0,0.0]
    D_H2=f_D_H2_liste()
    for i in  range (int(T/dt)):
        t.append(float(i+1)*dt)
        V0=((((x[i+1]-x[i])/dt)**2.0)+(((z[i+1]-z[i])/dt)**2.0))**0.5
        M.append(V0/(gamma*R/Masse_mol*T0(z[i+1])))
        x.append((dt**2.0)*(-1.0*g*np.sin(angle_alpha)-0.5*Cx*masse_vol(z[i+1])/Rmax*(V0**2)+Fsp(M[i+1],D_H2[i+1],z[i+1])*np.cos(angle_alpha)*D_air(M[i+1],z[i+1])/m_engin)-x[i]+2*x[i+1])
        z.append((dt**2.0)*(-1.0*g*np.cos(angle_alpha)+0.5*Cz*masse_vol(z[i+1])/Rmax*(V0**2)+Fsp(M[i+1],D_H2[i+1],z[i+1])*np.sin(angle_alpha)*D_air(M[i+1],z[i+1])/m_engin)-z[i]+2*z[i+1])
        if z[i+2]<0 :
            z[i+2]=0
    return(x,z,M,t)


delta_variation_D_H2=0.01

def variation_de_liste_D_H2(D_H2_liste,nombre_de_variations):
    res=[]
    for k in range (0,len(D_H2_liste)):
        res.append(D_H2_liste[k]+((random.random()-0.5)*2.0*delta_variation_D_H2*float(nombre_de_variations+1)))
        if res[k]<=0 :
           res[k]=0.01
    return(res)

nombre_de_variations_max=200
nombre_de_changements_de_liste_max=10


def optimization_D_H2_liste () :
    D_H2_liste1=f_D_H2_liste()
    z_f1=f_x_z_liste(D_H2_liste1)[1][int(T/dt)+1]
    M0_f1=f_x_z_liste(D_H2_liste1)[2][int(T/dt)+1]
    nombre_de_variations=0
    nombre_de_changements_de_liste=0
    while (nombre_de_variations_max>nombre_de_variations) and (nombre_de_changements_de_liste_max>nombre_de_changements_de_liste) :
        D_H2_liste2=variation_de_liste_D_H2(D_H2_liste1,nombre_de_variations)
        z_f2=f_x_z_liste(D_H2_liste2)[1][int(T/dt)+1]
        M0_f2=f_x_z_liste(D_H2_liste2)[2][int(T/dt)+1]
        if (abs(z_f1-27000.0)>abs(z_f2-27000.0)) and (abs(z_f1-27000.0)<100.0) :
            if (abs(M0_f1-6.0)>abs(M0_f2-6.0)) and (abs(M0_f1-6.0)<=0.1) :
                if (sum(D_H2_liste1)<sum(D_H2_liste2)) :
                    D_H2_liste1=D_H2_liste2
                    nombre_de_changements_de_liste+=1
                else :
                    nombre_de_variations+=1
            elif (abs(M0_f1-6.0)>abs(M0_f2-6.0)) and (abs(M0_f1-6.0)>0.1) :
                D_H2_liste1=D_H2_liste2
                nombre_de_changements_de_liste+=1
            else :
                nombre_de_variations+=1
        elif (abs(z_f1-27000.0)>abs(z_f2-27000.0)) and (abs(z_f1-27000.0)>100.0) :
            D_H2_liste1=D_H2_liste2
            nombre_de_changements_de_liste+=1
        else :
            nombre_de_variations+=1
        print("nombre_de_variations :",nombre_de_variations)
        print("nombre_de_changements_de_liste :",nombre_de_changements_de_liste)
    return(D_H2_liste1)
        


    

D_H2_liste=optimization_D_H2_liste ()
x_liste=f_x_z_liste(D_H2_liste)[0]
z_liste=f_x_z_liste(D_H2_liste)[1]
M0_liste=f_x_z_liste(D_H2_liste)[2]
t_liste=f_x_z_liste(D_H2_liste)[3]
print(D_H2_liste)
plt.plot(t_liste,z_liste)
plt.show()
