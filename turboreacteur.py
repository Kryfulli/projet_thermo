import numpy as np
import matplotlib.pyplot as plt
#import csv

#variables
Tt4=1600
T0=217
gamma=1.4
cp=1004 #cp massique
PCI=42.8*(10**6)
Patm=1
g=9.81
Pi=1
T_ref = 298.15 #Température de référence
r=287.058
Dn=180
R=8.314
P_ref=100000
M0  =1
z0  = 1000
mode = 'Ts'

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

def P0 (z=10000) :
    return(Patm*np.exp(-7*g*z/(2*cp*T0)))

def Pt0 (M0) :
    return(P0()*((1+(gamma-1)/2*M0*M0)**(gamma/(gamma-1))))

def Pt1 (M0) :
    if M0<1 :
        Pt1=Pt0(M0)*0.95
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
    return((gamma*r*T0)**0.5)

def Fsp (M0,PIc) :
    res = Tt3(M0,PIc)/Tt2(M0)
    res = res -1
    res = Tt0(M0)/Tt4*res
    res = 1-res
    res = res*Tt3(M0,PIc)/Tt2(M0)*Tt0(M0)/T0
    res = res -1
    res = res*Tt4/Tt3(M0,PIc)
    res = res*2/(gamma-1)
    if res<0 :
        return (0)
    res = res**(0.5)
    res = res - M0
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
    r=((2*M0)/(v9(M0,PIc)/a0(M0)+M0))
    if r<1.0:
        return r
    else:
        return 1.0

def rendement_th(M0,PIc) :
    return (1-(1/(Tt3(M0,PIc)/Tt2(M0)*Tt0(M0)/T0)))

def rendement_global(M0,PIc) :
    return (rendement_p(M0,PIc)*rendement_th(M0,PIc))


def diagramme_Fsp_f_Mo () :
    nb_points=200
    M0=np.linspace(0,6,nb_points)
    Fsp_PIc_1=[Fsp(float(k)/float(nb_points-1)*6,1) for k in range (nb_points)]
    Fsp_PIc_3=[Fsp(float(k)/float(nb_points-1)*6,3) for k in range (nb_points)]
    Fsp_PIc_10=[Fsp(float(k)/float(nb_points-1)*6,10) for k in range (nb_points)]
    Fsp_PIc_20=[Fsp(float(k)/float(nb_points-1)*6,20) for k in range (nb_points)]
    plt.plot(M0,Fsp_PIc_1)
    plt.plot(M0,Fsp_PIc_3)
    plt.plot(M0,Fsp_PIc_10)
    plt.plot(M0,Fsp_PIc_20)
    max_fsp=1200
    plt.legend(["PIc=1", "PIc=3","PIc=10","PIc=20"])
    plt.axis([0, 6, 0, max_fsp])
    plt.xlabel('Nombre de Mach en amont')
    plt.ylabel('Poussee specifique N.kg-1.s')
    plt.grid()
    plt.show()


def diagramme_Fsp_f_PIc () :
    nb_points=200
    PIc=np.linspace(0,50,nb_points)
    Fsp_M0_05=[Fsp(0.5,float(k)/float(nb_points-1)*50) for k in range (nb_points)]
    Fsp_M0_1=[Fsp(1,float(k)/float(nb_points-1)*50) for k in range (nb_points)]
    Fsp_M0_2=[Fsp(2,float(k)/float(nb_points-1)*50) for k in range (nb_points)]
    Fsp_M0_3=[Fsp(3,float(k)/float(nb_points-1)*50) for k in range (nb_points)]
    plt.plot(PIc,Fsp_M0_05)
    plt.plot(PIc,Fsp_M0_1)
    plt.plot(PIc,Fsp_M0_2)
    plt.plot(PIc,Fsp_M0_3)
    max_fsp=1200
    plt.legend(["M0=0.5", "M0=1","M0=2","M0=3"])
    plt.axis([0, 50, 0, max_fsp])
    plt.xlabel('Taux de compression')
    plt.ylabel('Poussee specifique N.kg-1.s')
    plt.grid()
    plt.show()


def diagramme_f_f_Mo () :
    nb_points=200
    M0=np.linspace(0,6,nb_points)
    f_PIc_1=[f(float(k)/float(nb_points-1)*6,1) for k in range (nb_points)]
    f_PIc_3=[f(float(k)/float(nb_points-1)*6,3) for k in range (nb_points)]
    f_PIc_10=[f(float(k)/float(nb_points-1)*6,10) for k in range (nb_points)]
    f_PIc_20=[f(float(k)/float(nb_points-1)*6,20) for k in range (nb_points)]
    plt.plot(M0,f_PIc_1)
    plt.plot(M0,f_PIc_3)
    plt.plot(M0,f_PIc_10)
    plt.plot(M0,f_PIc_20)
    max_f=0.035
    plt.legend(["PIc=1", "PIc=3","PIc=10","PIc=20"])
    plt.axis([0, 6, 0, max_f])
    plt.xlabel('Nombre de Mach en amont')
    plt.ylabel('Rapport de melange')
    plt.grid()
    plt.show()

def diagramme_f_f_PIc () :
    nb_points=200
    PIc=np.linspace(0,50,nb_points)
    f_M0_05=[f(0.5,float(k)/float(nb_points-1)*50) for k in range (nb_points)]
    f_M0_1=[f(1,float(k)/float(nb_points-1)*50) for k in range (nb_points)]
    f_M0_2=[f(2,float(k)/float(nb_points-1)*50) for k in range (nb_points)]
    f_M0_3=[f(3,float(k)/float(nb_points-1)*50) for k in range (nb_points)]
    plt.plot(PIc,f_M0_05)
    plt.plot(PIc,f_M0_1)
    plt.plot(PIc,f_M0_2)
    plt.plot(PIc,f_M0_3)
    max_f=0.035
    plt.legend(["M0=0.5", "M0=1","M0=2","M0=3"])
    plt.axis([0, 50, 0, max_f])
    plt.xlabel('Taux de compression')
    plt.ylabel('Rapport de melange')
    plt.grid()
    plt.show()


def diagramme_rendement_f_Mo () :
    nb_points=200
    M0=np.linspace(0,6,nb_points)
    rendement_p_PIc_1=[rendement_p(float(k)/float(nb_points-1)*6,1) for k in range (nb_points)]
    rendement_p_PIc_20=[rendement_p(float(k)/float(nb_points-1)*6,20) for k in range (nb_points)]
    rendement_th_PIc_1=[rendement_th(float(k)/float(nb_points-1)*6,1) for k in range (nb_points)]
    rendement_th_PIc_20=[rendement_th(float(k)/float(nb_points-1)*6,20) for k in range (nb_points)]
    rendement_global_PIc_1=[rendement_global(float(k)/float(nb_points-1)*6,1) for k in range (nb_points)]
    rendement_global_PIc_20=[rendement_global(float(k)/float(nb_points-1)*6,20) for k in range (nb_points)]
    plt.plot(M0,rendement_p_PIc_1)
    plt.plot(M0,rendement_p_PIc_20)
    plt.plot(M0,rendement_th_PIc_1)
    plt.plot(M0,rendement_th_PIc_20)
    plt.plot(M0,rendement_global_PIc_1)
    plt.plot(M0,rendement_global_PIc_20)
    max_rendement=1
    plt.legend(["rendement_p_PIc=1", "rendement_p_PIc=20","rendement_th_PIc=1","rendement_th_PIc=20","rendement_global_PIc=1","rendement_global_PIc=20"])
    plt.axis([0, 6, 0, max_rendement])
    plt.xlabel('Nombre de Mach en amont')
    plt.ylabel('rendement')
    plt.grid()
    plt.show()
    
def P_altitude(z):
    return 100*1013.25*(1-(0.0065*z/288.15))**5.255

def S_k(k,X_k,P,T):
    S=0
    if (k=='O2'):
        a=[3.28253784E+00,1.48308754E-03,-7.57966669E-07,2.09470555E-10,-2.16717794E-14,-1.08845772E+03,5.45323129E+00,3.78245636E+00,-2.99673416E-03,9.84730201E-06,-9.68129509E-09,3.24372837E-12,-1.06394356E+03,3.65767573E+00]     
        if (T>1000):
            S+= R*(a[0]*np.log(T) + a[1]*T + (a[2]*T**2)/2 + (a[3]*T**3)/3 + (a[4]*T**4)/4) + a[6]
        else:
            S+= R*(a[0+7]*np.log(T) + a[1+7]*T + (a[2+7]*T**2)/2 + (a[3+7]*T**3)/3 + (a[4+7]*T**4)/4) + a[6+7]
        return S - R*np.log(P/P_ref) - R*np.log(X_k)
    elif (k=='N2'):
        a=[0.02926640E+02,0.14879768E-02,-0.05684760E-05,0.10097038E-09,-0.06753351E-13,-0.09227977E+04,0.05980528E+02,0.03298677E+02,0.14082404E-02,-0.03963222E-04,0.05641515E-07,-0.02444854E-10,-0.10208999E+04,0.03950372E+02]
        if (T>1000):
            S+= R*(a[0]*np.log(T) + a[1]*T + (a[2]*T**2)/2 + (a[3]*T**3)/3 + (a[4]*T**4)/4) + a[6]
        else:
            S+= R*(a[0+7]*np.log(T) + a[1+7]*T + (a[2+7]*T**2)/2 + (a[3+7]*T**3)/3 + (a[4+7]*T**4)/4) + a[6+7]
        return S - R*np.log(P/P_ref) - R*np.log(X_k)
    elif (k=='air'):
        return 0.21*S_k('O2',0.21*X_k,P,T) + 0.79*S_k('N2',0.79*X_k,P,T)


def S(P,T):
        return S_k('air',1,P,T)
    
def graph_ts(display=False):

    T = T0
    P = P_altitude(z0)

    #A l'entrée, on a T0 et une entropie correspondante à la pression à l'altitude souhaitée
    p_0 = [T,S(P,T)]

    #Le fluide est arrêté de manière isentropique au point 1 avant d'être compressé
    P   = P * np.exp(R/cp*np.log((0.2*(M0**2))))
    T   = T*(1+(0.2*(M0**2)))
    p_1 = [T,p_0[1]]

    #Compresseur isentropique du turboréacteur au point 2
    P   = P*Pi
    T = T*Pi**(r/cp)
    p_3 = [T,p_1[1]]
    
    #Combustion isobare du carburant au point 3
    T   = Tt4
    p_4 = [T,S(P,T)]

    x=[p_0[1],p_1[1],p_3[1]]
    y=[p_0[0],p_1[0],p_3[0]]
    [x.append(i) for i in np.linspace(p_3[1],p_4[1],100)]
    B = np.log(p_4[0]/p_3[0])/(p_4[1]-p_3[1])
    A = p_4[0]/(np.exp(B*p_4[1]))
    [y.append(A*np.exp(B*s)) for s in np.linspace(p_3[1],p_4[1],100)]

    #Détente isentropique au point 4 : même travail que celui du compresseur au point 2
    # on considère que cp = constant --> même longueur de segments
    T = T - (p_3[0]-p_1[0])
    P = P*(T/(Tt4))**(cp/r)
    p_5 = [T,p_4[1]]
    x.append(p_5[1])
    y.append(p_5[0])
    
    #Tuyère au point 9 adapté à la pression atmo
    P = P_altitude(z0)
    T = T0*(P/P_ref)**(R/cp)
    p_9 = [T,p_4[1]]
    x.append(p_9[1])
    y.append(p_9[0])
    if (display):
        fig, ax = plt.subplots()
        ax.set_ylabel('Température (K)')
        ax.set_xlabel('Entropie (J/K)')
        ax.grid(True)
        ax.plot(x,y)
        plt.show()
    return [x,y]

def graph_n(display=False):
    x=np.linspace(0.05,6,200)
    y1=[rendement_th(i,Pi) for i in x]
    y2=[rendement_p(i,Pi) for i in x]
    y3=[rendement_global(i,Pi) for i in x]
    
    return (x,y1,y2,y3)


def get_data():
    if (mode=='Ts'):
        return graph_ts()
    elif (mode=='n'):
        return graph_n()
