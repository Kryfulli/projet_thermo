# -*- coding: cp1252 -*-
import numpy
import matplotlib.pyplot as plt

Tt = [] #Températures totales 
Ts = [] #Températures statiques
Pt = [] #Pressions totales
Ps = [] #Pressions statiques

T_ref = 298.15 #Température de référence
T0 = 217 #Température statique en entrée d'air
gamma=1.4
cp=1004
PCI=42.8e6
r=287.058
Dn=180
R=8.314
P_ref=100000

def Pt1_Pt0(M): #Rapport Pt1/Pt0, pertes de charges à l'entrée d'air
    if (M<1):
        return 0.95
    else:
        return 0.95-(0.075*(M-1)**1.35)

def Tt(M0,Tt4=1600):
    Tt0=T0*(1+((gamma-1)/2)*M**2)
    Tt1=Tt0*(Pt1_Pt0(M))**(gamma/(gamma-1))
    return [Tt0,Tt1,Tt1,Tt1,Tt4,Tt4,Tt4,Tt4,Tt4,Tt4]

def np(M0,Tt4=1600): #Rendement propulsif
    V9_a0=(2/(gamma-1)*Tt4/T0/((1+((gamma-1)/2)*M0**2))*((1+((gamma-1)/2)*M0**2)-1))**0.5
    return (2*M0)/(V9_a0+M0)

def nth(M0): #Rendement thermique
    return 1-(1/(1+((gamma-1)/2)*M0**2))

def ng(M0,Tt4): #Rendement global
    return np(M0,Tt4)*nth(M0)

def Fsp(M0,Tt4=1600): #poussée spécifique
    V9_a0=(2/(gamma-1)*Tt4/T0/((1+((gamma-1)/2)*M0**2))*((1+((gamma-1)/2)*M0**2)-1))**0.5
    a0=(gamma*r*T0)**0.5
    return a0*(V9_a0-M0)

def f(M0,Tt4=1600): #rapport de mélange
    return (cp*T0*((Tt4/T0)-((Pt1_Pt0(M0))**(gamma/(gamma-1)))))/PCI

def csp(M0,Tt4=1600): #consommation spécifique
    return f(M0,Tt4)/Fsp(M0,Tt4)
    
def graph_Fsp():
    fig, ax = plt.subplots()
    x=numpy.linspace(0.1,7,200)
    y1=Fsp(x,1600)
    y2=Fsp(x,1900)
    y3=Fsp(x,2200)
    ax.plot(x,y1,label='T=1600K')
    ax.plot(x,y2,label='T=1900K')
    ax.plot(x,y3,label='T=2200K')
    ax.legend(loc='upper right', shadow=True, fontsize='x-large')
    ax.set_ylabel('Poussee specifique (s)')
    ax.set_xlabel('Nombre de Mach')
    ax.grid(True)
    plt.show()

def graph_csp():
    x=numpy.linspace(0.1,7,200)
    y1=[csp(i,1600) for i in x]
    y2=[csp(i,1900) for i in x]
    y3=[csp(i,2200) for i in x]
    plt.plot(x,y1,label='T=1600K')
    plt.plot(x,y2)
    plt.plot(x,y3)
    plt.legend("T=1600K","T=1900K","T=2200K")
    plt.show()

def graph_f():
    x=numpy.linspace(0.1,7,200)
    y1=[f(i,1600) for i in x]
    y2=[f(i,1900) for i in x]
    y3=[f(i,2200) for i in x]
    plt.plot(x,y1)
    plt.plot(x,y2)
    plt.plot(x,y3)
    plt.legend("T=1600K","T=1900K","T=2200K")
    plt.show()

def graph_n():
    fig, ax = plt.subplots()
    x=numpy.linspace(0.1,7,200)
    y1_1=[np(i,1600) for i in x]
    y1_2=[nth(i) for i in x]
    y1_3=[ng(i,1600) for i in x]
    y2_1=[np(i,1900) for i in x]
    y2_2=[nth(i) for i in x]
    y2_3=[ng(i,1900) for i in x]
    y3_1=[np(i,2200) for i in x]
    y3_2=[nth(i) for i in x]
    y3_3=[ng(i,2200) for i in x]
    ax.plot(x,y1_1,label='Rendement propulsif, T=1600K')
    ax.plot(x,y1_2,label='Rendement thermique, T=1600K')
    ax.plot(x,y1_3,label='Rendement global, T=1600K')
    ax.plot(x,y2_1,label='Rendement propulsif, T=1900K')
    ax.plot(x,y2_2,label='Rendement thermique, T=1900K')
    ax.plot(x,y2_3,label='Rendement global, T=1900K')
    ax.plot(x,y3_1,label='Rendement propulsif, T=2200K')
    ax.plot(x,y3_2,label='Rendement thermique, T=2200K')
    ax.plot(x,y3_3,label='Rendement global, T=2200K')
    ax.legend(loc='upper left', shadow=True, fontsize='x-large')
    ax.set_ylabel('Rendements')
    ax.set_xlabel('Nombre de Mach')
    ax.grid(True)
    plt.show()



def P_altitude(z):
    return 100*1013.25*(1-(0.0065*z/288.15))**5.255

def S_k(k,X_k,P,T):
    S=0
    if (k=='O2'):
        a=[3.28253784E+00,1.48308754E-03,-7.57966669E-07,2.09470555E-10,-2.16717794E-14,-1.08845772E+03,5.45323129E+00,3.78245636E+00,-2.99673416E-03,9.84730201E-06,-9.68129509E-09,3.24372837E-12,-1.06394356E+03,3.65767573E+00]     
        if (T>1000):
            S+= R*(a[0]*numpy.log(T) + a[1]*T + (a[2]*T**2)/2 + (a[3]*T**3)/3 + (a[4]*T**4)/4) + a[6]
        else:
            S+= R*(a[0+7]*numpy.log(T) + a[1+7]*T + (a[2+7]*T**2)/2 + (a[3+7]*T**3)/3 + (a[4+7]*T**4)/4) + a[6+7]
        return S - R*numpy.log(P/P_ref) - R*numpy.log(X_k)
    elif (k=='N2'):
        a=[0.02926640E+02,0.14879768E-02,-0.05684760E-05,0.10097038E-09,-0.06753351E-13,-0.09227977E+04,0.05980528E+02,0.03298677E+02,0.14082404E-02,-0.03963222E-04,0.05641515E-07,-0.02444854E-10,-0.10208999E+04,0.03950372E+02]
        if (T>1000):
            S+= R*(a[0]*numpy.log(T) + a[1]*T + (a[2]*T**2)/2 + (a[3]*T**3)/3 + (a[4]*T**4)/4) + a[6]
        else:
            S+= R*(a[0+7]*numpy.log(T) + a[1+7]*T + (a[2+7]*T**2)/2 + (a[3+7]*T**3)/3 + (a[4+7]*T**4)/4) + a[6+7]
        return S - R*numpy.log(P/P_ref) - R*numpy.log(X_k)
    elif (k=='air'):
        return 0.21*S_k('O2',0.21*X_k,P,T) + 0.79*S_k('N2',0.79*X_k,P,T)


def S(P,T):
        return cp*numpy.log(T/T0) - R*numpy.log(P/P_altitude(20000))
    
def graph_ts(M0=2,Tt4=1600,display=True):

    T = T0
    P = P_altitude(20000)

    #A l'entrée, on a T0 et une entropie correspondante à la pression à l'altitude souhaitée
    p_0 = [T,S(P,T)]

    #Le fluide est arrêté de manière isentropique au point 3 pour pouvoir faire la combustion.
    P   = P * numpy.exp(R/cp*numpy.log((0.2*(M0**2))))
    T   = T*(1+(0.2*(M0**2)))
    p_3 = [T,p_0[1]]

    #Combustion isobare du carburant.
    T   = Tt4
    p_4 = [T,S(P,T)]

    x=[p_0[1],p_3[1]]
    y=[p_0[0],p_3[0]]
    [x.append(i) for i in numpy.linspace(0,p_4[1],100)]
    [y.append(T0+T0*numpy.exp((s + R*numpy.log(P/P_ref))/cp)) for s in numpy.linspace(0,p_4[1],100)]

    #Tuyère au point 9 adapté à la pression atmo
    P = P_altitude(20000)
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

def graph_hs(M0=2,Tt4=1600):
    l = graph_ts(M0,Tt4,False)
    x = l[0]
    y = [cp*i/1000 for i in l[1]]
    fig, ax = plt.subplots()
    ax.set_ylabel('Enthalpie (kJ)')
    ax.set_xlabel('Entropie (J/K)')
    ax.grid(True)
    ax.plot(x,y)
    plt.show()
    return [x,y]


#hs=cp*(T-T_ref)
#Ds = cp ln(T/T_ref) - R ln(P/P_ref)
#lors d'une isentrope, ln(T2/T1) = R/Cp * ln(P2/P1)
#S_k = S_k_chimique - R*ln(P/P_ref) - R*ln(X_k)
#S = somme(S_k * X_k) -> s = S/M_mélange
#14 coeffs : 7 premiers pour calculer les HT, de 1000 à 3500K
# 7 seconds pour de 200 à 1000K.
#rajouter un échangeur à l'entrée du compresseur du turboréacteur
#webplotanalyser

#ecrire bilan sur chambre de combustion
