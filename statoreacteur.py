# -*- coding: cp1252 -*-
import numpy
import matplotlib.pyplot as plt

Tt = [] #Températures totales 
Ts = [] #Températures statiques
Pt = [] #Pressions totales
Ps = [] #Pressions statiques

T0 = 217 #Température statique en entrée d'air
gamma=1.4
cp=1004
PCI=42.8e6
r=287.058
Dn=180

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

def isobare(T,s,s_depart=0):
    C1=-cp*numpy.log(T)
    return cp*(numpy.exp(s/cp)-numpy.exp(s_depart/cp))*numpy.exp(-C1/cp)


def graph_hs(M0,Tt4=1600):
    #P_i : s,h
    M0=1.2
    Tt4=1600
    C1=-cp*numpy.log(T0)
    p_0=[0,cp*T0]
    p_3=[0,cp*T0+(0.5*(M0**2)*gamma*r*T0)]
    h4=p_3[1]+(f(M0,Tt4)*Dn*PCI)
    s4=numpy.log((numpy.exp(C1/cp)/cp*h4)+1)*cp
    p_4=[s4,h4]
    p_9=[s4,isobare(T0,s4)]
    x=[0]
    x=x+numpy.linspace(0,s4,100)
    x=x+[s4]
    print(len(x))
    y=[p_0[1]]
    y=y+[isobare(Tt4,i) for i in numpy.linspace(0,s4,100)]
    print(len(y))
    y.append(p_9[1])
    fig, ax = plt.subplots()
    ax.plot(x,y)
    plt.show()


M0=1.2
Tt4=1600
p_0=[0,cp*T0]
p_3=[0,cp*T0+(0.5*(M0**2)*gamma*r*T0)]
C1=-cp*numpy.log(T0)
h4=p_3[1]+(f(M0,Tt4)*Dn*PCI)
x=[p_0[0],p_3[0]]
y=[p_0[1],p_3[1]]
s4=numpy.log((numpy.exp(C1/cp)/cp*h4)+1)*cp
[x.append(i) for i in numpy.linspace(0,s4,100)]
[y.append(isobare(Tt4,i)+p_3[1]) for i in numpy.linspace(0,s4,100)]
fig, ax = plt.subplots()
ax.plot([i/1000 for i in x],[i/1000 for i in y])
plt.show()




    
