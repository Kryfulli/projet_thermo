import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, Cursor, MultiCursor
import importlib

tr   = importlib.import_module('turboreacteur')
sr   = importlib.import_module('statoreacteur')
#s_tr = importlib.import_module('superturboreacteur')

fig, ax = plt.subplots(figsize=(16,9))
ax.autoscale(True)
l,    = plt.plot([],[]) #Diagramme T,s

ln1,  = plt.plot([],[], label='Rendement propulsif', linestyle='dashed')
ln2,  = plt.plot([],[], label='Rendement thermique', linestyle='dashed') # Courbes de rendement
ln3,  = plt.plot([],[], label='Rendement global')

ls1,  = plt.plot([],[], label='Consommation Spécifique')
ls2,  = plt.plot([],[], label='Poussée Spécifique')

le1    = plt.legend((ls1,ls2),('Consommation Spécifique','Poussée Spécifique'),loc='upper left')
le1    = plt.legend((ln1,ln2,ln3),('Rendement propulsif','Rendement thermique','Rendement global'),loc='upper left')

le1.set_visible(False)

ax.set_xlim(150,220)
ax.set_ylim(100,2350)
plt.xlabel('Entropie (J)')
plt.ylabel('Température statique (K)')
plt.grid(True)
plt.subplots_adjust(left=0.25,bottom=0.5,right=0.80,top=0.95)

fig.patch.set_facecolor('lightgray')



# Conditions à pouvoir modifier :
# • Boutons :
#       • Modèle (idéal vs réaliste)
#       • Moteur (turbo, stato, turbostato)
#       • Diagramme Ts, ou courbes de rendement
#
# • Sliders :
#       • Altitude
#       • Tt4
#       • Taux de compression (turbo only)
#       • Nombre de Mach



axAlt = plt.axes([0.250, 0.10, 0.55, 0.03], facecolor='white')
axM   = plt.axes([0.250, 0.15, 0.55, 0.03], facecolor='white')
axPi  = plt.axes([0.250, 0.20, 0.55, 0.03], facecolor='white')
axTt4 = plt.axes([0.250, 0.25, 0.55, 0.03], facecolor='white')

axMod = plt.axes([0.025, 0.60, 0.15, 0.15], facecolor='silver')
axMot = plt.axes([0.025, 0.76, 0.15, 0.15], facecolor='silver')
axGra = plt.axes([0.025, 0.44, 0.15, 0.15], facecolor='silver')

sAlt  = Slider(axAlt,'Altitude (km) :',0,27)
sM    = Slider(axM,'Mach :',0,6)
sPi   = Slider(axPi,'Taux de compression :',1,20)
sTt4  = Slider(axTt4,'Température de combustion (K) :',1600,2200)

cAx   = Cursor(ax,useblit=True,color='red')

rMod  = RadioButtons(axMod, ('Modèle Idéal', 'Modèle Réaliste'), active=0)
rMot  = RadioButtons(axMot, ('Turboréacteur', 'Statoréacteur', 'Superturboréacteur'), active=1)
rGra  = RadioButtons(axGra, ('Diagramme T,s', 'Courbes de rendement', 'Poussée Spécifique','Conso spécifique'), active=0)



def update_mod(m):
    m.Tt4  = sTt4.val
    m.M0   = sM.val
    m.z0   = sAlt.val*1000
    m.Pi   = sPi.val
    if (rGra.value_selected=='Diagramme T,s'):
        m.mode = 'Ts'
        x,y    = m.get_data()
        l.set_xdata(x)
        l.set_ydata(y)
        ax.set_xlim(l.get_xdata()[0]-10,l.get_xdata()[-1]+10)
        ax.set_ylim(100,max(2350,l.get_ydata()[2]+100))
        le1.set_visible(False)
        fig.canvas.draw_idle()
    elif (rGra.value_selected=='Courbes de rendement'):
        m.mode = 'n'
        x, y1, y2, y3 = m.get_data()
        ln1.set_xdata(x)
        ln1.set_ydata(y1)
        ln2.set_xdata(x)
        ln2.set_ydata(y2)
        ln3.set_xdata(x)
        ln3.set_ydata(y3)
        le1.set_visible(True)
        fig.canvas.draw_idle()
    elif (rGra.value_selected=='Poussée Spécifique'):
        m.mode = 'P_s'
        x, y1 = m.get_data()
        ls1.set_xdata(x)
        ls1.set_ydata(y1)
        le1.set_visible(False)
        fig.canvas.draw_idle()
    else:
        m.mode = 'C_s'
        x, y1 = m.get_data()
        ls1.set_xdata(x)
        ls1.set_ydata(y1)
        le1.set_visible(False)
        fig.canvas.draw_idle()
        
def update(val):
    if (rMot.value_selected=='Turboréacteur'):
        update_mod(tr)
    elif (rMot.value_selected=='Statoréacteur'):
        update_mod(sr)
    elif (rMot.value_selected=='Superturboréacteur'):
        update_mod(s_tr)

def switch_mode(val):
    if (rGra.value_selected=='Diagramme T,s'):
        ln1.set_visible(False)
        ln2.set_visible(False)
        ln3.set_visible(False)
        l.set_visible(True)
        sM.set_active(True)
        sM.set_val(0.5)
        axM.set_facecolor('white')
        ax.set_xlim(l.get_xdata()[0]-10,l.get_xdata()[-1]+10)
        ax.set_ylim(100,max(2350,l.get_ydata()[2]+100))
        ax.set_xlabel('Entropie (J)')
        ax.set_ylabel('Température statique (K)')
        update(1)
        
    elif (rGra.value_selected=='Courbes de rendement'):
        ln1.set_visible(True)
        ln2.set_visible(True)
        ln3.set_visible(True)
        l.set_visible(False)
        sM.set_active(False)
        sM.set_val(0)
        axM.set_facecolor('gray')
        ax.set_xlabel('Nombre de Mach')
        ax.set_ylabel('Rendements')
        ax.set_xlim(0,6)
        ax.set_ylim(0,1)
        update(1)
    elif (rGra.value_selected=='Poussée Spécifique'):
        ln1.set_visible(False)
        ln2.set_visible(False)
        ln3.set_visible(False)
        l.set_visible(False)
        ax.set_xlabel('Nombre de Mach')
        ax.set_ylabel('Poussée spécifique (s)')
        ax.set_xlim(0,6)
        ax.set_ylim(0,800)
        update(1)
    else:
        ln1.set_visible(False)
        ln2.set_visible(False)
        ln3.set_visible(False)
        l.set_visible(False)
        ax.set_xlabel('Nombre de Mach')
        ax.set_ylabel('Consommation spécifique (kg.s-1)')
        ax.set_xlim(0,6)
        ax.set_ylim(0,.001)
        update(1)
        
def switch_mot(val):
    if (rMot.value_selected=='Statoréacteur'):
        sPi.set_active(False)
        sPi.set_val(1)
        axPi.set_facecolor('gray')
        update(1)
    else:
        sPi.set_active(True)
        sPi.set_val(5)
        axPi.set_facecolor('white')
        update(1)

sM.on_changed(update)
sAlt.on_changed(update)
sPi.on_changed(update)
sTt4.on_changed(update)
rMod.on_clicked(update)
rMot.on_clicked(switch_mot)
rGra.on_clicked(switch_mode)
sPi.set_active(False)
sPi.set_val(1)
axPi.set_facecolor('gray')
update(1)
plt.show()
