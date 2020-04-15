import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, Cursor, MultiCursor,TextBox
import importlib

tr   = importlib.import_module('turboreacteur')
sr   = importlib.import_module('statoreacteur')
#s_tr = importlib.import_module('superturboreacteur')

fig, ax = plt.subplots(figsize=(16,9))

l,    = plt.plot([0,1,2],[0,1.1,3.5])
plt.xlim(-100,2500)
plt.ylim(100,2350)
plt.xlabel('Entropie (kJ)')
plt.ylabel('Température statique (K)')

plt.subplots_adjust(left=0.25,bottom=0.5,right=0.80,top=0.95)

fig.patch.set_facecolor('lightgray')



# Conditions à pouvoir modifier :
# • Boutons :
#       • Modèle (idéal vs réaliste)
#       • Moteur (turbo, stato, turbostato)
#       • Configuration ailaire de l'engin (RadioButton)
#       • Diagramme Ts, ou courbes de rendement
#
# • Sliders :
#       • Altitude
#       • Tt4
#       • Taux de compression (turbo only)
#       • Nombre de Mach



axAlt = plt.axes([0.25, 0.10, 0.55, 0.03], facecolor='white')
axM   = plt.axes([0.25, 0.15, 0.55, 0.03], facecolor='white')
axPi  = plt.axes([0.25, 0.20, 0.55, 0.03], facecolor='white')
axTt4 = plt.axes([0.25, 0.25, 0.55, 0.03], facecolor='white')

axMod = plt.axes([0.025, 0.6, 0.15, 0.15], facecolor='silver')
axMot = plt.axes([0.025, 0.76, 0.15, 0.15], facecolor='silver')

sAlt  = Slider(axAlt,'Altitude (km) :',0,27)
sM    = Slider(axM,'Mach :',0,6)
sPi   = Slider(axPi,'Taux de compression :',1,20)
sTt4  = Slider(axTt4,'Température de combustion (K) :',1600,2200)

cAx   = Cursor(ax,useblit=True,color='red')

rMod  = RadioButtons(axMod, ('Modèle Idéal', 'Modèle Réaliste'), active=0)
rMot  = RadioButtons(axMot, ('Turboréacteur', 'Statoréacteur', 'Superturboréacteur'), active=1)

props = dict(boxstyle='round', facecolor='silver', alpha=0.5)
#tCur  = ax.text(1.01,0.75,'Résultats :\n\nPoussée spécifique :\nConsommation spécifique :    \nRapport de mélange :\nRendements :',bbox=props)


def update_mod(m):
    m.Tt4 = sTt4.val
    m.M0  = sM.val
    m.z0  = sAlt.val*1000
    m.Pi  = sPi.val
    x,y   = m.get_data()
    l.set_xdata(x)
    l.set_ydata(y)
    fig.canvas.draw_idle()

def update(val):
    if (rMot.value_selected=='Turboréacteur'):
        update_mod(tr)
    elif (rMot.value_selected=='Statoréacteur'):
        update_mod(sr)
    elif (rMot.value_selected=='Superturboréacteur'):
        update_mod(s_tr)

sM.on_changed(update)
sAlt.on_changed(update)
sPi.on_changed(update)
sTt4.on_changed(update)
rMod.on_clicked(update)
rMot.on_clicked(update)

update(1)
plt.show()
