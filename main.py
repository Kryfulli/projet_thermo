import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, Cursor, MultiCursor,TextBox
import importlib

tr = importlib.import_module('turboreacteur')
sr = importlib.import_module('statoreacteur')

fig, ax = plt.subplots(figsize=(16,9))
plt.subplots_adjust(left=0.25,bottom=0.5,right=0.80,top=0.95)

axcolor = 'lightgoldenrodyellow'




# Conditions à pouvoir modifier :
# • Boutons :
#       • Modèle (idéal vs réaliste)
#       • Moteur (turbo, stato, turbostato)
#       • Configuration ailaire de l'engin (RadioButton)
# • Sliders :
#       • Altitude
#       • Tt4
#       • Taux de compression (turbo only)
#       • Nombre de Mach



axAlt = plt.axes([0.25, 0.10, 0.65, 0.03], facecolor=axcolor)
axM   = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
axPi  = plt.axes([0.25, 0.20, 0.65, 0.03], facecolor=axcolor)
axTt4 = plt.axes([0.25, 0.25, 0.65, 0.03], facecolor=axcolor)

axMod = plt.axes([0.025, 0.6, 0.15, 0.15], facecolor=axcolor)
axMot = plt.axes([0.025, 0.76, 0.15, 0.15], facecolor=axcolor)

sAlt  = Slider(axAlt,'Altitude (km) :',0,27)
sM    = Slider(axM,'Mach :',0,6)
sPi   = Slider(axPi,'Taux de compression :',1,20)
sTt4  = Slider(axTt4,'Température de combustion (K) :',1600,2200)

cAx   = Cursor(ax,useblit=True,color='red')

rMod  = RadioButtons(axMod, ('Modèle Idéal', 'Modèle Réaliste'), active=0)
rMot  = RadioButtons(axMot, ('Turboréacteur', 'Statoréacteur', 'Turbostatoréacteur'), active=0)

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
tCur  = ax.text(1.01,0.75,'Résultats :\n\nPoussée spécifique :\nConsommation spécifique :    \nRapport de mélange :\nRendements :',bbox=props)



plt.show()
