'''
ExoSpin run script - Diagram Plot


@authors : I. Abdoulwahab & P. Palma-Bifani & G. Chauvin & A. Simonnin

'''

# ------------------------------------------------------------------------------------------------------------------------------------------------
## Imports

import matplotlib.pyplot as plt
import numpy as np


def diagram_plot(exoplanet_name, io, ip, omega_o = 0):
    star = plt.Circle((0.5,0.5), 0.2, color = '#ffd319')
    star_stroke = plt.Circle((0.5,0.5), 0.3, color = '##bb6f1e')