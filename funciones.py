""" 
ATMÓSFERAS ESTELARES - ANÁLISIS DE DOS MODELOS DE ATMÓSFERA

Módulo con las funciones utilizadas para el análisis y la elaboración
de figuras.

Autores: Víctor Alonso, Alejandro Cerón
"""

### Imports:
#############################################################################################
import os
import shutil
import sys

import pdb

import numpy as np
from astropy.constants import M_sun, R_sun, G, R, h, m_e, m_p, c, sigma_sb, k_B
from astropy.io.ascii import read
from astropy import units as u
from astropy.table import QTable

import pandas as pd
pd.set_option('mode.chained_assignment', None)

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import matplotlib.pyplot as plt
import matplotlib.lines as mlines

### Utilizar LaTeX en las figuras:
#############################################################################################
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif') 
plt.rc('font', size=16)  

plt.rcParams['text.latex.preamble'] = r'''
            \usepackage{siunitx}
            \sisetup{
              detect-family,
              separate-uncertainty=true,
              output-decimal-marker={.},
              exponent-product=\cdot,
              inter-unit-product=\cdot,
            }
            \DeclareSIUnit{\cts}{cts}
            \DeclareSIUnit{\year}{yr}
            \DeclareSIUnit{\dyn}{dyn}
            \DeclareSIUnit{\mag}{mag}
            \usepackage{sansmath}  % Allows sans-serif in math mode
            \sansmath
            '''
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Computer Modern Serif",
})

#############################################################################################
def plot_lgtauR(tablas:list, axis, figure, save_path=None, show=True, im_format="pdf"):
    #for tabla in tablas:
        #axis.plot(depth, plot_lgtauR)
    pass
