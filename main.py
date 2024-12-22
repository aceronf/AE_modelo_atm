""" 
ATMÓSFERAS ESTELARES - ANÁLISIS DE DOS MODELOS DE ATMÓSFERA

Descripción general del programa. Información acerca de los modelos en 
https://marcs.astro.uu.se/docs.html

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

### Ficheros con los modelos:
#############################################################################################
t_5000 = "t5000.dat"
t_8000 = "t8000.dat"

### Función para leer las tablas:
#############################################################################################
def read_table(model_name:str):
    # Leemos la tabla y obtenemos una Qtable para tener unidades en las columnas:
    table = read (model_name, header_start = 24, data_start = 25) 
    table = QTable(table)
    # Adjudicamos unidades:
    table["k"].unit = u.dimensionless_unscaled # Numeral del punto
    table["lgTauR"].unit = u.dimensionless_unscaled # Log de la opacidad de Rosseland
    table["lgTau5"].unit = u.dimensionless_unscaled # Log de la opacidad a 5000 A
    table["Depth"].unit = u.cm # Profundidad en cm. Valor de 0 cuando TauR = 1
    table["T"].unit = u.Kelvin # Temperatura en Kelvin
    table["Pe"].unit = u.dyn / u.cm**2 # Presión electrónica
    table["Pg"].unit = u.dyn / u.cm**2 # Presión del gas
    table["Prad"].unit = u.dyn / u.cm**2 # Presión de la radiación
    table["Pturb"].unit = u.dyn / u.cm**2 # Presión turbulenta
    return table

### Programa principal:
#############################################################################################
if __name__ == "__main__":

    ### Se leen los datos y se convierten a tablas de astropy
    #########################################################################################
    t_5000_table = read_table(t_5000)
    t_8000_table = read_table(t_8000)
    #print(t_5000_table)

    ### Se crea un directorio donde guardar los resultados:
    #########################################################################################
    results_dir = "resultados"
    # Creamos la carpeta con los plots:
    if os.path.exists(results_dir):
            shutil.rmtree(results_dir)    
    os.makedirs(results_dir)

    ### Plot de lgTauR frente a la profundidad:
    #########################################################################################
    fig1, ax1 = plt.subplots(figsize=(15, 8))
    ax1.plot(t_5000_table["Depth"].value, t_5000_table["lgTauR"].value, color="purple")
    ax1.plot(t_8000_table["Depth"].value, t_8000_table["lgTauR"].value, color="orange")
    plt.show()





