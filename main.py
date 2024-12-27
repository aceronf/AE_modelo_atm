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
from astropy.constants import Ryd, h, m_e, m_p, c ,k_B
from astropy.io.ascii import read
from astropy import units as u
from astropy.table import QTable

import pandas as pd
pd.set_option('mode.chained_assignment', None)

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from figuras import plot_lgtauR, plot_T, plot_Pe, plot_Pe_Pg, plot_Prad_Pg
from poblaciones import poblaciones

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

plot_params = {
     "label":[r"$T_{\mathrm{eff}}=5000 \, \unit{\kelvin}$", 
              r"$T_{\mathrm{eff}}=8000 \, \unit{\kelvin}$" ],
     "color":["red", "blue"],
     "linestyle":["-", "--"],
     "Teff":[5000, 8000]
                }


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

    ### 1) Plot de lgTauR frente a la profundidad:
    #########################################################################################
    fig1, ax1 = plt.subplots(figsize=(15, 8))
    plot_lgtauR((t_5000_table, t_8000_table), params=plot_params, axis=ax1, figure=fig1,
                save_path=os.path.join(results_dir,"lgTauR_R.pdf"))

    ### 2) Plot T, Pe, Pe/Pg y Prad/Pg frente a lgTauR
    #########################################################################################
    fig2_1, ax2_1 = plt.subplots(figsize=(15, 8))
    plot_T((t_5000_table, t_8000_table), params=plot_params, axis=ax2_1, figure=fig2_1,
                save_path=os.path.join(results_dir,"T_lgTauR.pdf"))
    

    fig2_2, ax2_2 = plt.subplots(figsize=(15, 8))
    plot_Pe((t_5000_table, t_8000_table), params=plot_params, axis=ax2_2, figure=fig2_2,
                save_path=os.path.join(results_dir,"Pe_lgTauR.pdf"))
    
    fig2_3, ax2_3 = plt.subplots(figsize=(15, 8))
    plot_Pe_Pg((t_5000_table, t_8000_table), params=plot_params, axis=ax2_3, figure=fig2_3,
                save_path=os.path.join(results_dir,"Pe_Pg_lgTauR.pdf"))
    
    fig2_4, ax2_4 = plt.subplots(figsize=(15, 8))
    plot_Prad_Pg((t_5000_table, t_8000_table), params=plot_params, axis=ax2_4, figure=fig2_4,
                save_path=os.path.join(results_dir,"Prad_Pg_lgTauR.pdf"))
    
    ### 3) Plot T frente a lgTauR comparando con cuerpo gris
    #########################################################################################
    fig3, ax3 = plt.subplots(figsize=(15, 8))
    plot_T((t_5000_table, t_8000_table), params=plot_params, axis=ax3, figure=fig3, grey_atmos=True,
                save_path=os.path.join(results_dir,"T_lgTauR_gris.pdf"))


    ### 4) Poblaciones en tau=0.5 y tau=5
    #########################################################################################

    # Líneas de las tablas en las que tau=0.5 y tau=5:
    tau_05 = np.abs(10**t_5000_table["lgTauR"] - 0.5).argmin() 
    tau_5 = np.abs(10**t_5000_table["lgTauR"] - 5).argmin() 

    # Tablas en las que guardar las poblaciones calculadas:
    poblaciones_t_5000 = QTable(
                                names=("Ne", "HI", "HII", "Hmenos", "HI_n1", "HI_n2", "HI_n3"),
                                units=[u.cm**-3] * 7 ) # Set all columns to have units of m**-3   
    poblaciones_t_8000 = QTable(
                                names=("Ne", "HI", "HII", "Hmenos", "HI_n1", "HI_n2", "HI_n3"),
                                units=[u.cm**-3] * 7)  # Set all columns to have units of m**-3 
    for row in t_5000_table:
         T = row["T"]
         Pe = row["Pe"]
         poblaciones_t_5000.add_row(poblaciones(Pe, T))
    poblaciones_t_5000.add_column(t_5000_table["lgTauR"], index=0)
    poblaciones_t_5000[[tau_05, tau_5]].write(os.path.join(results_dir,"poblaciones_5000.dat"), format="ascii.fixed_width", overwrite=True)

    for row in t_8000_table:
         T = row["T"]
         Pe = row["Pe"]
         poblaciones_t_8000.add_row(poblaciones(Pe, T))
    poblaciones_t_8000.add_column(t_8000_table["lgTauR"], index=0)
    poblaciones_t_8000[[tau_05, tau_5]].write(os.path.join(results_dir,"poblaciones_8000.dat"), format="ascii.fixed_width", overwrite=True)

