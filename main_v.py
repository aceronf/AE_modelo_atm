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

import sympy as sp

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from figuras_v import plot_gen
from poblaciones import poblaciones

### Utilizar LaTeX en las figuras:
#############################################################################################
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=20)  

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
            '''


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


# Primera parte de plotear
def plotear():

    plot_params = {
        "label":[r"$T_{\mathrm{eff}}=5000\ \mathrm{\unit{\kelvin}}$", 
                r"$T_{\mathrm{eff}}=8000\ \mathrm{\unit{\kelvin}}$" ],
        "color":["red", "blue"],
        "linestyle":["-", "--"],
        "Teff":[5000, 8000]
                    }
    
    plot_number = 1
    
    ### 1) Plot de lgTauR frente a la profundidad:
    #########################################################################################
    x_column = "Depth" # Opcional pasarlo así o dentro del for
    y_column = "lgTauR" # Opcional pasarlo así o dentro del for
    
    x_data = []
    y_data = []
    
    # Obteniendo los datos a representar
    for table in table_list:
        x_data.append(table[x_column].to(u.km))
        y_data.append(table[y_column])
    
    # Representando los datos
    if '-noplot' not in sys.argv:
        plot_gen(x_data,y_data,
                label_list=plot_params["label"],
                fig_name=os.path.join(results_dir_path,f"{plot_number}_lgTauR_R"),
                x_axis_label = r"$r\ [\mathrm{\unit{\kilo\meter}}]$",
                y_axis_label = r"$\log_{10}(\tau_R)$",
                guide_lines=[True,(0,0)])
    
    # Aumentando el número identificador del plot
    plot_number += 1
    
    # Olvidando variables para evitar errores
    del x_column,y_column,x_data,y_data
    
    
    ### 2) Plot T, Pe, Pe/Pg y Prad/Pg frente a lgTauR
    #########################################################################################
    
    # Valor fijo para x
    x_column = "lgTauR"
    
    # Magnitudes a plotear en y
    mag_list = ['T','Pe','Pe/Pg','Prad/Pg']
    mag_list_fig_names = ['T','Pe','Pe_Pg','Prad_Pg']

    # Eje logaritmico para cada magnitud
    mag_log_y = [False,True,True,True]
    mag_yax_lab = [r"$T\ [\unit{\kelvin}$]",
                   r"$P_{\mathrm{e}}\ [\unit{\dyn\cdot cm^{-2}}$]",
                   r"$P_{\mathrm{e}}/P_{\mathrm{g}}$",
                   r"$P_{\mathrm{rad}}/P_{\mathrm{g}}$",
    ]
    
    # Plotear para cada magnitud
    for mag_pos,mag in enumerate(mag_list):
        
        x_data = []
        y_data = []
        
        # Para las dos tablas
        for table in table_list:
            x_data.append(table[x_column])
            
            # Si se tiene un cociente se identifica como tal
            if '/' in mag:
                mag1 = mag.split('/')[0]
                mag2 = mag.split('/')[1]
                
                y_data.append(table[mag1].value/table[mag2].value)
            
            else:
                y_data.append(table[mag])
    
        # Representando los datos
        fig_name = os.path.join(results_dir_path,f"{plot_number}_{mag_list_fig_names[mag_pos]}_{x_column}")
        
        if '-noplot' not in sys.argv:
            plot_gen(x_data,y_data,
                label_list=plot_params["label"],
                fig_name=fig_name,
                x_axis_label = r"$\log_{10}(\tau_R)$",
                y_axis_label= mag_yax_lab[mag_pos],
                y_log_scale=mag_log_y[mag_pos],
                y_ticks_dec=0,
                guide_lines=[True,(0,None)])
        
        # Aumentando el número identificador del plot
        plot_number += 1

    ### 3) Plot T frente a lgTauR comparando con cuerpo gris
    #########################################################################################
    
    x_column = "lgTauR" 
    y_column = "T" 
    
    x_data = []
    y_data = []
    
    # Obteniendo los datos a representar
    for table_pos,table in enumerate(table_list):
        x_data.append(table[x_column])
        y_data.append(table[y_column])
        
        # Calculando para la atmósfera gris
        x_data.append(table[x_column])
        y_data_grey = (3/4*(plot_params["Teff"][table_pos]**4)*(10**table[x_column].value+2/3))**(1/4)
        y_data.append(y_data_grey)

    
    # Representando los datos
    if '-noplot' not in sys.argv:
        plot_gen(x_data,y_data,
                label_list=[plot_params["label"][0]]+['']+[plot_params["label"][1]]+[''],
                fig_name=os.path.join(results_dir_path,f"{plot_number}_T_lgTauR_gris"),
                x_axis_label = r"$\log_{10}(\tau_R)$",
                y_axis_label = r"$T\ [\unit{\kelvin}$]",
                line_color = ['red','grey','dodgerblue','grey'],
                line_style = ['-','-','--','--'],
                guide_lines=[True,(0,None)])
    
    # Aumentando el número identificador del plot
    plot_number += 1
    


# Función de partición
def partition_function(T,g_HI,chi_HI):
    """
    Calcula la función de partición basada en la temperatura y los niveles de energía.
    """
    
    # Tema 3 diapositiva 8
    return sum(g * np.exp(-chi / (k_b * T)) for g, chi in zip(g_HI, chi_HI))

def saha_equation(T, n_e, chi):
    """
    Aplica la ecuación de Saha para calcular la densidad de átomos ionizados.
    
    Args:
        T (float): Temperatura en Kelvin.
        n_e (float): Densidad electrónica (1/m^3).
        ionization_energy (float): Energía de ionización en eV.
        
    Returns:
        float: Densidad de partículas ionizadas.
    """
    
    # https://www.britannica.com/science/Saha-equation
    nj_nj1 = (2/n_e) * np.exp(-chi / (k_b * T)) * (2 * np.pi * mass_e * k_b * T / h_cgs**2)**(3/2)
    
    return nj_nj1

def boltzmann_distribution(T,g_HI,chi_HI):
    """
    Calcula las fracciones de población relativa en los niveles de HI según Boltzmann.
    """
    
    # Tema 3, diapositiva 8
    U = partition_function(T,g_HI,chi_HI)    
    populations = [(g/U) * np.exp(-chi / (k_b * T)) for g, chi in zip(g_HI, chi_HI)]
    
    return populations
    


### Programa principal:
#############################################################################################
if __name__ == "__main__":
    
    ### Ficheros con los modelos:
    #############################################################################################
    t_5000 = "t5000.dat"
    t_8000 = "t8000.dat"
    
    ### Se leen los datos y se convierten a tablas de astropy
    #########################################################################################
    t_5000_table = read_table(t_5000)
    t_8000_table = read_table(t_8000)
    table_list = [t_5000_table,t_8000_table]
    #print(t_5000_table)

    ### Se crea un directorio donde guardar los resultados:
    #########################################################################################
    
    cwd = os.getcwd() # Obtenemos el cwd
    results_dir_name = "resultados"
    results_dir_path = f"{cwd}/{results_dir_name}" # Path asoluto a los resultados
    
    # Creamos la carpeta con los plots:
    if os.path.exists(results_dir_path):
            shutil.rmtree(results_dir_path)    
    os.makedirs(results_dir_path)

    plotear()
    
    ### 4) Poblaciones en tau=0.5 y tau=5
    #########################################################################################

    # Líneas de las tablas en las que tau=0.5 y tau=5:
    tau_05 = np.abs(10**t_5000_table["lgTauR"] - 1).argmin() 
    tau_5 = np.abs(10**t_5000_table["lgTauR"] - 10).argmin() 

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
    column_formats = {col: ".3e" for col in poblaciones_t_5000.colnames[1:]}
    poblaciones_t_5000[[tau_05, tau_5]].write(os.path.join(results_dir_path,"poblaciones_5000.dat"), format="ascii.fixed_width", overwrite=True,
                                            formats=column_formats)


    for row in t_8000_table:
        T = row["T"]
        Pe = row["Pe"]
        poblaciones_t_8000.add_row(poblaciones(Pe, T))
    poblaciones_t_8000.add_column(t_8000_table["lgTauR"], index=0)
    column_formats = {col: ".3e" for col in poblaciones_t_8000.colnames[1:]}
    poblaciones_t_8000[[tau_05, tau_5]].write(os.path.join(results_dir_path,"poblaciones_8000.dat"), format="ascii.fixed_width", overwrite=True,
                                            formats=column_formats)

    
    
    



    
    
    