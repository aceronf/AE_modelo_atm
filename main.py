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
warnings.filterwarnings("ignore")


import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from figuras import plot_gen
from poblaciones import poblaciones, energias
from opacidades import opacidades

### Utilizar LaTeX en las figuras:
#############################################################################################
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=26)  

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

def search_line(mag,value,tabla):
    
    # Nos devuelve la linea que tiene el valor más cercano a
    # al valor de la magnitud buscada en una tabla
    line = np.abs(10**tabla[mag] - value).argmin() 
    
    return line

# Primera parte de plotear
plot_params = {
        "label":[r"$T_{\mathrm{eff}}=5000\ \mathrm{\unit{\kelvin}}$", 
                r"$T_{\mathrm{eff}}=8000\ \mathrm{\unit{\kelvin}}$" ],
        "color":["red", "blue"],
        "linestyle":["-", "--"],
        "Teff":[5000, 8000]
                    }
def plotear():
    
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
                x_ticks_dec = 0,
                y_axis_label = r"$\log(\tau_\mathrm{R})$",
                aspect_ratio=[1,3],
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
    
    leg_pos_list = ['upper left',
                    'upper left',
                    'upper left',
                    'upper right']
    
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
                x_axis_label = r"$\log(\tau_\mathrm{R})$",
                y_axis_label= mag_yax_lab[mag_pos],
                y_log_scale=mag_log_y[mag_pos],
                y_ticks_dec=0,
                legend_pos=leg_pos_list[mag_pos],
                aspect_ratio=[1,3],
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
                x_axis_label = r"$\log(\tau_\mathrm{R})$",
                y_axis_label = r"$T\ [\unit{\kelvin}$]",
                y_ticks_dec = 0,
                line_color = ['red','grey','dodgerblue','grey'],
                line_style = ['-','-','--','--'],
                aspect_ratio=[1,3],
                guide_lines=[True,(0,None)])
    
    # Aumentando el número identificador del plot
    plot_number += 1
    
def poblaciones_main():
    
    ### 4) Poblaciones en tau=0.5 y tau=5
    #########################################################################################

    # Líneas de las tablas en las que tau=0.5 y tau=5:
    # Estas tau indican las condiciones de Pe y T que se van a utilizar
    tau_05 = search_line("lgTauR",0.5,t_5000_table)
    tau_5 = search_line("lgTauR",5,t_5000_table)

    # Tablas en las que guardar las poblaciones calculadas:
    poblaciones_t_5000 = QTable(
                                names=("Ne", "HI", "HII", "Hmenos", "HI_n1", "HI_n2", "HI_n3"),
                                units=[u.cm**-3] * 7 ) # Set all columns to have units of m**-3   
    poblaciones_t_8000 = QTable(
                                names=("Ne", "HI", "HII", "Hmenos", "HI_n1", "HI_n2", "HI_n3"),
                                units=[u.cm**-3] * 7)  # Set all columns to have units of m**-3 
    
    print('\nPoblaciones para el modelo de 5000 K')
    for row in t_5000_table:
        T = row["T"]
        Pe = row["Pe"]
        poblaciones_t_5000.add_row(poblaciones(Pe, T))
        
    poblaciones_t_5000.add_column(t_5000_table["lgTauR"], index=0)
    column_formats = {col: ".3e" for col in poblaciones_t_5000.colnames[1:]}
    poblaciones_t_5000[[tau_05, tau_5]].write(os.path.join(results_dir_path,"poblaciones_5000.dat"),
                                              format="ascii.fixed_width", 
                                              overwrite=True,
                                              formats=column_formats)

    print('\nPoblaciones para el modelo de 8000 K')
    for row in t_8000_table:
        T = row["T"]
        Pe = row["Pe"]
        poblaciones_t_8000.add_row(poblaciones(Pe, T))
        
    poblaciones_t_8000.add_column(t_8000_table["lgTauR"], index=0)
    column_formats = {col: ".3e" for col in poblaciones_t_8000.colnames[1:]}
    poblaciones_t_8000[[tau_05, tau_5]].write(os.path.join(results_dir_path,"poblaciones_8000.dat"), 
                                              format="ascii.fixed_width", 
                                              overwrite=True,
                                              formats=column_formats)
    
    

def opacidades_main(tauR=1,
                    lambdas=None,
                    niveles=None,
                    plot_k_tot = True,
                    table_name='opacidades',
                    table_export=True):
    
    ### 5) Calculo de las opacidades a tau=1 para los cantos
    #########################################################################################
    tablas_opacidad = []

    if lambdas == None:
        
        if niveles == None:
            
            print('No es posible calcular las opaciades. Proporcione una de las opciones:')
            print('- Longitudes de onda a estudiar')
            print('- Niveles atómicos en los que calcular')
        
        else:
            
            # Longitudes de onda en la que se hace el estudio    
            # Para las longitudes de onda de ionización
            wl_energias = energias(1,niveles)[0]
            wl_io = ((h*c) / -wl_energias).to(u.AA).value
            
            # Seleccionamos un rango alrededor de la longitud
            # de onda de ionización
            wl_cantos = []
            wl_cantos_rango = 1e-3
            for wl in wl_io:
                wl_cantos.append(wl*(1-wl_cantos_rango))
                #wl_cantos.append(wl.value)
                wl_cantos.append(wl*(1+wl_cantos_rango))
                
            lambdas = (np.array(wl_cantos)*u.AA).to(u.m)
        
    else:
        
        lambdas = ((lambdas)*u.AA).to(u.m)
        
    # Estudiando cada modelo de tmeperaturas
    for model_index, model_table in enumerate(table_list):
        
        # Innformamos por pantalla del proceso actual
        print(f'\nOpacidades para el modelo de Teff={plot_params["Teff"][model_index]} K')
        print(f'Generando tabla: {table_name+"_"+str(plot_params["Teff"][model_index])}')
        
        
        # Creamos una tabla
        table = QTable(names=("wl","k_e", "k_ff_Hmenos", 
                                "k_ff_HI", "k_bf_Hmenos_n1", "k_bf_HI_n1", 
                                "k_bf_HI_n2", "k_bf_HI_n3","k_bf_HI_sum","k_tot"),
                                units=['Angstrom',None, None, 
                                    None, None, None,
                                    None, None, None, None]
                                )
            
        # Opteniendo la Pe a tau=1
        tau_index = search_line("lgTauR",tauR,model_table)
        Pe = (model_table["Pe"][tau_index]).si
        # Obteniendo la T a tau=1
        T = (model_table["T"][tau_index]).si
            
        # Calculando poblaciones para cada temperatura
        # Todas las poblaciones en 1/m**3
        Ne, HI, HII, Hmenos, HI_n1, HI_n2, HI_n3 = poblaciones(Pe, T)

        # Agrupando poblaciones de niveles de HI en un array
        niv_pop = u.Quantity([HI_n1, HI_n2, HI_n3])
            
        # Para cada longitud de onda se calculan las opacidades
        # Longitudes de onda en metros
        for wl in lambdas:

            # Calculando las opacidades
            # Todas las magnitudes esán en el sistema internacional
            k_bf_H, k_ff_H, k_bf_Hmenos, k_ff_Hmenos, k_e, k_bf_H_sum, k_tot_val = opacidades(Pe=Pe,
                                                                        T=T,Ne=Ne,niv_pop=niv_pop,
                                                                        HI_pop=HI,
                                                                        HII_pop=HII,
                                                                        wl=wl,Z=1)
                            
            # Mostrar los resultados por pantalla si se solicita
            if '--screen-info' in sys.argv:
                    print(f'T={T} y wl={wl:.2E}')
                    print(f'Opacidad H bf: {k_bf_H}')
                    print(f'Opacidad H ff: {k_ff_H}')
                    print(f'Opacidad H- bf: {k_bf_Hmenos:.2E}')
                    print(f'Opacidad H- ff: {k_ff_Hmenos}')
                    print(f'Opacidad debido a los electrones: {k_e:.4f}')
                    print()
              
            # Añadiendo los datos para cada longitud de onda a las tablas
            table.add_row(((wl.to(u.AA)),
                                                        k_e,
                                                        k_ff_Hmenos,
                                                        k_ff_H,
                                                        k_bf_Hmenos,
                                                        k_bf_H[0],
                                                        k_bf_H[1],
                                                        k_bf_H[2],
                                                        k_bf_H_sum,
                                                        k_tot_val))
                    
        if table_export == True:
            # Exportando las tablas para cada temperatura
            column_formats = {col: ".3e" for col in table.colnames[1:]}
            table_filename = f"{table_name}_{str(plot_params['Teff'][model_index])}.dat"
            table.write(os.path.join(results_dir_path,table_filename), 
                                                    format="ascii.fixed_width", 
                                                    overwrite=True,
                                                    formats=column_formats)

        tablas_opacidad.append(table)

        if plot_k_tot == True:
            
            x_data = [table['wl'].value]*8
            y_data = [table['k_tot'], table['k_e'], table['k_ff_Hmenos'], table['k_ff_HI'], 
                      table['k_bf_Hmenos_n1'], table['k_bf_HI_n1'], table['k_bf_HI_n2'], table['k_bf_HI_n3']]

            plot_gen(x_data=x_data,y_data=y_data,
                    label_list=[r"Total",
                                r"Dispersión de $e^-$",
                                r"f-f, H$^-$",
                                r"f-f, HI",
                                r"b-f, H$^-$",
                                r"b-f, HI, $n=1$",
                                r"b-f, HI, $n=2$",
                                r"b-f, HI, $n=3$"],
                    fig_name=os.path.join(results_dir_path,f"opacidad_desglosada_{plot_params['Teff'][model_index]}"),
                    x_axis_label = r"$\lambda \ [\mathrm{\unit{\angstrom}}]$",
                    x_ticks_dec = 0,
                    y_axis_label = r"$\kappa \ [\mathrm{cm^{-1}}]$",
                    legend_pos='best',
                    leg_col=2,
                    guide_lines=[False,(0,0)],
                    #line_color=["black", "gold", "red", "darkgreen", "orange", "navy", "blue", "cyan"],
                    line_color=["black", "navy", "red", "lime", "orange", "gold", "blue", "cyan"],
                    anchuras=[2]+[1.3]*7,
                    zorder=[10]+[1]*7,
                    #line_style=["solid", "dashdot", "dotted", "dashed", "dotted", "dashed", "dashed", "dashed"],
                    line_style=["solid", "dotted", "dashdot", "dashed", "dashdot", "dashed", "dashed", "dashed"],
                    y_log_scale=True,
                    x_log_scale=False,
                    aspect_ratio=[1.4,2])

    # Plot de las opacidades totales de ambas estrellas
    x_data = [table['wl'].value for table in tablas_opacidad]   
    y_data = [table['k_tot'] for table in tablas_opacidad]

    plot_gen(x_data,y_data,
                    label_list=[r"$T_{\mathrm{eff}}$"+ rf"$\ = {plot_params['Teff'][0]}$"+ r"$\ \unit{\kelvin}$",
                                r"$T_{\mathrm{eff}}$"+ rf"$\ = {plot_params['Teff'][1]}$"+ r"$\ \unit{\kelvin}$"],
                    fig_name=os.path.join(results_dir_path,f"opacidad_total"),
                    x_axis_label = r"$\lambda\ [\mathrm{\unit{\angstrom}}]$",
                    x_ticks_dec = 0,
                    y_axis_label = r"$\kappa_{\mathrm{tot}}\ [\mathrm{cm^{-1}}]$",
                    legend_pos='best',
                    guide_lines=[False,(0,0)],
                    y_log_scale=True,
                    x_log_scale=False,
                    aspect_ratio=[1.4,2])

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

    # Informando del proceso que se está llevando a cabo
    print('\nGenerando los gráficos de los modelos\n')
    plotear()
    

    ### 4) Calculo de las poblaciones
    poblaciones_main()
    
    ### 5) Calculo de las opacidades
    opacidades_main(niveles=[1,2,3],
                    tauR=1,
                    table_name='opacidades_tabla2')
                    
    opacidades_main(lambdas=list(np.arange(500,20010,10)),
                    tauR=1,
                    plot_k_tot=True,
                    table_name='opacidades_rango',
                    table_export=True)
                    
    print('\nEl estudio de poblaciones y opacidades ha finalizado')
    print('Las gráficas y valores obtenidos se pueden consultar en\n./resultados')
