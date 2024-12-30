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
from matplotlib.ticker import FormatStrFormatter
import matplotlib.font_manager as font_manager

### Utilizar LaTeX en las figuras:
#############################################################################################


# Creating a general function to plot
def plot_gen(x_data=[],
             y_data=[],
             con_factor=1,
             label_list=[],
             legend_pos = 'upper left',
             leg_col=1,
             fig_title = None,
             fig_name='figure_withnoname',
             x_axis_label = 'X label',
             x_log_scale = False,
             x_invert=False,
             x_lim=None,
             x_ticks_num = 7,
             x_ticks_dec = 2,
             x_secaxix_ticks = True,
             y_axis_label = 'Y label',
             y_log_scale = False,
             y_invert = False,
             y_lim = None,
             y_ticks_num = 7,
             y_ticks_dec = 2,
             y_secaxix_ticks = True,
             line_color = None,
             line_alpha = None,
             line_style = None,
             anchuras = None,
             mark_style = None,
             mark_size = None,
             zorder = None,
             scatter_plot = False,
             point_lab = None,
             guide_lines = [False,(None,None)]):
    
    # Determinando el tamaño de las figuras
    
    # Cuánto miden de alto
    plot_cols = 1
    size_factor_cols = 6
    fig_size_cols = size_factor_cols * plot_cols
    
    # Cuánto miden de ancho
    plot_rows = 1
    size_factor_rows = 6*2
    fig_size_rows = size_factor_rows * plot_rows

    # Tamaño de ls figuras
    fig_size_global = (fig_size_rows,fig_size_cols)

    # Tamaño de la letra de la leyenda
    legend_font_size = 20
    
    # Definiendo algunos colores
    if line_color is None:
        line_color = ['red','dodgerblue','darkorange','lime','gold','deeppink','blueviolet','black']

    # Definiendo alpha:
    if line_alpha is None:
        line_alpha = [1]*len(x_data)
    # zorder
    if zorder is None:
        zorder = [0]*len(x_data)

    # Estilos de linea
    if line_style is None:
        line_style = ['-', '--', '-', '--','-', '--','-', '--','-', '--','-', '--','-', '--']
    
    # Grosor de linea
    if anchuras is None:
        anchuras = [2]*len(x_data)

    # Estilos de marcadores de scatter
    if mark_style is None:
        mark_style = ['.', '<', '>', 'v','*']
    
    if mark_size is None:
        mark_size = [20,10,10,10,10,10]
    
    # Informando en la terminal sobre lo que se está ploteando
    print(f'Creando el grafico {fig_name.split("/")[-1]}')
    
    # Creando la figura como un solo subplot con ejes 'ax'
    fig = plt.figure(figsize=fig_size_global)
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    ax = plt.subplot(plot_rows, plot_cols, 1)
    
    # Si se indica un título se pone. Por defecto desactivado
    if fig_title is not None:
        fig.suptitle(f'{fig_title}')
    
    # Representando los datos
    for data_pos in range(len(x_data)):
        plt.plot(np.array(x_data[data_pos])/con_factor,
                 np.array(y_data[data_pos]),
                 label=f'{label_list[data_pos]}',
                 color=line_color[data_pos],
                 linestyle=line_style[data_pos],
                 linewidth=anchuras[data_pos],
                 zorder=zorder[data_pos],
                 alpha=line_alpha[data_pos])
    
    # Para añadir lineas verticales de referencia
    # Por defecto desactivadas
    if guide_lines[0]==True:
        if guide_lines[1][0] is not None:
            ax.axvline(guide_lines[1][0], color="black",linestyle="dashdot", linewidth=1, 
                       zorder=0, alpha=0.9)
        if guide_lines[1][1] is not None:
            ax.axhline(guide_lines[1][1], color="black",linestyle="dashdot", linewidth=1, 
                       zorder=0, alpha=0.9)
        
        
    ##### Personalizando los ejes    
    
    ### Eje X inferior
    ax.set_xlabel(f'{x_axis_label}',fontsize=25)

    # Si se dan los limites se ponen sino se calculan
    #if x_lim is not None:
        #ax.set_xticks(np.linspace(x_lim[0], x_lim[1], x_ticks_num))
    #else:
        #ax.set_xticks(np.linspace(np.nanmin(x_data)/con_factor, np.nanmax(x_data)/con_factor, x_ticks_num))
    
    ax.set_xmargin(0.05)
    ax.minorticks_on()
    
    # Numero de decimales de las etiquetas de los ticks
    ax.xaxis.set_major_formatter(FormatStrFormatter(f'%.{x_ticks_dec}f'))
    
    # Escala logarítmica si se indica
    if x_log_scale == True:
        ax.set_xscale('log')
        
    # Invertir el eje si se indica
    if x_invert==True:
        ax.invert_xaxis()


    ### Eje X superior
    axxtop = ax.secondary_xaxis('top')
    axxbot_ticks = ax.get_xticks()
    axxtop.set_xticks(axxbot_ticks)
    
    axxtop.xaxis.set_major_formatter(FormatStrFormatter(f'%.{x_ticks_dec}f'))
    
    # Para poner np los labels sobre los ticks, por defecto se ponen
    if x_secaxix_ticks == False:
        axxtop.xaxis.set_major_formatter(FormatStrFormatter(''))

    axxtop.xaxis.set_minor_formatter(plt.NullFormatter())
    axxtop.minorticks_on()

    axxtop.tick_params(axis='x', which='major')
    
    # Escala logarítmica si se indica
    if x_log_scale == True:
        axxtop.set_xscale('log')
        
    # Invertir el eje si se indica
    if x_invert==True:
        axxtop.invert_xaxis()


    ### Eje Y izquierdo
    ax.set_ylabel(f'{y_axis_label}',fontsize=25,labelpad=10)
    
    # Si se dan los limites se ponen sino se calculan
    #if y_lim is not None:
        #ax.set_yticks(np.linspace(y_lim[0], y_lim[1], y_ticks_num))
    #else:
        #ax.set_yticks(np.linspace(np.nanmin(y_data), np.nanmax(y_data), y_ticks_num))
    
    
    ax.set_ymargin(0.1)
    ax.minorticks_on()
    ax.yaxis.set_major_formatter(FormatStrFormatter(f'%.{y_ticks_dec}f'))

    # Escala logarítmica si se indica
    if y_log_scale == True:
        ax.set_yscale('log')
    
    # Invertir el eje si se indica
    if y_invert==True:
        ax.invert_yaxis()


    ### Eje Y derecho
    axyrig = ax.secondary_yaxis('right')
    
    # Obteniendo los ticks
    axylef_ticks = ax.get_yticks()
    axyrig.set_yticks(axylef_ticks)
    
    axyrig.yaxis.set_major_formatter(FormatStrFormatter(f'%.{y_ticks_dec}f'))
    
    # Para no poner los labels sobre los ticks, por defecto se ponen
    if y_secaxix_ticks == False:
        axyrig.yaxis.set_major_formatter(FormatStrFormatter(''))
 
    axyrig.yaxis.set_minor_formatter(plt.NullFormatter())   
    axyrig.minorticks_on()
        
    axyrig.tick_params(axis='y', which='major')
    
    # Escala logarítmica si se indica
    if y_log_scale == True:
        axyrig.set_yscale('log')
    
    # Invertir el eje si se indica
    if y_invert==True:
        axyrig.invert_yaxis()


    #### Otros aspectos a personalizar
    
    # Generando un grid sutil
    ax.grid(True,alpha=0.2)  

    # Propiedades de la leyenda
    font = font_manager.FontProperties(size=legend_font_size)
    plt.legend(loc=legend_pos,ncol=leg_col,prop=font)

    # Dando espacio si se pone el título
    fig.tight_layout(pad=0.5)
        
    # Guardando los archivos en pdf o png desde la terminal.
    # Por defecto en pdf, para guardar en png usar -png en terminal
    if '-png' in sys.argv:
        plt.savefig(fig_name+'.png',format='png', dpi=1000, bbox_inches='tight')   
        
    else:
        plt.savefig(fig_name+'.pdf',format='pdf', dpi=1000, bbox_inches='tight')   
        
    # Cerrando el plot
    plt.close()   