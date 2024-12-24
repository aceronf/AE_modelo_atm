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


#############################################################################################
def plot_lgtauR(tablas:tuple, axis, figure, params, 
                save_path=None, show=False, im_format="pdf"):

    for i, tabla in enumerate(tablas):
        axis.plot(tabla["Depth"].to(u.km), tabla["lgTauR"], color = params["color"][i], 
                  linestyle = params["linestyle"][i], linewidth = 2, label=params["label"][i])

    axis.set_xlim(min([tabla["Depth"][0].to(u.km).value for tabla in tablas]),
                   max([tabla["Depth"][-1].to(u.km).value for tabla in tablas]))
    
    axis.axvline(0, color="black",linestyle="dashdot", linewidth=1, 
        zorder=0, alpha=0.9)
    axis.axhline(0, color="black",linestyle="dashdot", linewidth=1, 
        zorder=0, alpha=0.9)
    
    axis.set_xlabel(r"$r\ [\mathrm{\unit{\kilo\meter}}]$", fontsize=30)
    axis.set_ylabel(r"$\log_{10}(\tau_R)$", fontsize=30)
    
    axis.tick_params(axis='both', which='major', labelsize=24)
    axis.legend(fontsize=24)
    
    if save_path is not None:
        
        # Guardando los archivos en pdf o png desde la terminal.
        # Por defecto en pdf, para guardar en png usar -png
        if '-png' in sys.argv:
            plt.savefig(save_path+'.png',format='png', dpi=1000, bbox_inches='tight')   
            
        else:
            plt.savefig(save_path+'.pdf',format='pdf', dpi=1000, bbox_inches='tight')   
        
    # Cerrando el plot
    plt.close()
    
    

#############################################################################################
def plot_T(tablas:tuple, axis, figure, params, grey_atmos=False,
                save_path=None, show=False, im_format="pdf"):

    for i, tabla in enumerate(tablas):
        axis.plot(tabla["lgTauR"], tabla["T"], color = params["color"][i], 
                  linestyle = params["linestyle"][i], linewidth = 2, label=params["label"][i])
        if grey_atmos:
            T_grey = (3/4*(params["Teff"][i]**4)*(10**tabla["lgTauR"].value+2/3))**(1/4)
            axis.plot(tabla["lgTauR"], T_grey, color = "black", 
                  linestyle = params["linestyle"][i], linewidth = 2, alpha=0.3)

    axis.set_xlim(min([tabla["lgTauR"][0].value for tabla in tablas]),
                   max([tabla["lgTauR"][-1].value for tabla in tablas]))
    
    axis.axvline(0, color="black",linestyle="dashdot", linewidth=1, 
        zorder=0, alpha=0.9)
    
    axis.set_xlabel(r"$\log_{10}(\tau_{\mathrm{R}})$", fontsize=30)
    axis.set_ylabel(r"$T$ [$\unit{\kelvin}$]", fontsize=30)
    
    axis.tick_params(axis='both', which='major', labelsize=24)
    axis.legend(fontsize=24)

    if show:
        plt.show()
    if save_path is not None:
        
        # Guardando los archivos en pdf o png desde la terminal.
        # Por defecto en pdf, para guardar en png usar -png
        if '-png' in sys.argv:
            plt.savefig(save_path+'.png',format='png', dpi=1000, bbox_inches='tight')   
            
        else:
            plt.savefig(save_path+'.pdf',format='pdf', dpi=1000, bbox_inches='tight')   
        
    # Cerrando el plot
    plt.close()
        
#############################################################################################
def plot_Pe(tablas:tuple, axis, figure, params, 
                save_path=None, show=False, im_format="pdf"):

    for i, tabla in enumerate(tablas):
        axis.plot(tabla["lgTauR"], tabla["Pe"], color = params["color"][i], 
                  linestyle = params["linestyle"][i], linewidth = 2, label=params["label"][i])

    axis.set_xlim(min([tabla["lgTauR"][0].value for tabla in tablas]),
                   max([tabla["lgTauR"][-1].value for tabla in tablas]))
    
    axis.axvline(0, color="black",linestyle="dashdot", linewidth=1, 
        zorder=0, alpha=0.9)
    
    axis.set_xlabel(r"$\log_{10}(\tau_{\mathrm{R}})$", fontsize=30)
    axis.set_ylabel(r"$P_{\mathrm{e}}$ [$\unit{\dyn\per\centi\meter\squared}$]", fontsize=30)
    
    axis.tick_params(axis='both', which='major', labelsize=24)
    axis.set_yscale("log")
    axis.legend(fontsize=24)

    if show:
        plt.show()
    if save_path is not None:
        figure.savefig(save_path, format=im_format, bbox_inches='tight')    

#############################################################################################
def plot_Pe_Pg(tablas:tuple, axis, figure, params, 
                save_path=None, show=False, im_format="pdf"):

    for i, tabla in enumerate(tablas):
        axis.plot(tabla["lgTauR"], tabla["Pe"]/tabla["Pg"], color = params["color"][i], 
                  linestyle = params["linestyle"][i], linewidth = 2, label=params["label"][i])

    axis.set_xlim(min([tabla["lgTauR"][0].value for tabla in tablas]),
                   max([tabla["lgTauR"][-1].value for tabla in tablas]))
    
    axis.axvline(0, color="black",linestyle="dashdot", linewidth=1, 
        zorder=0, alpha=0.9)
    
    axis.set_xlabel(r"$\log_{10}(\tau_{\mathrm{R}})$", fontsize=30)
    axis.set_ylabel(r"$P_{\mathrm{e}}/P_{\mathrm{g}}$", fontsize=30)
    
    axis.tick_params(axis='both', which='major', labelsize=24)
    axis.set_yscale("log")
    axis.legend(fontsize=24)

    if show:
        plt.show()
    if save_path is not None:
        figure.savefig(save_path, format=im_format, bbox_inches='tight')    

#############################################################################################
def plot_Prad_Pg(tablas:tuple, axis, figure, params, 
                save_path=None, show=False, im_format="pdf"):

    for i, tabla in enumerate(tablas):
        axis.plot(tabla["lgTauR"], tabla["Prad"]/tabla["Pg"], color = params["color"][i], 
                  linestyle = params["linestyle"][i], linewidth = 2, label=params["label"][i])

    axis.set_xlim(min([tabla["lgTauR"][0].value for tabla in tablas]),
                   max([tabla["lgTauR"][-1].value for tabla in tablas]))
    
    axis.axvline(0, color="black",linestyle="dashdot", linewidth=1, 
        zorder=0, alpha=0.9)
    
    axis.set_xlabel(r"$\log_{10}(\tau_{\mathrm{R}})$", fontsize=30)
    axis.set_ylabel(r"$P_{\mathrm{rad}}/P_{\mathrm{g}}$", fontsize=30)
    
    axis.tick_params(axis='both', which='major', labelsize=24)
    axis.set_yscale("log")
    axis.legend(fontsize=24)

    if show:
        plt.show()
    if save_path is not None:
        figure.savefig(save_path, format=im_format, bbox_inches='tight')    