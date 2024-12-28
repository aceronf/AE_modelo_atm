""" 
ATMÓSFERAS ESTELARES - ANÁLISIS DE DOS MODELOS DE ATMÓSFERA

Código para calcular poblaciones

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

niveles_HI = np.array([1,2,3])
#############################################################################################
def energias(Z, niveles):
    """
    Calcula los niveles de energía de un átomo hidrogenoide con número atómico Z

    Parameters
    ----------
    Z : int
        número atómico
    niveles : array
        Números cuánticos principales de los niveles de energía deseados

    Returns
    -------
    energias
        Array con las energías en eV
    """
    
    lambdas = np.zeros(len(niveles))*u.m
    energias = np.zeros(len(niveles))*u.eV
    for i, n in enumerate(niveles):
        lambdas[i] = (Ryd*Z**2/n**2)**(-1)
        energias[i] = -(h*c/lambdas[i]).to(u.eV)
        
    return energias

#############################################################################################
def deg_HI(niveles):
    """
    Calcula la degeneración de los 

    Parameters
    ----------
    niveles : array
        Números cuánticos principales de los niveles de energía deseados

    Returns
    -------
    Array con los valores de la degeneración en cada nivel
    """
    return 2*niveles**2

#############################################################################################
def Boltzmann_HI(T, niveles):

    g_HI = deg_HI(niveles)  # Degeneración de los niveles de HI
    E_HI = energias(1,niveles) # Energías de los niveles de HI
    chi_HI = E_HI - E_HI[0] # Energías de excitación
    # Función de partición U para el hidrógeno neutro:
    U = sum(g * np.exp(-chi / (k_B * T)) for g, chi in zip(g_HI, chi_HI)).value
    # Poblaciones de cada nivel relativas al número total de átomos de HI:
    poblaciones = [((g/U) * np.exp(-chi / (k_B * T))).value for g, chi in zip(g_HI, chi_HI)]

    return U, poblaciones

#############################################################################################
def Saha(T, Ne, U1, U0, I):

    exponent = (-I/(k_B*T)).to(u.dimensionless_unscaled)
    N1_N0 = ((2*np.pi*m_e*k_B*T)/(h**2))**(3/2) * 2 * U1/U0 * 1/Ne * np.exp(exponent)

    return N1_N0.to(u.dimensionless_unscaled).value

#############################################################################################
def poblaciones(Pe, T):
    ### Estados de ionización:
    # Densidad de electrones:
    Ne = (Pe/(k_B*T)).to(u.m**(-3))
    # Función de partición del HI:
    # Con boltzmann sacamos la distribución de niveles
    U_HI, poblaciones_relativas_HI = Boltzmann_HI(T, niveles_HI)
    # Cociente entre HII y HI
    # Con Saha hacemos el más ionizado respecto del menos
    HII_HI = Saha(T, Ne, 1, U_HI, 13.6*u.eV)
    # Cociente entre HI y el ion H-
    HI_Hmenos = Saha(T, Ne, U_HI, 1, 0.755*u.eV)
    # Cociente entre HII y el ion H-
    # El HII tiene 0 electrones y el H- 2 eleectrones
    HII_Hmenos = HII_HI*HI_Hmenos
    # Densidad del ion H- por conservacion de carga:
    Hmenos = Ne/(HII_Hmenos-1)
    # Densidad de HII
    HII = Ne + Hmenos
    # Densidad de HI
    HI = HI_Hmenos*Hmenos

    ### Estados de excitación para HI:
    HI_n1 = HI * poblaciones_relativas_HI[0]
    HI_n2 = HI * poblaciones_relativas_HI[1]
    HI_n3 = HI * poblaciones_relativas_HI[2]
    
    if '--screen-info' in sys.argv:
        
        print(f'Pe={Pe},T={T}')
        
        print('Mediante Saha')
        print(f'HII/HI: {HII_HI:.2f} atomos mas de HII que de HI')
        print(f'HI/H-:  {HI_Hmenos:.2f} atomos mas de HI que de H-')
        print(f'HII/H-: {HII_Hmenos:.2f} atomos mas de HII que de H-')
        
        print(f'Conociendo Ne={Ne:.2E} (gas ideal)')
        print(f'HI total: {HI:.2E} atomos')
        print(f'HII total:{HII:.2E} atomos')
        print(f'H- total: {Hmenos:.2E} atomos')
        
        print('Mediante Boltzmann')
        print(f'Poblaciones relativas HI, n=0, n=1, n=2:\n{poblaciones_relativas_HI}')
        print(f'HI n=1: {HI_n1:.2E} atomos')
        print(f'HI n=2: {HI_n2:.2E} atomos')
        print(f'HI n=2: {HI_n3:.2E} atomos')
        
        print()

    return Ne, HI, HII, Hmenos, HI_n1, HI_n2, HI_n3

