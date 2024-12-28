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

def gaunt_coeff_H(T,wl,niveles):
    
    # Coeficiente de Gaunt para las transiciones ligado-libre (bound-free): lineas
    g_bf = 1 - (0.3456/(wl * Ryd)**(1/3))*((wl*Ryd/niveles**2)-(1/2))
    
    # Coeficiente de Gaunt para las transiciones libre-libre (free-free): continuo
    g_ff = 1 + (0.3456/(wl * Ryd)**(1/3))*((wl*k_B*T/(h*c))+(1/2))

    # Estos valores son adimensionales
    return g_bf, g_ff


# Sección eficaz y opacidad para el HI
# para sus diferentes niveles
def cross_section_H(T,Z,niveles,wl):
    
    g_bf, g_ff = gaunt_coeff_H(T=T,wl=wl,niveles=niveles_HI)
    
    # Frecuencia calculada a partir de la longitud de onda
    # con unidades de 1 / s
    freq = c / wl
    
    # No sé de donde salen las unidades de sigma
    sigma_bf = (2.815e29 * (Z**4 / (niveles**5 * freq ** 3)) * g_bf)*u.cm**2/u.s**3
    
    sigma_ff = 3.7e8 * (Z**2 / T**(1/2) * freq**3) * g_ff
    
    return sigma_bf, sigma_ff

def opacity_H(T,Z,Ne,niv_pop,ion_pop,wl,niveles):
    
    sigma_bf, sigma_ff = cross_section_H(T=T,
                                         Z=Z,
                                         wl=wl,
                                         niveles=niveles_HI)
    
    freq = c / wl
    
    exponent = (-(h * freq)/(k_B*T)).to(u.dimensionless_unscaled)
    
    k_bf = sigma_bf*niv_pop*(1-np.exp(exponent))
    
    k_ff = sigma_ff*Ne*ion_pop*(1-np.exp(exponent))
    
    return k_bf, k_ff


# Sección eficaz y opacidad para el H-
def cross_section_Hmenos(T,wl):
    
    a0 = 1.99654
    a1 = -1.18267e-5
    a2 = 2.64243e-6
    a3 = -4.40524e-10
    a4 = 3.23992e-14
    a5 = -1.395683e-18
    a6 = 2.78701e-23
    
    wl = wl.value
    sigma_bf = (a0 + a1*wl + 
                a2*wl**2 + 
                a3*wl**3 + 
                a4*wl**4 + 
                a5*wl**5 + 
                a6*wl**6)*10e-18
    
    f0 = -2.2763 -1.6850*np.log10(wl) +0.76661*np.log10(wl)**2 -0.053346*np.log10(wl)**3
    f1 = +15.2827 -9.2846*np.log10(wl) +1.99381*np.log10(wl)**2 -0.142631*np.log10(wl)**3
    f2 = -197.789 +190.226*np.log10(wl) -67.9775*np.log10(wl)**2 +10.6913*np.log10(wl)**3 -0.625151*np.log10(wl)**4
    
    theta = (5040 / T.value)
    
    sigma_ff = 10e-26*10**(f0+f1*np.log10(theta)+f2*np.log(theta)**2)
    
    return sigma_bf, sigma_ff

def opacity_Hmenos(Pe,T,wl,HI):
    
    sigma_bf, sigma_ff = cross_section_Hmenos(T,wl)
    
    theta = 5040 / T.value
    
    k_bf = 4.158e-10 * sigma_bf * Pe * theta**(5/2) * 10**(0.754*theta)
    
    k_ff = Pe * sigma_ff * HI
    
    return k_bf, k_ff


# Opacidad de los electrones
def opacity_e(Ne):
    
    k_e = 6.25e-25 * Ne
    
    return k_e

def opacidades(Pe,T,Ne,Z,niv_pop,ion_pop,wl):
    
    # Opacidad para el H
    k_bf_H, k_ff_H = opacity_H(T=T,
                               Z=Z,
                               Ne=Ne,
                               niv_pop=niv_pop,
                               ion_pop=ion_pop,
                               wl=wl,
                               niveles=niveles_HI)
    
    # Opacidad para el H-
    k_bf_Hmenos, k_ff_Hmenos = opacity_Hmenos(Pe,T,wl,niv_pop[0])
    
    # Opacidad de los electrones
    k_e = opacity_e(Ne)
    
    return k_bf_H, k_ff_H, k_bf_Hmenos, k_ff_Hmenos, k_e