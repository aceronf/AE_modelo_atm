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
from astropy.constants import Ryd, h, m_e, m_p, c ,k_B, e
from astropy.io.ascii import read
from astropy import units as u
from astropy.table import QTable

import pandas as pd
pd.set_option('mode.chained_assignment', None)

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

niveles_HI = np.array([1,2,3])

# Coefficientes de Gaunt, sección eficaz y opacidad para el HI
# para sus diferentes niveles
def gaunt_coeff_H(T,wl,niveles):
    
    # Coeficiente de Gaunt para las transiciones ligado-libre (bound-free): lineas
    g_bf = 1 - (0.3456/(wl * Ryd)**(1/3))*(((wl*Ryd)/niveles**2)-(1/2))

    # Coeficiente de Gaunt para las transiciones libre-libre (free-free): continuo
    g_ff = 1 + (0.3456/(wl * Ryd)**(1/3))*((wl*k_B*T/(h*c))+(1/2))

    # Estos valores son adimensionales
    return g_bf, g_ff

def cross_section_H(T,Z,niveles,wl):
    
    g_bf, g_ff = gaunt_coeff_H(T=T,wl=wl,niveles=niveles_HI)
    
    # Frecuencia calculada a partir de la longitud de onda
    # con unidades de 1 / s
    freq = c / wl
    
    # No sé de donde salen las unidades de sigma
    # Suponiendo qeu se retorna en cm^2 con las constantes
    sigma_bf = ((2.815e29 * (Z**4 / (niveles**5 * freq ** 3)) * g_bf))*u.cm**2
    
    sigma_ff = (((3.7e8 * (Z**2 / (T**(1/2) * freq**3)) * g_ff)).value)*u.cm**2
    
    return sigma_bf, sigma_ff

def opacity_H(T,Z,Ne,niv_pop,ion_pop,wl,niveles):
    
    sigma_bf, sigma_ff = cross_section_H(T=T,
                                         Z=Z,
                                         wl=wl,
                                         niveles=niveles_HI)
    
    # Esto sale en 1 / s
    freq = (c / wl)
    
    # Esto sale adimensional
    exponent = (-(h * freq)/(k_B*T))
    
    # Pasamos a cgs porque las sigmas estan en cgs
    
    # Esto queda en s3 / cm
    k_bf = (sigma_bf*niv_pop.cgs*(1-np.exp(exponent)))
    
    # Esto queda en 1 / cm4
    k_ff = (sigma_ff*Ne.cgs*ion_pop.cgs*(1-np.exp(exponent)))    
    
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
    
    # Longitud de onda tiene que estar en angstroms
    wl_a = wl.to(u.AA).value
    
    # Suponiendo qeu se retorna en cm^2 con las constantes
    sigma_bf = ((a0 + a1*wl_a + 
                a2*wl_a**2 + 
                a3*wl_a**3 + 
                a4*wl_a**4 + 
                a5*wl_a**5 + 
                a6*wl_a**6)*10e-18)*u.cm**2
    
    f0 = -2.2763 -1.6850*np.log10(wl_a) +0.76661*np.log10(wl_a)**2 -0.053346*np.log10(wl_a)**3
    f1 = +15.2827 -9.2846*np.log10(wl_a) +1.99381*np.log10(wl_a)**2 -0.142631*np.log10(wl_a)**3
    f2 = -197.789 +190.226*np.log10(wl_a) -67.9775*np.log10(wl_a)**2 +10.6913*np.log10(wl_a)**3 -0.625151*np.log10(wl_a)**4
    
    # Adimensional para los logaritmos
    theta = (5040 / T.value)
    
    # Suponiendo qeu se retorna en cm^2 con las constantes
    sigma_ff = 10e-26*10**(f0+f1*np.log10(theta)+f2*np.log(theta)**2)*u.cm**2
    
    
    return sigma_bf, sigma_ff

def opacity_Hmenos(Pe,T,wl,HI):
    
    # En cm**2
    sigma_bf, sigma_ff = cross_section_Hmenos(T,wl)
    
    # Adimensional
    theta = 5040 / T.value
    
    # Pasando a sistema cgs sale dyn no 1 / cm
    k_bf = (4.158e-10 * sigma_bf * Pe.cgs * theta**(5/2) * 10**(0.754*theta)).cgs
    
    # Pasando a sistema cgs sale Ba / cm (casi)
    k_ff = (Pe.cgs * sigma_ff.cgs * HI.cgs).cgs
    
    return k_bf, k_ff


# Opacidad de los electrones
def opacity_e(Ne):
    
    # Pasando a cgs sale  1 / cm3
    k_e = 6.25e-25 * Ne.cgs
    
    return k_e

def opacidades(Pe,T,Ne,Z,niv_pop,ion_pop,wl):
    # Todo está en sistema internacional de entrada
    
    # Opacidad para el H
    #  s3 / cm ----- 1 / cm4
    k_bf_H, k_ff_H = opacity_H(T=T,
                               Z=Z,
                               Ne=Ne,
                               niv_pop=niv_pop,
                               ion_pop=ion_pop,
                               wl=wl,
                               niveles=niveles_HI)
    
    # Opacidad para el H-
    #  dyn ----- Ba / cm
    k_bf_Hmenos, k_ff_Hmenos = opacity_Hmenos(Pe,T,wl,niv_pop[0])
    
    # Opacidad de los electrones
    # 1 / cm3
    k_e = opacity_e(Ne)
    
    # s3 / cm
    k_bf_H_sum = np.sum(k_bf_H)
    
    k_tot = 0

    # Sin dimensiones porque sino no me dejaba sumarlas
    k_tot =  np.sum(k_bf_H).value + np.sum(k_ff_H).value + k_bf_Hmenos.value + k_ff_Hmenos.value + k_e.value
    
    return k_bf_H.value, k_ff_H.value, k_bf_Hmenos.value, k_ff_Hmenos.value, k_e.value, k_bf_H_sum.value, k_tot