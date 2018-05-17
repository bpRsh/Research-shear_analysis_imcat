#!python
# -*- coding: utf-8 -*-#
"""
Calculate and plot reduced shear vs other quantities.

@author: Bhishan Poudel

@date:  Feb 2, 2018

"""
# Imports
from __future__ import print_function, division,with_statement,unicode_literals,absolute_import
import sys
import matplotlib.pyplot as plt
import scipy
import numpy as np


# Physics
c            = 3e5    # Speed of light km/h
H0           = 67.80  # Hubble constant (km/h)/Mpc
Omega_m      = 0.315  # Current mass density parameter
Omega_lambda = 0.685  # Effective density of dark energy

sigma        = 1000   # Lens Velocity dispersion km/h (lens.txt)
z_source     = 0.7    # Source galaxy redshift
z_lens       = 0.3    # Lens redshift
pix_scale    = 0.2    # final pixel scale (arcsecond per pixel)


# Hubble evolution function E(z)
# https://en.wikipedia.org/wiki/Distance_measures_(cosmology)
#
# When omega_r  and omega_kappa = 0
#
# E(z) = sqrt(     omega_m (1+z)**3  + omega_lambda     ) 
def E_z(z,Omega_m,Omega_lambda):
    return np.sqrt(np.reciprocal((Omega_m * (1+z)**3.0)+ Omega_lambda))

# Comoving distance is Hubble distance * integral of inverse E(z)
# d_C(z) =   int_0_to_z  1/E(z)
def comoving_dist(E_z,z1,z2):
    integral_E_z_inv = scipy.integrate.quad(E_z,z1,z2,args=(Omega_m,Omega_lambda))
    d_H = c/H0
    d_C = d_H * integral_E_z_inv[0]
    return d_C

# NOTE: In lensing equation distances D_s, D_d, D_ds all are ang-diam-distances.
# Angular diameter distance of source from observer
# angular_diameter_distance = comoving_distance / (1+z)  when curvature is zero.
D_s  = comoving_dist(E,0,z_source) / (1+z_source)

# angular diameter distance of souce from lens   
D_ds = comoving_dist(E,z_lens,z_source) / (1+z_source)

# kappa_constant
kappa_constant = 206264.8062471 * (2*np.pi * (sigma*sigma) / (c*c) ) * (D_ds/D_s) / pix_scale



# ellipticity
color_ellip = 'color_galshear_ellip.dat'
mono_ellip =  'mono_galshear_ellip.dat'

# shear
color_shear = 'color_galshear_shear.dat'
mono_shear = 'mono_galshear_shear.dat'

# read values from data files
r    = np.genfromtxt(color_ellip, usecols=(1),unpack=True)  # radius
et_c = np.genfromtxt(color_ellip, usecols=(3),unpack=True)  # et for color
et_m = np.genfromtxt(mono_ellip,  usecols=(3),unpack=True)  # et for mono
gt_c = np.genfromtxt(color_shear, usecols=(3),unpack=True)  # gt for color
gt_m = np.genfromtxt(mono_shear,  usecols=(3),unpack=True)  # gt for mono

# ellipticity and shear ratios
erat = et_c / et_m
grat = gt_c / gt_m

# kappa (convergence)
# kappa = 1/xi    sigma^2/2/G    1/Sigma_crit
#
# kappa is inversely proportional to transverse_distance_xi
kappa = kappa_constant / r




# topic: shear and reduced shear (near to the bottom)
# https://en.wikipedia.org/wiki/Gravitational_lensing_formalism
#
# If source is circular with radius R, lensing will make it elliptical
# a = R / (1-kappa-gamma)
# b = R / (1-kappa+gamma)
#
# For weak lensing: gamma ~ kappa
gamma = kappa
g = gamma / ( 1- kappa)

# difference and ratio
e_diff = et_c - et_m
e_rat  = et_c / et_m

g_diff = gt_c - gt_m
g_rat  = gt_c / gt_m

# et and em
et_gt = ['et_c', 'et_m', 'e_diff','e_rat', 'gt_c', 'gt_m', 'g_diff', 'g_rat']

# plot
for y in et_gt:
    fig, ax = plt.subplots(figsize=(12,8))
    ax.plot(g,eval(y),ls='-',label=y, lw=1)
    plt.xlabel('g')
    plt.ylabel(y)
    plt.title('Plot of {} vs reduced shear'.format(y))
    plt.savefig('plots/{}_vs_g.pdf'.format(y))

    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.close()
