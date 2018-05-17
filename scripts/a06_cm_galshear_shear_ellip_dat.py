#!python
# -*- coding: utf-8 -*-
#
# Author      : Bhishan Poudel; Physics PhD Student, Ohio University
# Date        : Sep 15, 2016
# Last update : Sep 11, 2017 Mon
"""
.. note::

   This program creates a new dat file with variables
   shear and ellipticity
   by reading only some columns from two input files for
   chromatic and monochromatic cases.

:Depends: This program depends on following:
    color_galshear_shear.dat
    mono_galshear_shear.dat
    color_galshear_ellip.dat
    mono_galshear_ellip.dat

:Outputs: The outputs are in the folder galshear/galshear_cat_z0.5/:
    color_mono_galshear_shear.dat
    color_mono_galshear_ellip.dat

:Runtime: 0.3 seconds.

"""
# Imports
import numpy as np
import subprocess
import os
import sys
import shutil

def cm_shear(z):
    """Create cm shear file with variables:
      r       rkappa     ngals   etmono eterrormono etcolor  eterrorcolor

    Args:
      z (float): redshift e.g. 0.5, 0.7, 1.0, 1.5

    :Inputs: color_galshear_shear.dat and mono_galshear_shear.dat

    :Outputs:
	  color_mono_galshear_shear.dat
	  color_mono_galshear_shear.cat

    """
    # bin   r   ngals   et  eterror   rkappa   kappa  kappaerror  nu
    # 0     1   2       3   4         5        6      7           8
    infile1 = 'galshear/galshear_cat_z{0}/color_galshear_shear.dat'.format(z)
    infile2 = 'galshear/galshear_cat_z{0}/mono_galshear_shear.dat'.format(z)
    r,rkappa,ngals,etmono ,eterrormono  = np.genfromtxt(infile1,delimiter=None,usecols=(1,5,2,3,4),dtype=float,unpack=True)
    r,rkappa,ngals,etcolor,eterrorcolor = np.genfromtxt(infile2,delimiter=None,usecols=(1,5,2,3,4),dtype=float,unpack=True)

    # write to a file
    outfile = 'galshear/galshear_cat_z{0}/color_mono_galshear_shear.dat'.format(z)
    print('Creating : ', outfile)
    with open(outfile,'w') as f:
        # write header
        header = '# r       rkappa     ngals   etmono eterrormono etcolor  eterrorcolor '
        print(header,file=f)

        # write data
        for i in range(len(r)):
            print(r[i],rkappa[i],ngals[i],etmono[i],eterrormono[i],etcolor[i],eterrorcolor[i],sep='   ', file=f)


    # convert dat to cat
    outfile = 'galshear/galshear_cat_z{0}/color_mono_galshear_shear.cat'.format(z)
    cmd = 'lc -C -n r -n rkappa -n ngals -n etmono -n eterrormono -n etcolor -n eterrorcolor < galshear/galshear_cat_z{0}/color_mono_galshear_shear.dat > galshear/galshear_cat_z{0}/color_mono_galshear_shear.cat'.format(z)
    print('Creating : ', outfile)
    subprocess.call(cmd, shell=True)


def cm_ellip(z):
    """Create cm ellip file with variables:
      r       rkappa     ngals   etmono eterrormono etcolor  eterrorcolor

    Args:
      z (float): redshift e.g. 0.5, 0.7, 1.0, 1.5

    :Inputs: color_galshear_ellip.dat and mono_galshear_ellip.dat

    :Outputs:
	  color_mono_galshear_ellip.dat
	  color_mono_galshear_ellip.cat

    :Commands: The command used is::

      lc -C -n r -n rkappa -n ngals -n etmono -n eterrormono -n etcolor -n eterrorcolor < galshear/galshear_cat_z{0}/color_mono_galshear_ellip.dat > galshear/galshear_cat_z{0}/color_mono_galshear_ellip.cat

    """
    infile1 = 'galshear/galshear_cat_z{0}/color_galshear_ellip.dat'.format(z)
    infile2 = 'galshear/galshear_cat_z{0}/mono_galshear_ellip.dat'.format(z)
    r,rkappa,ngals,etmono,eterrormono = np.genfromtxt(infile1,delimiter=None,usecols=(1,5,2,3,4),dtype=float,unpack=True)
    r,rkappa,ngals,etcolor,eterrorcolor = np.genfromtxt(infile2,delimiter=None,usecols=(1,5,2,3,4),dtype=float,unpack=True)



    # write to a file
    outfile = 'galshear/galshear_cat_z{0}/color_mono_galshear_ellip.dat'.format(z)
    print('Creating : ', outfile)
    with open(outfile,'w') as f:

        # write header
        header = '# r       rkappa     ngals   etmono eterrormono etcolor  eterrorcolor '
        print(header,file=f)

        # write data
        for i in range(len(r)):
            print(r[i],rkappa[i],ngals[i],etmono[i],eterrormono[i],etcolor[i],eterrorcolor[i],sep='   ', file=f)


    # convert dat to cat
    outfile = 'galshear/galshear_cat_z{0}/color_mono_galshear_shear.cat'.format(z)
    print('Creating : ', outfile)
    cmd = 'lc -C -n r -n rkappa -n ngals -n etmono -n eterrormono -n etcolor -n eterrorcolor < galshear/galshear_cat_z{0}/color_mono_galshear_ellip.dat > galshear/galshear_cat_z{0}/color_mono_galshear_ellip.cat'.format(z)
    subprocess.call(cmd, shell=True)

def main():
    """Run main function."""
    z = float(sys.argv[1])
    cm_shear(z)
    cm_ellip(z)


if __name__ == "__main__":
    import time

    # Beginning time
    program_begin_time = time.time()
    begin_ctime        = time.ctime()

    #  Run the main program
    main()

    # # Print the time taken
    # program_end_time = time.time()
    # end_ctime        = time.ctime()
    # seconds          = program_end_time - program_begin_time
    # m, s             = divmod(seconds, 60)
    # h, m             = divmod(m, 60)
    # d, h             = divmod(h, 24)
    # print("\nBegin time: ", begin_ctime)
    # print("End   time: ", end_ctime, "\n")
    # print("Time taken: {0: .0f} days, {1: .0f} hours, \
    #   {2: .0f} minutes, {3: f} seconds.".format(d, h, m, s))
