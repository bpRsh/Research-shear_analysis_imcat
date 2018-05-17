#!python
# -*- coding: utf-8 -*-
#
# Author      : Bhishan Poudel; Physics PhD Student, Ohio University
# Date        : Sep 27, 2016
# Last update : Jan 18, 2018
"""
.. note::

    This program finds ellipticity profile from galshear/galshear_shear.cat.
    This does not change input file, only creates 4 output files.

:Depends: galshear_shear.cat

:Outputs: The outputs are given below::

  color_galshear_shear.dat
  mono_galshear_shear.dat
  color_galshear_ellip.dat
  mono_galshear_ellip.dat

:etprofile: Short description is given below::

	NAME
		etprofile --- calculates tangential alignment profile

	SYNOPSIS
	  etprofile [option...] < catfile > asciifile
	      -o io jo	# origin about which we do profile (2048, 2048)
	      -d dlnr		# log bin size 0.25
	      -r rmin rmax	# min and max radii (200, 2000)
	      -l lossfactor	# multiply e by 1/ lossfactor
	      -e ename	# name for 2-vector ellipticity (e)
	      -x xname	# name for 2-vector spatial coordinate (x)

:Runtime: 1 minutes for redshift 1.0 and 216 files.

"""
# Imports
import os,sys
import time

def etprofile_(z):
    """Run etprofile on combined shear cat file and create FOUR dat files for c/m ellp/shr.

    Args:
      z (float): redshift e.g. 0.5, 0.7, 1.0, 1.5

    :Inputs: galshear_shear.cat

    :Outputs:
	  color_galshear_shear.dat
	  mono_galshear_shear.dat
	  color_galshear_ellip.dat
	  mono_galshear_ellip.dat

    """
    # Variables
    pwd = "galshear/galshear_cat_z{0}".format(z)

    # Commands to run
    commands = """
    cd {0}
    etprofile -o 1700 1700 -d 0.2 -r 100 1200 -e cg_avg < galshear_shear.cat | lc -O > color_galshear_shear.dat ;
    etprofile -o 1700 1700 -d 0.2 -r 100 1200 -e mg_avg < galshear_shear.cat | lc -O > mono_galshear_shear.dat ;
    lc +all 'ce_avg = %ce %c9e vadd 0.5 vscale' < galshear_shear.cat | etprofile -o 1700 1700 -d 0.2 -r 100 1200 -e ce_avg | lc -O > color_galshear_ellip.dat   ;
    lc +all 'me_avg = %me %m9e vadd 0.5 vscale' < galshear_shear.cat | etprofile -o 1700 1700 -d 0.2 -r 100 1200 -e me_avg | lc -O > mono_galshear_ellip.dat
    cd -
    """.format(pwd)

    print("\nPwd: {} ".format(pwd))
    print("""Creating:
    color_galshear_shear.dat
    mono_galshear_shear.dat
    color_galshear_ellip.dat
    mono_galshear_ellip.dat """)

    # print(commands)
    os.system(commands)

def main():
    """Run main function."""
    z = float(sys.argv[1])
    etprofile_(z)



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
