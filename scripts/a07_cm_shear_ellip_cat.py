#!python
# -*- coding: utf-8 -*-
#
# Author      : Bhishan Poudel; Physics PhD Student, Ohio University
# Date        : Sep 27, 2016
# Last update : Sep 11, 2017 Mon
"""
.. note::

   This program creates the shear and ellip cat files from dat files.

:Depends: This program depends on following:
    color_mono_galshear_shear.dat
    color_mono_galshear_ellip.dat

:Outputs: The outputs are in the folder galshear/galshear_cat_z0.5/:
    color_mono_galshear_shear.cat
    color_mono_galshear_ellip.cat

:Runtime: 0.1 seconds.

"""
# Imports
import os,sys
import time


def cm_shear_ellip(z):
    """Create TWO color-mono cat files for shear and ellip from dat files.

    Args:
      z (float): redshift e.g. 0.5, 0.7, 1.0, 1.5

    :Inputs: The outputs are following:
      color_mono_galshear_shear.dat
      color_mono_galshear_ellip.dat


    :Outputs: The outputs are following:
      color_mono_galshear_shear.cat
      color_mono_galshear_ellip.cat

    """
    # Variables
    pwd = "galshear/galshear_cat_z{0}".format(z)

    # Commands to run
    commands = """
    cd {0}
    lc -C -n r -n rkappa -n ngals -n gm -n gmerr -n gc -n gcerr < color_mono_galshear_shear.dat > color_mono_galshear_shear.cat ;
    lc -C -n r -n rkappa -n ngals -n em -n emerr -n ec -n ecerr < color_mono_galshear_ellip.dat > color_mono_galshear_ellip.cat
    """.format(pwd)

    print("\nCreating galshear/color_mono_galshear_shear.cat for redshift {} ".format(z))
    print("Creating galshear/color_mono_galshear_ellip.cat for redshift {} ".format(z))

    # print(commands)
    os.system(commands)


def main():
    """Run main function."""
    z = float(sys.argv[1])
    cm_shear_ellip(z)


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
