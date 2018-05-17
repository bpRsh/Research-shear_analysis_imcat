#!python
# -*- coding: utf-8 -*-
#
# Author      : Bhishan Poudel; Physics Graduate Student, Ohio University
# Date        : Aug 24, 2016
# Last update : Jan 24, 2018
# Runtime     : 2 mins (redshift 1.0 and 216 files)
#
# Updates
# on jan 24, 2018 changed rg cut to 2.9.
"""
.. note::

    This program creates P_gamma values (Pg0 and Pg1)
    for chromatic and monochromatic fitsfiles. (c,c9,m,m9)

    The suffix 9 is for 90 degree rotated case.

    In total it creates 2*4 = 8 par files inside galshear.
    Also, in the end it creates two more cat files::

      galshear/galshear_big.cat
      galshear/galshear_cut.cat

    This means in total 2*4 = 8 + 2 = 10 files are created.

:Depends:
  galshear/galshear_*.cat

:Outputs:
  inside galshear 8 par files::

    galshear_c9pg0.par  galshear_cpg0.par   galshear_m9pg0.par  galshear_mpg0.par
    galshear_c9pg1.par  galshear_cpg1.par   galshear_m9pg1.par  galshear_mpg1.par

  Also, galshear_big_cat and galshear_cut.cat


:Runtime: 2 mins (redshift 1.0 and 216 files)

"""
# Imports
import os,sys
import time


def bigcat_cutcat(z):
    """Create big cat file and also create modified cut of it.

    :Inputs: All cat files for the given redshift. e.g. galshear_z0.5_0.fits

    :Ouputs: Outputs are following:
      1. galshear_big.cat
      2. galshear_cut.cat

    """

    # Variables
    pwd = "galshear/galshear_cat_z{0}".format(z)

    # Commands to Run
    commands= """
    cd {0}
    catcats galshear_z{1}_[0123456789].cat galshear_z{2}_[123456789].cat galshear_z{3}_[0123456789].cat > galshear_big.cat
    lc -c < galshear_big.cat
    lc -i '%rg 2.9 > %ce %ce dot 1 < and %me %me dot 1 < and %c9e %c9e dot 1 < and %m9e %m9e dot 1 < and %x[0] 20 > %x[0] 3376 < and %x[1] 20 > and %x[1] 3376 < and and %dx %dx dot sqrt 0.078 < and %mag 3 < and' < galshear_big.cat > galshear_cut.cat
    """.format(pwd,z,z,z)

    print("\nRunning catcats to get big cat file for redshift {} :\n".format(z))
    print(commands)
    os.system(commands)



def Pgamma(z,x):
    """Create P_gamma par files for chromatic, monochromatic and their rotated cat files.

    Args:
      pwd (str): present working directory e.g. galshear/galshear_cat_z0.5
      z (float): redshift e.g. 0.5, 0.7, 1.0, 1.5
      x (str): either of letter 'c', 'c9', 'm', or 'm9'  eg. x = 'c'

    :Inputs: Inputs are following:
      1. galshear_cut.cat

    :Outputs: Outputs are following:
      galshear_cpg0.par galshear_cpg1.par
      galshear_c9pg0.par galshear_c9pg1.par
      galshear_mpg0.par galshear_mpg1.par
      galshear_m9pg0.par galshear_m9pg1.par

    """
    # Variables
    pwd = "galshear/galshear_cat_z{0}".format(z)

    # Commands to run
    commands = """
    cd {0}
    lc +all 'x = %rg %{1}e[0] 2 vector' '{1}Pg0 = %{1}Pg[0][0]' < galshear_cut.cat | fit2Dpolymodel2 x 4 1 {1}Pg0 > galshear_{1}pg0.par
    lc +all 'x = %rg %{1}e[1] 2 vector' '{1}Pg1 = %{1}Pg[1][1]' < galshear_cut.cat | fit2Dpolymodel2 x 4 1 {1}Pg1 > galshear_{1}pg1.par
    """.format(pwd,x)

    print("\nCreating Pgamma par files for {0}pg0 and {0}pg1 for redshift {1}:\n".format(x,z))
    print(commands)
    os.system(commands)

def create_pars(z):
    for x in ['c','c9','m','m9']:
        Pgamma(z,x)

def main():
    """Run main function."""
    # First create big_cat and cut_cat
    z = float(sys.argv[1])
    bigcat_cutcat(z)

    # Then create Pgamma par files.
    create_pars(z)

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
