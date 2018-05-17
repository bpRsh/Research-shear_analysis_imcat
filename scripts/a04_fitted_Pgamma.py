#!python
# -*- coding: utf-8 -*-
#
# Author      : Bhishan Poudel; Physics Graduate Student, Ohio University
# Date        : Aug 24, 2016
# Last update : Jan 18, 2018
#
# Runtime:  30 secs for redshift 1.0 and 216 files.
"""
.. note::

    1. This program creates fitted P gamma values (i.e. galshear_fpg.cat)
       from galshear_cut.cat and 8 other par files.
    2. It will also create shear catalog file.

:Depends: Depends of following files:
  galshear/galshear_cut.cat
  galshear/galshear_*.par  # 8 par files for monochromatic and chromatic

  i.e. ::

    galshear_c9pg0.par  galshear_cpg0.par   galshear_m9pg0.par  galshear_mpg0.par
    galshear_c9pg1.par  galshear_cpg1.par   galshear_m9pg1.par  galshear_mpg1.par

    Also, galshear_big_cat and galshear_cut.cat

:Output: The Outputs are:
  galshear/galshear_fpg.cat
  galshear/galshear_shear.cat


:Runtime:  2.5 mins (for 4 redshifts and 100 files)
"""
# Imports
import os,sys
import time


def fitted_Pgamma(z):
    """Create fitted Pgamma cat file for given redshift.

    Args:
      z (float): redshift e.g. 0.5, 0.7, 1.0, 1.5

    :Inputs: Inputs are galshear_cut.cat and following 8 Pgamma par files:
        1. galshear_mpg0.par
        2. galshear_m9pg0.par
        3. galshear_cpg0.par
        4. galshear_c9pg0.par
        5. galshear_mpg1.par
        6. galshear_m9pg1.par
        7. galshear_cpg1.par
        8. galshear_c9pg1.par

    :Outputs: galshear_fpg.cat

    :Commands: gen2Dpolymodel galshear_mpg0.par 
    :Commands: lc +all 'x= ox' > galshear_fpg.cat AND OTHERS!!

    """
    # Variables
    pwd = "galshear/galshear_cat_z{0}".format(z)

    # Commands to run
    commands = """
    cd {0}
    lc +all 'ox = %x' 'x = %rg %e[0] 2 vector' < galshear_cut.cat | gen2Dpolymodel galshear_mpg0.par | gen2Dpolymodel galshear_m9pg0.par | gen2Dpolymodel galshear_cpg0.par | gen2Dpolymodel galshear_c9pg0.par | lc +all 'x = %rg %e[1] 2 vector' | gen2Dpolymodel galshear_mpg1.par | gen2Dpolymodel galshear_m9pg1.par | gen2Dpolymodel galshear_cpg1.par | gen2Dpolymodel galshear_c9pg1.par | lc +all 'x = %ox' > galshear_fpg.cat
    """.format(pwd)

    print("\nCreating fitted Pgamma cat file galshear_fpg.cat for redshift {0}:\n".format(z))
    print(commands)
    os.system(commands)




def shear(z):
    """Create shear cat file from fitted Pgamma file.

    Args:
      z (float): redshift e.g. 0.5, 0.7, 1.0, 1.5

    :Inputs: galshear_fpg.cat

    :Outputs: galshear_shear.cat

    :Commands: lc +all 'mg = %me[0] %mPg0mod / %me[1] %mPg1mod / 2 vector'

    """
    # Variables
    pwd = "galshear/galshear_cat_z{0}".format(z)

    # Commands to run
    commands = """
    cd {0}
    lc +all 'mg = %me[0] %mPg0mod / %me[1] %mPg1mod / 2 vector' 'm9g = %m9e[0] %m9Pg0mod / %m9e[1] %m9Pg1mod / 2 vector' 'cg = %ce[0] %cPg0mod / %ce[1] %cPg1mod / 2 vector' 'c9g = %c9e[0] %c9Pg0mod / %c9e[1] %c9Pg1mod / 2 vector' < galshear_fpg.cat | lc +all 'mg_avg = %mg %m9g vadd 0.5 vscale' 'cg_avg = %cg %c9g vadd 0.5 vscale' | lc -i '%dx %dx dot sqrt 0.078 < %mag 3 < and' > galshear_shear.cat
    """.format(pwd)

    print("\nCreating fitted Pgamma cat file galshear_shear.cat for redshift {0}:\n".format(z))
    print(commands)
    os.system(commands)



def main():
    """Run main function."""
    z = float(sys.argv[1])
    fitted_Pgamma(z)
    shear(z)


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
