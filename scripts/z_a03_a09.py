#!python
# -*- coding: utf-8 -*-
"""
Run python scripts from a01 to a07.

:Author:  Bhishan Poudel; Physics Graduate Student, Ohio University

:Date: Aug 01, 2016

:Last update: Oct 2, 2017

"""
# Imports
from __future__ import print_function, unicode_literals, division, absolute_import,with_statement
import os,sys,time

# start time
start_time = time.time()
start_ctime = time.ctime()

def main():
    """Run main function."""
    z = float(sys.argv[1])
    os.system('python a03_Pgamma_cat.py {z:.1f} '.format(z=z)) # 2 min
    os.system('python a04_fitted_Pgamma.py {z:.1f} '.format(z=z)) # 30 secs
    os.system('python a05_etprofile_cm_shear_ellip.py {z:.1f} '.format(z=z)) # 1 min.
    os.system('python a06_cm_galshear_shear_ellip_dat.py {z:.1f} '.format(z=z)) # 1 sec.
    os.system('python a07_cm_shear_ellip_cat.py {z:.1f} '.format(z=z)) # 1 seconds.
    os.system('python a08_create_plots.py {z:.1f} '.format(z=z)) # 1 secs.
    os.system('python a09_make_pdf.py {z:.1f} '.format(z=z)) # 1 secs.


if __name__ == "__main__":
        import time, os

        # Beginning time
        program_begin_time = time.time()
        begin_ctime        = time.ctime()

        #  Run the main program
        main()

        # Print the time taken
        program_end_time = time.time()
        end_ctime        = time.ctime()
        seconds          = program_end_time - program_begin_time
        m, s             = divmod(seconds, 60)
        h, m             = divmod(m, 60)
        d, h             = divmod(h, 24)
        print("\n\nBegin time: ", begin_ctime)
        print("End   time: ", end_ctime, "\n")
        print("Time taken: {0: .0f} days, {1: .0f} hours, \
              {2: .0f} minutes, {3: f} seconds.".format(d, h, m, s))
        print("End of Program: {}".format(os.path.basename(__file__)))
        print("\n")
