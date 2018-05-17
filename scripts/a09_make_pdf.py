#!python
# -*- coding: utf-8 -*-
#
# Author      : Bhishan Poudel; Physics PhD Student, Ohio University
# Date        : Jan 18, 2018
#  done;
"""
.. note::

   This program creates the plots for final shear analysis.

:Depends: This program depends on following:
    imageMagick command montage
    
    plots/galshear_plots_z1.0/8 ps files.

:Outputs: The outputs are in the folder plots/galshear_plots_z1.0/:
    shear_z1_0.pdf

:Runtime: 2 sec

"""
import os,sys


def create_pdf(z):
    z0,z1 = str(z).split('.')
    montage = r'montage *.ps -tile 2x4 -rotate 90 -geometry 1000x1000+20+20  \
    -title "ellipticity and shear for z = {}"  shear.pdf;'.format(z)
    
    command = """
# cd to plot directory containing ps files
cd plots/galshear_plots_z{z0:.0f}.{z1:.0f}

# convert ps images to pdf
{montage:s}

# rename pdf file
mv shear.pdf shear_z{z0:.0f}_{z1:.0f}.pdf

# open final pdf file
open shear_z{z0:.0f}_{z1:.0f}.pdf
    """.format(z0=float(z0),z1=float(z1),montage=montage)

    cmd = command.strip()
    os.system(cmd)


def main():
    """Run main function."""
    z = float(sys.argv[1])
    create_pdf(z)

if __name__ == "__main__":
    main()
