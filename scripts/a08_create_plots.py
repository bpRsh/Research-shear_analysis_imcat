#!python
# -*- coding: utf-8 -*-
#
# Author      : Bhishan Poudel; Physics PhD Student, Ohio University
# Date        : Sep 27, 2016
# Last update : Sep 11, 2017 Mon
"""
.. note::

   This program creates the pdf of the plots for shear analysis.
   It will first create postscript files using imcat tool plotcat.
   Then will transform ps to pdf using the tool ps2pdf.
   Then it will combine pdfs using ghostscript command.

   Note that plotcat uses pgplot and we should have installed it first befre
   using this program.
   `pgplot website <http://www.astro.caltech.edu/~tjp/pgplot/macos/>`

:Depends: This program depends on following:
    color_mono_galshear_shear.cat
    color_mono_galshear_ellip.cat

:Outputs: The outputs are in the folder plots/galshear_plots_z0.5/:
    There are 6 output plots for ps, pdf and finally only one combined pdf.
    r_gm_shear.ps     r_gm_shear.pdf
    r_gc_shear.ps     r_gc_shear.pdf
    r_erat_shear.ps   r_erat_shear.pdf
    r_em_ellip.ps     r_em_ellip.pdf
    r_erat_ellip.ps   r_erat_ellip.pdf
    r_erat_ellipp.ps  r_erat_ellipp.pdf

    shear_analysis_z0.5.pdf

:Commands: The commands used are::

    plotcat r gm -w 3  -T 'shear analysis for default'  -d 'r_gm_shear.ps/ps'  <  color_mono_galshear_shear.cat
    ps2pdf r_gm_shear.ps r_gm_shear.pdf


"""

# Imports
import time
import os,sys


def plots(z):
    """Create the final pdf of plots for shear analysis of given redshift.

    Args:
      z (float): redshift e.g. 0.5, 0.7, 1.0, 1.5

    :Inputs: This program depends on following:
        color_mono_galshear_shear.cat
        color_mono_galshear_ellip.cat


    :Outputs: The outputs are in the folder plots/galshear_plots_z0.5/:
        There are 6 output plots for ps, pdf and finally only one combined pdf.
        r_gm_shear.ps     r_gm_shear.pdf
        r_gc_shear.ps     r_gc_shear.pdf
        r_grat_shear.ps   r_grat_shear.pdf
        r_em_ellip.ps     r_em_ellip.pdf
        r_erat_ellip.ps   r_erat_ellip.pdf
        r_erat_ellipp.ps  r_erat_ellipp.pdf

        shear_analysis_z0.5.pdf

    """
    pwd = "galshear/galshear_cat_z{}".format(z)

    # input catalogs
    incat1 = '{}/color_mono_galshear_shear.cat'.format(pwd)
    incat2 = '{}/color_mono_galshear_ellip.cat'.format(pwd)

    # Plot path
    plot_path = 'plots/galshear_plots_z{}'.format(z)
    if not os.path.isdir(plot_path):
        os.makedirs(plot_path)


    # for shear (grat is gamma ratio i.e. shear ratio of chromatic to monocr.)
    cmd1 = "plotcat r gm -w 3                       -T 'shear at redshift {}: plotcat r gm -w 3 ' ".format(z)                           + " -d '{}/r_gm_shear.ps/ps' ".format(plot_path)        + " <  " + incat1
    cmd2 = "plotcat r gc -w 3                       -T 'shear at redshift {}: plotcat r gc -w 3 ' ".format(z)                           + " -d '{}/r_gc_shear.ps/ps' ".format(plot_path)        + " <  " + incat1
    cmd3 = "plotcat r 'grat = %gc %gm /' -w 3       -T 'shear at redshift {}: plotcat r grat = %gc %gm / -w 3 ' ".format(z)             + " -d '{}/r_grat_shear.ps/ps' ".format(plot_path)      + " <  " + incat1
    cmd4 = "plotcat r 'grat = %gc %gm - %gm /' -w 3 -T 'shear at redshift {}: plotcat r grat = %gc %gm - %gm / -w 3 ' ".format(z)       + " -d '{}/r_grat_sheardiff.ps/ps' ".format(plot_path)  + " <  " + incat1

    # for ellipticity
    cmd5 = "plotcat r em -w 3                       -T 'ellipticity at redshift {}: plotcat r em -w 3 ' ".format(z)                     + " -d '{}/r_em_ellip.ps/ps' ".format(plot_path)       + " <  " + incat2
    cmd6 = "plotcat r ec -w 3                       -T 'ellipticity at redshift {}: plotcat r ec -w 3 ' ".format(z)                     + " -d '{}/r_ec_ellip.ps/ps' ".format(plot_path)       + " <  " + incat2
    cmd7 = "plotcat r 'erat = %ec %em /' -w 3       -T 'ellipticity at redshift {}: plotcat r erat = %ec %em / -w 3 ' ".format(z)       + " -d '{}/r_erat_ellip.ps/ps' ".format(plot_path)     + " <  " + incat2
    cmd8 = "plotcat r 'erat = %ec %em - %em /' -w 3 -T 'ellipticity at redshift {}: plotcat r erat = %ec %em - %em / -w 3 ' ".format(z) + " -d '{}/r_erat_ellipdiff.ps/ps' ".format(plot_path) + " <  " + incat2


    # cmds to plot
    commands = [cmd1, cmd2, cmd3, cmd4, cmd5, cmd6,cmd7,cmd8]
    for cmd in commands:
        os.system(cmd)

    print('\nOutput directory: {}'.format(plot_path))

def main():
    """Run main function."""
    z = float(sys.argv[1])
    plots(z)


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
