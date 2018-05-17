#!python
# -*- coding: utf-8 -*-
#
# Author      : Douglas Clowe; Associate Professor, Ohio University
# Author      : Bhishan Poudel; Physics Graduate Student, Ohio University
# Date        : Aug 22, 2016
# Last update : Jan 18, 2018
"""
.. note::

    1.  This program runs imcat commands on all the jedisim output files.
        For example, it reads in the following files:
            a. lsst_z0.5_0.fits
            b. lsst90_z0.5_0.fits
            c. mono_z0.5_0.fits
            d. mono90_z0.5_0.fits

        And, also reads other paramters files:
            a. psf10.par

        Then finally create catalogs for each input galaxies like galshear_z0.5_0.cat
        in the output folder galshear/galshear_cat_z0.5 and so on.


    2.  The parameter files are created using a1_psf10_par.py.


    3.  The psf 10 is chosen as the middle of 21 normalized psf, created by
        Phosim software for narrowband_10.icat and narrowband_10.sed for a given seed.
        The Phosim gives unnormalized psf and we normalize all the psf
        so that the sum of all the pixels in these psfs are equal to that of psf10.fits.

        We also note that, psf10.fits is same for both normalized and unnormalized cases.

:Inputs: The inputs are given below:
    All jedisim_output fitsfiles
    psf10.par

:Outputs: galshear/galshear_cat_zREDSHIFT/galshear_*.cat


:Runtime:  15 mins for 216 galaxies with redshift 1.0 only ( Dec 15, 2017 Pisces)

..note::

    This program will read four folders lsst,lsst_mono,lsst90, and lsst_mono90.

"""
# Imports
import subprocess
import os
import time
import shutil
import re
import sys

# beginning time
program_begin_time = time.time()
begin_ctime        = time.ctime()
print('Begin time: ', begin_ctime)

def galshear_cats(z,start,end,indir):
    """This program will create galaxy catalog files.
    
    It will create output folders if they do not exists previously.
    
    """
    
    # Strip Last '/' char in indir
    if indir[-1] == '/':
        indir = indir[0:-1]
    
    
    # Error check (file existence of psf10.par)
    if not os.path.isfile('psf/psf10.par'):
        print('Error: FILE NOT FOUND psf/psf10.par ')
        sys.exit(1)        
        
    # output catalog file
    outdir   = 'galshear/galshear_cat_z{}'.format(z)
    
    # Do not overwrite outdir
    if os.path.isdir(outdir):
        print('ERROR: Output folder exists already.')
        print(outdir)
        sys.exit(1)
    
    # create outdir if not exist.
    if not os.path.isdir(outdir):
            os.makedirs(outdir)
           
        
    for i in range(start,end+1):

        # output catalog file.
        ofile    = outdir + '/galshear_z{}_{:d}.cat'.format(z,i)

        # chromatic files
        # /Users/poudel/Rsh_out/jedisim_v2_outputs/z0.5/lsst/lsst_z0.5_0.fits
        cfile    = indir + "/z{}/".format(z) + "lsst/lsst_z{}_{:d}.fits".format(z,i)
        c9file   = indir + "/z{}/".format(z) + "lsst90/lsst90_z{}_{:d}.fits".format(z,i)
        cparfile = 'psf/psf10.par' # psf10.par

        # monochromatic files
        mfile    = indir + "/z{}/".format(z) + "lsst_mono/lsst_mono_z{}_{:d}.fits".format(z,i)
        m9file   = indir + "/z{}/".format(z) + "lsst_mono90/lsst_mono90_z{}_{:d}.fits".format(z,i)
        mparfile = 'psf/psf10.par'

        # Error check for four files, lsst,lsst90,lsst_mono,lsst_mono90
        for f in [cfile,c9file,mfile,m9file]:
            if not os.path.isfile(f):
                print('Error: FILE NOT FOUND {} '.format(f))
                sys.exit(1)
                
        # After error check, run the bash commands.            
        # commands to run
        commands = "hfindpeaks " + cfile + " -r 0.5 20 | "                                  + \
        "getsky -Z rg 3 | "                                                                 + \
        "apphot -z 30 -M 30 | "                                                             + \
        "getshapes | "                                                                      + \
        "lc +all 'ox = %x' | "                                                              + \
        "cleancat 5 |  "                                                                    + \
        "apphot -z 30 -M 30 | "                                                             + \
        "getshapes | "                                                                      + \
        "lc +all 'x = %x %d vadd' |  "                                                      + \
        "apphot -z 30 -M 30 | "                                                             + \
        "getshapes | "                                                                      + \
        "lc +all 'x = %x %d vadd' |  "                                                      + \
        "apphot -z 30 -M 30 | "                                                             + \
        "getshapes | "                                                                      + \
        "lc +all 'dx = %x %ox vsub' | "                                                     + \
        "gen2Dpolymodel " + cparfile + " | "                                                + \
        "lc +all 'Pg = %psh %psm %stmod[0] %stmod[1] 2 vector "                             + \
                          "%stmod[2] %stmod[3] 2 vector 2 vector "                          + \
                          "%stmod[4] %stmod[5] 2 vector %stmod[6] "                         + \
                          "%stmod[7] 2 vector 2 vector inverse dot "                        + \
                          "dot msub' 'e = %e %psm %stmod[4] "                               + \
        "%stmod[5] 2 vector %stmod[6] "                                                     + \
        "%stmod[7] 2 vector 2 vector inverse dot "                                          + \
        "%stmod[8] %stmod[9] 2 vector dot vsub' | "                                         + \
        "lc +all 'ce = %e' 'cPg = %Pg' 'cmag = %mag' | "                                    + \
        "apphot -z 30 -M 30 -f " + c9file + " | "                                           + \
        "getshapes -f  "+ c9file + " | "                                                    + \
        "lc +all 'Pg = %psh %psm %stmod[0] %stmod[1] 2 vector %stmod[2] "                   + \
                          "%stmod[3] 2 vector 2 vector %stmod[4] "                          + \
                          "%stmod[5] 2 vector %stmod[6] "                                   + \
                          "%stmod[7] 2 vector 2 vector inverse dot dot "                    + \
                          "msub' 'e = %e %psm %stmod[4] %stmod[5] 2 vector "                + \
                          "%stmod[6] %stmod[7] 2 vector 2 vector inverse dot "              + \
                          "%stmod[8] %stmod[9] 2 vector dot vsub' | "                       + \
        "lc +all 'c9e = %e' 'c9Pg = %Pg' 'c9mag = %mag' | "                                 + \
        "apphot -z 30 -M 30 -f " + mfile + " | "                                            + \
        "getshapes -f "+ mfile + " | "                                                      + \
        "gen2Dpolymodel " + mparfile + " | "                                                + \
        "lc +all 'Pg = %psh %psm %stmod[0] %stmod[1] 2 vector %stmod[2] "                   + \
                          "%stmod[3] 2 vector 2 vector %stmod[4] "                          + \
                          "%stmod[5] 2 vector %stmod[6] "                                   + \
                          "%stmod[7] 2 vector 2 vector inverse dot dot msub' 'e = %e %psm " + \
                          "%stmod[4] %stmod[5] 2 vector %stmod[6] "                         + \
                          "%stmod[7] 2 vector 2 vector inverse dot %stmod[8] "              + \
                          "%stmod[9] 2 vector dot vsub' | "                                 + \
        "lc +all 'me = %e' 'mPg = %Pg' 'mmag = %mag' | "                                    + \
        "apphot -z 30 -M 30 -f " + m9file + "| "                                            + \
        "getshapes -f " + m9file + " | "                                                    + \
        "lc +all 'Pg = %psh %psm %stmod[0] %stmod[1] 2 vector %stmod[2] "                   + \
                          "%stmod[3] 2 vector 2 vector %stmod[4] "                          + \
                          "%stmod[5] 2 vector %stmod[6] "                                   + \
                          "%stmod[7] 2 vector 2 vector inverse dot dot msub' 'e = %e %psm " + \
                          "%stmod[4] %stmod[5] 2 vector %stmod[6] "                         + \
                          "%stmod[7] 2 vector 2 vector inverse dot %stmod[8] "              + \
                          "%stmod[9] 2 vector dot vsub' | "                                 + \
        "lc +all 'm9e = %e' 'm9Pg = %Pg' 'm9mag = %mag' > "                                 + \
        ofile

        # run the program
        if not os.path.isfile(ofile):
            print('\nCreating the cat file :', ofile)
            os.system(commands)

##=============================================================================    
def main():
    """Run main function."""
    z = float(sys.argv[1])
    
    # Need: 
    #     1. indir/z0.7/lsst/lsst_z0.7_0.fits
    #     2. indir/z0.7/lsst90/lsst_z0.7_0.fits
    #     3. indir/z0.7/lsst_mono/lsst_mono_z0.7_0.fits
    #     4. indir/z0.7/lsst_mono90/lsst_mono90_z0.7_0.fits
    
    
    # Input directory should have folder z0.7 (or so)
    # indir  = '/Users/poudel/Rsh_out/jedisim_v3_outputs'  # XXX change    
    # indir = '/Volumes/BPWD1/jedisim_v3_outputs/jout_z0.7_2018_Feb02_12_14/' # XXX
    
    indir = '/Users/poudel/Research/a4_jedisim/jedisim/jedisim_output/jout_z0.7_2018_Feb02_18_08' + '/' # XXX
    start,end = 0,295 # inclusive # XXX change
    
    # After changing above parameters, run this.
    galshear_cats(z,start,end,indir)

if __name__ == "__main__":
    import time

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
    print("\nBegin time: ", begin_ctime)
    print("End   time: ", end_ctime, "\n")
    print("Time taken: {0: .0f} days, {1: .0f} hours, \
      {2: .0f} minutes, {3: f} seconds.".format(d, h, m, s))
