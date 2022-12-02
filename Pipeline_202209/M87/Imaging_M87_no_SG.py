#-------------------------------------------------------------------------------
# Modules
#-------------------------------------------------------------------------------
#import matplotlib
#matplotlib.use('Agg')

import matplotlib.pyplot as plt

import os
import argparse
import ehtim as eh
import numpy as np

import argparse

from ehtim.const_def import *
from ehtim.statistics.dataframes import *
from ehtim.statistics.stats import *
from ehtim.observing.obs_helpers import *
import itertools as it

from astropy.utils import iers
iers.conf.auto_download = False 


#added by yuwei
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

# Load the image and the array


p_size=140
im = eh.image.load_txt_yuwei('../../models/rowan_m87_test.txt', p_size, p_size)

im.xdim=p_size
im.ydim=p_size
im.display()
res2=1.1304955044788069e-10
res3=res2/2
im=im.blur_circ(res3)
im.display()

eht = eh.array.load_txt('../../arrays_new_version_more_sites/EHT_2022_no_SG.txt')


# Look at the image
im.display()



tint_sec = 10
tadv_sec = 600

tstart_hr = 0
tstop_hr = 12
bw_hz = 8e9

#obs = im.observe(eht, tint_sec, tadv_sec, tstart_hr, tstop_hr, bw_hz,timetype='GMST', elevmin=10, sgrscat=False, ampcal=True, phasecal=False)

#obs.save_txt('simulated_M87_no_SG.txt') # exports a text file with the visibilities
#obs.save_uvfits('simulated_M87_no_SG.uvfits') # exports a UVFITS file modeled on template.UVP


#beamparams = obs.fit_beam() # fitted beam parameters (fwhm_maj, fwhm_min, theta) in radians

    

#####changed by yuwei
parser = argparse.ArgumentParser(description="Fiducial eht-imaging script for M87")
parser.add_argument('-i', '--infile',    default="simulated_M87_no_SG.uvfits", help="input UVFITS file")
parser.add_argument('-i2', '--infile2',  default="",          help="optional 2nd input file (different band) for imaging")
parser.add_argument('-o', '--outfile',   default='simulated_out_M87_no_SG.fits',  help='output FITS image')
parser.add_argument('--savepdf',         default='simulated_out_M87_no_SG.pdf',       help='saves image pdf (True or False)?',action='store_true')
parser.add_argument('--imgsum',          default='summary_M87_no_SG.pdf',       help='generate image summary pdf',action='store_true')
args = parser.parse_args()


pdfout1 = 'M87_source' + '.pdf'
#im.display(cbar_unit=['Tb'],label_type='scale',export_pdf=pdfout1)
im.display(label_type='scale',export_pdf=pdfout1)



zbl = im.total_flux() # total flux
prior_fwhm = 40.0*eh.RADPERUAS  # Gaussian prior FWHM (radians)


# constant regularization weights



reg_term  = {'simple' : 10,    # Maximum-Entropy
              'tv'     : 1.0,    # Total Variation
              'tv2'    : 1.0,    # Total Squared Variation
              'l1'     : 0.0,    # L1 sparsity prior
              'cm'     : 20, #added by yuwei
              'flux'   : 100}    # compact flux constraint


data_term = {'amp'    : 1,    # visibility amplitudes
              'cphase' : 10,    # closure phases
              'logcamp': 10}    # log closure amplitudes



obsfile   = args.infile         # Pre-processed observation file
ttype     = 'nfft'              # Type of Fourier transform ('direct', 'nfft', or 'fast')
#npix      = 64                  # Number of pixels across the reconstructed image deleted by yuwei
npix      = 32                  # Number of pixels across the reconstructed image added by yuwei
fov       = im.fovx() # slightly enlarge the field of view
maxit     = 200                 # Maximum number of convergence iterations for imaging
stop      = 1e-4                # Imager stopping criterion
gain_tol  = [0.02,0.2]          # Asymmetric gain tolerance for self-cal; we expect larger values
                                # for unaccounted sensitivity loss
                                # than for unaccounted sensitivity improvement
uv_zblcut = 0.1e9               # uv-distance that separates the inter-site "zero"-baselines
                                # from intra-site baselines
reverse_taper_uas = 5.0         # Finest resolution of reconstructed features


#-------------------------------------------------------------------------------
# Define helper functions
#-------------------------------------------------------------------------------

# Rescale short baselines to excise contributions from extended flux.
# setting zbl < zbl_tot assumes there is an extended constant flux component of zbl_tot-zbl Jy


# repeat imaging with blurring to assure good convergence
def converge(major=3, blur_frac=1.0):
    for repeat in range(major):
        init = imgr.out_last().blur_circ(blur_frac*res)
        imgr.init_next = init
        #imgr.make_image_I(show_updates=False)
        imgr.make_image_I(show_updates=False)
        
#-------------------------------------------------------------------------------
# Prepare the data
#-------------------------------------------------------------------------------

# Load a single uvfits file
if args.infile2 == '':
    # load the uvfits file
    obs = eh.obsdata.load_uvfits(obsfile)

    
    
# Flag out sites in the obs.tarr table with no measurements


obs_orig = obs.copy() # save obs before any further modifications

obs.reorder_tarr_snr()

beamparams_M87_no_SG = obs.fit_beam()

#-------------------------------------------------------------------------------
# Pre-calibrate the data
#-------------------------------------------------------------------------------

obs_sc = obs.copy() # From here on out, don't change obs. Use obs_sc to track gain changes
res    = obs.res()  # The nominal array resolution: 1/(longest baseline)

# Make a Gaussian prior image for maximum entropy regularization
# This Gaussian is also the initial image
gaussprior = eh.image.make_square(obs_sc, npix, fov)
gaussprior = gaussprior.add_gauss(zbl, (prior_fwhm, prior_fwhm, 0, 0, 0))

# To avoid gradient singularities in the first step, add an additional small Gaussians
gaussprior = gaussprior.add_gauss(zbl*1e-3, (prior_fwhm, prior_fwhm, 0, prior_fwhm, prior_fwhm))

# Reverse taper the observation: this enforces a maximum resolution on reconstructed features
if reverse_taper_uas > 0:
    obs_sc = obs_sc.reverse_taper(reverse_taper_uas*eh.RADPERUAS)
    
# Make a copy of the initial data (before any self-calibration but after the taper)
obs_sc_init = obs_sc.copy()

# First  Round of Imaging
#-------------------------
print("Round 1: Imaging with visibility amplitudes and closure quantities...")

# Initialize imaging with a Gaussian image


#added by yuwei
imgr = eh.imager.Imager(obs_sc, gaussprior, prior_im=gaussprior,
                        flux=zbl, data_term=data_term, maxit=maxit,
                        norm_reg=True,
                        reg_term=reg_term, ttype=ttype, cp_uv_min=uv_zblcut, stop=stop)

# Imaging


#added by yuwei
imgr.make_image_I(show_updates=False)
converge()

# Self-calibrate to the previous model (phase-only);
# The solution_interval is 0 to align phases from high and low bands if needed
obs_sc = eh.selfcal(obs_sc, imgr.out_last(), method='phase', ttype=ttype, solution_interval=0.0)

# Second  Round of Imaging
#-------------------------
print("Round 2: Imaging with visibilities and closure quantities...")

# Blur the previous reconstruction to the intrinsic resolution of ~25 uas
init = imgr.out_last().blur_circ(res)



#added by yuwei
data_term_intermediate = {'vis':imgr.dat_terms_last()['amp']*1,
                          'cphase':imgr.dat_terms_last()['cphase']*1,
                          'logcamp':imgr.dat_terms_last()['logcamp']*1}



#added by yuwei
imgr = eh.imager.Imager(obs_sc, init, prior_im=gaussprior, flux=zbl,
                        data_term=data_term_intermediate, maxit=maxit, norm_reg=True,
                        reg_term = reg_term, ttype=ttype,
                        cp_uv_min=uv_zblcut, stop=stop)

# Imaging
imgr.make_image_I(show_updates=False)
converge()

# Self-calibrate to the previous model starting from scratch
# phase for all sites; amplitudes for LMT

# deleted by yuwei
obs_sc = eh.selfcal(obs_sc_init, imgr.out_last(), method='phase', ttype=ttype)
caltab = eh.selfcal(obs_sc, imgr.out_last(), sites=['LM'], method='both', gain_tol=gain_tol,
                    ttype=ttype, caltable=True)
obs_sc = caltab.applycal(obs_sc, interp='nearest',extrapolate=True)


# Third and Fourth Rounds of Imaging
#-----------------------------------
print("Rounds 3+4: Imaging with visibilities and closure quantities...")

# Increase the data weights before imaging again
data_term_final = {'vis':imgr.dat_terms_last()['vis']*5,
                   'cphase':imgr.dat_terms_last()['cphase']*2,
                   'logcamp':imgr.dat_terms_last()['logcamp']*2}

# Repeat imaging twice
for repeat_selfcal in range(2):
    # Blur the previous reconstruction to the intrinsic resolution of ~25 uas
    init = imgr.out_last().blur_circ(res)

    
    #added by yuwei
    imgr = eh.imager.Imager(obs_sc, init, prior_im=gaussprior, flux=zbl,
                            data_term=data_term_final, maxit=maxit, norm_reg=True,
                            reg_term=reg_term, ttype=ttype,
                            cp_uv_min=uv_zblcut, stop=stop)
        
    # Imaging
    imgr.make_image_I(show_updates=False)
    converge()

    # Self-calibrate

    
    #added by yuwei
    caltab = eh.selfcal(obs_sc_init, imgr.out_last(), method='both',
                        sites=['PV','KP','SPART','SMA','JCMT','ALMA','APEX','LMT','SPT','SMT','GLT','Yonsei','NOEMA','AMT'],
                        ttype=ttype, gain_tol=gain_tol, caltable=True)
    obs_sc = caltab.applycal(obs_sc_init, interp='nearest',extrapolate=True)
#-------------------------------------------------------------------------------
# Output the results
#-------------------------------------------------------------------------------

# Save the reconstructed image
im_out = imgr.out_last().copy()

(error_1, im1_pad_1, im2_shift_1) = im.compare_images(im_out, metric='nrmse',
                                                                  psize=None,
                                                                  target_fov=fov,
                                                                  blur_frac=0.0,
                                                                  beamparams=beamparams_M87_no_SG)

print("error")
print(error_1)
# If an inverse taper was used, restore the final image

# to be consistent with the original data
if reverse_taper_uas > 0.0:
    im_out = im_out.blur_circ(reverse_taper_uas*eh.RADPERUAS)

# Save the final image
im_out.save_fits(args.outfile)

# Optionally save a pdf of the final image
if args.savepdf:
    pdfout = os.path.splitext(args.outfile)[0] + '.pdf'
    im_out.display(label_type='scale', export_pdf=pdfout)
    #im_out.display(label_type='scale', beamparams=beamparams_M87_no_SG, export_pdf=pdfout)
