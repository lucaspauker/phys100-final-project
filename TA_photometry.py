'''
photometry.py

This is a stub file provided for lab 1.3 part I.

COPY THIS FILE TO A WORK DIRECTORY OF YOUR CHOOSING FIRST!!!!!!

Fill in the missing blanks in the following functions.
'''

###############################################

import numpy as np
import scipy.stats as stats
import scipy.optimize as opt
import astropy.io.fits as fits


################################################
def neglogprob(lambd, data):
    return -np.sum(stats.poisson.logpmf(data, lambd))


def aper_flux(image, x, y, source_radius=10, background_width=10,
             gain=3.0, nimages=1, errFlag=False,read_noise=None,
             avg_dark_current=None, exptime=None):
    """ Measures the flux (in adc counts), gaussian flux uncertainty,
    and the probability that the background is described by a constant

    Args:
        image (np.array): A numpy array containing the image
        x (int): The x position of the object, in pixels
        y (int): The y position of the object, in pixels
        source_radius (float): The radius of a circular region
            centered on the object representing the source.
        background_width (float): The width of the anulus around the
            source radius in which the background will be calculated
        gain (float): The gain of the detector
        nimages (float): The number of images that were coadded to make
            the image
        errFlag (bool): If true, then calculate the errors on the
            magnitude measurement. For use in lab 1.5
        read_noise (float): The average read noise for the CCD
        avg_dark_current (float): The average dark current in counts / pixel /s i.e. ADU / pixel / s
        exptime (float): The exposure time of each image in the coadd

    Returns:
        (float,float,float): 3 parameters, the flux, the 1 sigma error on the
            flux, and the probability that the source is background.

    """

    ysize, xsize = image.shape

    #create 2D arrays with the X,Y coordinates of every pixel
    X,Y = np.meshgrid(np.arange(xsize), np.arange(ysize))

    #calculate the distance of every pixel from the object
    #note: the x-1 accounts for the offset between sextractor and array indicies
    dR = np.sqrt((X - x - 1)**2 + (Y - y - 1)**2)

    #2D boolean arrays selecting the background and source regions
    inSource = dR <= source_radius
    inBackground = np.logical_and(dR > source_radius,
                                  dR <= (background_width + source_radius))

    #counting pixels
    nsourcepix = len(image[inSource])
    nbackgroundpix = len(image[inBackground])

    ### TODO: FINISH THESE LINES
    # calculate the flux of the source, the flux uncertainty, and the
    # probability that the background is described by a constant
    #
    # Feel free to add additional intermediate steps as needed. We will
    # need to calculate
    # the flux with Lab 1.3, but the uncertainties on the flux will wait
    # until Lab 1.5.  Set the fluxerr and background_prob to zero until
    # then.

    flux = np.sum(image[inSource]) - np.median(image[inBackground]) * nsourcepix # counts due to source i.e. ADU due to source
    flux /= exptime # counts / s i.e. ADU / s

    bkg = np.median(image[inBackground]) # counts / pixel i.e. ADU / pixel
    bkg /= exptime # counts / s / pixel i.e. ADU / s / pixel
    if (errFlag):
        # flux = Number of adc counts
        # flux * gain = Number of photoelectrons
        # 1 photoelectron <-> 1 photon
        # photon incidence is Poisson process
        # unc. on flux. This is Gaussian due to CLT because even though #sebastian said this point is wrong
        term1 = flux*gain*exptime # error on photoelectrons/photons # photons i.e. electrons
        term2 = bkg*gain*exptime*nsourcepix # photons i.e. electrons
        term3 = avg_dark_current*gain*exptime*nsourcepix # photons i.e. electrons
        term4 = ((gain*read_noise)**2)*nsourcepix
        # read_noise is in counts / pixel.
        # gain * read_noise is in units of electrons / pixel
        # (gain * read_noise)**2 is in units of electrons^2 / pixel^2
        # term4 is in units of electrons^2 / pixel
        # Sebastian says don't worry about it because term4 is the only one for which process is gaussian
        # terms 1 2 and 3 units dont make sense because we have set
        # are a poisson process that have mean in units of electrons and var in units of electrons
        # then we're assuming that process is approximately gaussian with the same mean and variance,
        # however, this means the standard deviation now has units of sqrt(electrons)
        # this is okay for a poisson process, because the stddev is the same as the mean
        # but for a gaussian this generates an inherent contradiction
        # however it is true that this approximation is valid in this scenario, even if it leads to mixed up units
        fluxerr= np.sqrt(term1 + term2 + term3 + term4) / exptime
        fluxerr /= gain # error on counts
        # the process by which individual measurements are made is Poisson
        # the properly normalized sum of the measurements is Gaussian, thus this uncertainty is Gaussian

        # find best fit uniform value
        # model: background level is uniform
        # if so, our data have poisson distribution with mean value equal to that uniform bkg level
        """
        best_fit = opt.minimize(neglogprob, np.mean(bkg), args=(image[inBackground]/exptime))['x']

        # evaluate chisq
        ndof = nbackgroundpix - 1
        unc = np.sqrt(term2 + term3 +term4) / exptime # error on photoelectrons/photons
        chisq = np.sum((image[inBackground]*gain/exptime - best_fit*gain * np.ones(nbackgroundpix))**2 / unc**2)
        # evaluate reduced chisq
        nu_chisq = chisq / ndof

        # compare reduced chisq to threshold
        p_data = stats.chi2.sf(nu_chisq, ndof) # survival function = 1 - cdf of the chisq distribution with ndof degrees of freedom
        background_prob= p_data # probability that background is described by a constant
        """
        background_prob = None
    else:
        fluxerr = None
        background_prob = None


    return flux, fluxerr, background_prob

#############################################

def aper_mag(flux, fluxerr, errflag=False):
    """ Convert a measured flux and fluxerr into a magnitude and
    magnitude error. Assumes gaussian error propogation.

    Args:
        flux (float or np.array): The flux of the object or objects
        fluxerr (float or np.array): The flux error of the object or
            objects

    Returns
        (flot or np.array, float or np.array): The magnitude and
            magnitude error of the object or objects.

    """

    ### TODO: FINISH THESE LINES
    # calculate the magnitude and the gaussian magnitude uncertainty.
    #
    # Feel free to add additional intermediate steps as needed. Raw
    # magnitudes will be calculated starting in Lab 1.3, and errors will
    # be propagated in Lab 1.5. Set the magerr term to zero until then.

    mag = -2.5 * np.log10(flux)
    constant = np.log10(np.exp(1)) # not -2.5 * np.log(10)
    magerr = np.abs(-2.5 * constant * fluxerr / flux)

    return mag, magerr


###############################################

def create_phot_cat(image_file, detectcat, **otherparams):
    """ Given an image and a detection catalog, calculates the flux at every
    x,y position in the detection catalog and returns the results as a new
    catalog.

    Args:
        image_file (str): The path to the image file to pull the image from
        detectcat (np.array): A detection catalog containing the x and y
            positions on the sources
        otherparams (dict): A dcitionary containing any additional params
            to be passed into aper_flux.

    Returns:
        (np.array): A catalog with flux, fluxerr, mag, magerr, and
            background_prob for each source in the detectcat catalog.

    """

    ##### TODO:  Complete this function to do aperture photometry on the input catalog
    hdulist = fits.open(image_file)
    image = hdulist[0]
    exptime = image.header['EXPTIME']

    #define some arrays to hold our calculations
    fluxs = np.zeros(len(detectcat))
    fluxerrs = np.zeros(len(detectcat))
    backprobs = np.zeros(len(detectcat))
    mags = np.zeros(len(detectcat))
    magerrs = np.zeros(len(detectcat))

    #loop over the objects to measure
    for i, x,y in zip(np.arange(len(detectcat)),
                      detectcat['x'],  #detectcat['X_IMAGE'],
                      detectcat['y']): #detectcat['Y_IMAGE']):

        #calculate the flux; propagate parameters to the aperFlux function
        fluxs[i], fluxerrs[i], backprobs[i] = \
            aper_flux(image.data, x,y, exptime=exptime, **otherparams)

        mags[i], magerrs[i] = aper_mag(fluxs[i], fluxerrs[i])

    #create columns for a new catalog
    nrows = len(detectcat)
    photcat = np.zeros(nrows, dtype=[('flux','f8'),('fluxerr','f8'),
                                     ('backprob','f8'),
                                     ('mag','f8'),('magerr','f8')])
    photcat['flux'] = fluxs
    photcat['fluxerr'] = fluxerrs
    photcat['backprob'] = backprobs
    photcat['mag'] = mags
    photcat['magerr'] = magerrs

    # now write some data
    #fits = fitsio.FITS('test.fits','rw')
    #hdr = fitsio.FITSHDR(detectcat.hdu.header)
    # create a new table extension and write the data
    #fits.write(data, header=hdr)

    return photcat


##############
