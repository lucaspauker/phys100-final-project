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


def multi_ellipse_aper_flux(image, n_ellipses=0, a_list=None, b_list=None,
                            h_list=None, k_list=None, bgx=0, bgy= 0, background_width=10,
                            angles=None, gain=3.0, nimages1=1, nimages2=1, errFlag=False, read_noise1=None,
                            read_noise2=None, avg_dark_current1=None, avg_dark_current2=None, exptime1=None,
                            exptime2=None):
    """ Measures the flux (in adc counts), gaussian flux uncertainty,
    and the probability that the background is described by a constant

    Args:
        image (np.array): A numpy array containing the image, should be coadded and
                    not averaged
        h_list (iterable): A list of x positions of the ellipse centers, in pixels
        k_list (iterable): The y positions of the ellipse centers, in pixels
        background_width (float): The width of the anulus around the
            source radius in which the background will be calculated
        gain (float): The gain of the detector
        nimages1 (float): The number of images of type 1 that were coadded to make
            the image
        nimages2 (float): The number of images of type 2 that were coadded to make
            the image
        errFlag (bool): If true, then calculate the errors on the
            magnitude measurement. For use in lab 1.5
        read_noise1 (float): The average read noise for the CCD for type 1 calibrations
        read_noise2 (float): The average read noise for the CCD for type 2 calibrations
        avg_dark_current1 (float): The average dark current for type 1 calibrations
        avg_dark_current2 (float): The average dark current for type 2 calibrations
        exptime (float): The exposure time of each image in the coadd

    Returns:
        (float,float,float): 3 parameters, the flux, the 1 sigma error on the
            flux, and the probability that the source is background.

    """

    ysize, xsize = image.shape
    inSource_list = []
    exptime = exptime1 + exptime2

    # create 2D arrays with the X,Y coordinates of every pixel
    X,Y = np.meshgrid(np.arange(xsize), np.arange(ysize))

    for i in range(n_ellipses):
        # calculate the distance of every pixel from the object
        # note: the x-1 accounts for the offset between sextractor and array
        # indicies
        dR = ((X - h_list[i] - 1)*np.cos(angles[i]) + (Y - k_list[i] - 1)*np.sin(angles[i]))**2/a_list[i]**2 + \
             (-(X - h_list[i] - 1)*np.sin(angles[i]) + (Y - k_list[i] - 1)*np.cos(angles[i]))**2/b_list[i]**2

        # 2D boolean arrays selecting the background and source regions
        inSource = (dR <= 1)
        for inOtherSource in inSource_list:
            inSource = inSource & np.invert(inOtherSource)
        inSource_list.append(inSource)

    inAnySource = inSource_list[0]
    for inSource in inSource_list[1:]:
        inAnySource = inAnySource | inSource

    dR = (X - bgx - 1)**2 + (Y - bgy - 1)**2
    inBackground = (dR <= background_width**2) & np.invert(inAnySource)


    # calculate the flux of the source, the flux uncertainty, and the
    # probability that the background is described by a constant
    #
    # Feel free to add additional intermediate steps as needed. We will
    # need to calculate
    # the flux with Lab 1.3, but the uncertainties on the flux will wait
    # until Lab 1.5.  Set the fluxerr and background_prob to zero until
    # then.

    nbackgroundpix = np.sum(inBackground)
    B_bar = np.sum(image[inBackground]) / nbackgroundpix

    nsourcepix_list = []
    fluxes = []
    fluxerrs = []

    for inSource in inSource_list:
        # Counting pixels
        nsourcepix = np.sum(inSource) # N_A
        nsourcepix_list.append(nsourcepix)
        #nbackgroundpix = np.sum(inBackground) #N_B
        # Expected background
        #B_bar = np.sum(image[inBackground]) / nbackgroundpix

        flux = ((np.sum(image[inSource]) - (B_bar * nsourcepix)) * gain) / exptime
        fluxes.append(flux)

        if (errFlag):
            if flux < 0:
                fluxerr = 0
            else:
                Fobj = np.sum(image[inSource]) * gain
                Fsky = B_bar * nsourcepix * gain
                D1 = avg_dark_current1 * nsourcepix * gain * exptime1
                D2 = avg_dark_current2 * nsourcepix * gain * exptime2
                RN1 = (read_noise1 ** 2) * nsourcepix * nimages1 * gain
                RN2 = (read_noise2 ** 2) * nsourcepix * nimages2 * gain
                #signoise = (Fobj*t) / np.sqrt(Fobj*t + Fsky*t*n + D*t*n + RN**2*n)
                fluxvar = (Fobj + Fsky + D1 + D2 + RN1 + RN2)  # check
                fluxerr = np.sqrt(fluxvar/(exptime ** 2))
                #fluxerr = signoise
            background_prob = 0
            fluxerrs.append(fluxerr)
        else:
            fluxerr = 0
            background_prob = 0
            fluxerrs.append(fluxerr)

    return fluxes, fluxerrs, inSource_list, inBackground


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
        avg_dark_current (float): The average dark current
        exptime (float): The exposure time of each image in the coadd

    Returns:
        (float,float,float): 3 parameters, the flux, the 1 sigma error on the
            flux, and the probability that the source is background.

    """

    ysize, xsize = image.shape

    # create 2D arrays with the X,Y coordinates of every pixel
    X,Y = np.meshgrid(np.arange(xsize), np.arange(ysize))

    # calculate the distance of every pixel from the object
    # note: the x-1 accounts for the offset between sextractor and array
    # indicies
    dR = np.sqrt((X - x - 1)**2 + (Y - y - 1)**2)

    # 2D boolean arrays selecting the background and source regions
    inSource = dR <= source_radius
    inBackground = np.logical_and(dR > source_radius,
                                  dR <= (background_width + source_radius))

    # calculate the flux of the source, the flux uncertainty, and the
    # probability that the background is described by a constant
    #
    # Feel free to add additional intermediate steps as needed. We will
    # need to calculate
    # the flux with Lab 1.3, but the uncertainties on the flux will wait
    # until Lab 1.5.  Set the fluxerr and background_prob to zero until
    # then.

    # Counting pixels
    nsourcepix = np.sum(inSource) # N_A
    nbackgroundpix = np.sum(inBackground) #N_B
    # Expected background
    B_bar = np.sum(image[inBackground]) / nbackgroundpix
    print(nsourcepix)
    print(nbackgroundpix)
    print((B_bar * nsourcepix * gain) / exptime)

    flux = ((np.sum(image[inSource]) - (B_bar * nsourcepix)) * gain) / exptime

    if (errFlag):
        if flux < 0:
            fluxerr = 0
        else:
            Fobj = flux * gain
            t = exptime
            Fsky = B_bar * nsourcepix * gain
            n = nsourcepix
            D = avg_dark_current * nsourcepix * gain
            RN = read_noise * nsourcepix * gain
            signoise = (Fobj*t) / np.sqrt(Fobj*t + Fsky*t*n + D*t*n + RN**2*n)
            #fluxvar = (flux + (B_bar * gain / exptime) + ((read_noise * nimages * gain) / exptime) + (avg_dark_current * nimages * gain))  # check
            #fluxerr = np.sqrt(fluxvar)
            fluxerr = signoise
        background_prob = 0
    else:
        fluxerr = 0
        background_prob = 0

    return flux, fluxerr, background_prob

#############################################


def aper_mag(flux, fluxerr):
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

    # calculate the magnitude and the gaussian magnitude uncertainty.
    #
    # Feel free to add additional intermediate steps as needed. Raw
    # magnitudes will be calculated starting in Lab 1.3, and errors will
    # be propagated in Lab 1.5. Set the magerr term to zero until then.

    if flux <= 0: return np.nan, 0
    mag = -2.5 * np.log10(flux)
    magerr = np.abs(2.5 * (fluxerr / (flux * np.log(10))))

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

    # TODO:  Complete this function to do aperture photometry on the input
    # catalog
    hdulist = fits.open(image_file)
    image = hdulist[0]
    exptime = image.header['EXPTIME']

    # define some arrays to hold our calculations
    fluxs = np.zeros(len(detectcat))
    fluxerrs = np.zeros(len(detectcat))
    backprobs = np.zeros(len(detectcat))
    mags = np.zeros(len(detectcat))
    magerrs = np.zeros(len(detectcat))

    # loop over the objects to measure
    for i, x,y in zip(np.arange(len(detectcat)),
                      detectcat['x'],
                      detectcat['y']):

        # calculate the flux; propagate parameters to the aperFlux function
        fluxs[i], fluxerrs[i], backprobs[i] = \
            aper_flux(image.data, x,y, exptime=exptime, **otherparams)

        mags[i], magerrs[i] = aper_mag(fluxs[i], fluxerrs[i])

    # create columns for a new catalog
    nrows = len(detectcat)
    photcat = np.zeros(nrows, dtype=[('flux','f8'),('fluxerr','f8'),
                                     ('backprob','f8'),
                                     ('mag','f8'),('magerr','f8')])
    photcat['flux'] = fluxs
    photcat['fluxerr'] = fluxerrs
    photcat['backprob'] = backprobs
    photcat['mag'] = mags
    photcat['magerr'] = magerrs

    return photcat


##############
