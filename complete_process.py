'''
This is the complete version of the Celestial Quadrilaterals' calibration script
'''
# Here are the libraries we need. Note that wherever you see np, that
# stands for Numpy. 'import' loads an external library.

import numpy as np
import fitsio
from astropy.stats.sigma_clipping import sigma_clip


def average_bias(biasfiles):
    """ Produce an average bias image from a list of individual bias exposures.

    Args:
        biasfiles ([str,...]): A list of string specifying the path to the fits
            files with the bias frames.

    Returns:
        (np.array): The average bias frame. Units of counts.
    """

    # Open each file using fitsio. Fill this in!
    biasdata= [fitsio.read(filepath) for filepath in biasfiles]
    biascube= np.array(biasdata)  # Use numpy to turn the biasdata list into an array

    biascube = sigma_clip(biascube, 3)
    averagebias = np.mean(biascube, axis=0) # Use numpy to average over the first axis. What kind
    # of average should we use?

    return averagebias


def average_dark(darkfiles,averagebias):
    """ Produce an average dark image from a list of individual dark exposures.

    Args:
        darkfiles ([str,...]): A list of string specifying the path to the fits
            files with the dark frames.
        averagebias (np.array): The average bias array returned from
            average_bias.

    Returns:
        (np.array): The average dark frame. Units of counts per second.
    """
    # Read in the darks (very similar to the bias version)
    darkdata = [fitsio.read(i) for i in darkfiles]

    # Make a list of exposure times (hint: you'll need to use the header)
    darkexpo = [fitsio.read_header(i)['EXPTIME'] for i in darkfiles]

    # Turn both lists into numpy arrays
    darkcube = np.array(darkdata)
    darkexpocube = np.array(darkexpo)

    darklist=[]
    # The zip command here loops over two lists at the same time. Iterating
    # over dark cube should return individual darks and iterating over
    # darkexpocube should return times.
    for (image,time) in zip(darkcube,darkexpocube):

        # Correct the dark for the bias. How do you do this?
        cleandark = image - averagebias
        # Normalize the dark for the exposure time. How do you do this?
        normdark = cleandark / time

        # We'll append this normalized dark to our list.
        darklist.append(normdark)

    # As before, we're going to want to turn our list into an array
    cleandarkcube = np.array(darklist)

    # How do we want to average these darks?
    cleandarkcube = sigma_clip(cleandarkcube, 5)
    averagedark = np.mean(cleandarkcube, axis=0)
    return averagedark


def average_flat(flatfiles,averagebias,averagedark):
    """ Produce an average dark image from a list of individual dark exposures.

    Args:
        flatfiles ([str,...]): A list of string specifying the path to the fits
            files with the flat frames.
        averagebias (np.array): The average bias array returned from
            average_bias.
        averagedark (np.array): The average dark array returned from
            average_dark.

    Returns:
        (np.array): The average flat frame. No units.
    """
    # As before, read in the flats from the fits
    flatdata=[fitsio.read(i) for i in flatfiles]

    # We'll also need the exposure time for each flat.
    flatexpo=[fitsio.read_header(i)['EXPTIME'] for i in flatfiles]

    # Make them into numpy arrays.
    flatcube=np.array(flatdata)
    flatexpocube=np.array(flatexpo)

    # Now we'll iterate through the flats like before and make our list
    # of cleaned and normalized flats.
    cleanflatlist=[]
    for (image,time) in zip(flatcube,flatexpocube):
        # First we need to correct the flat for the bias and dark. How?
        cleanflat = image - averagebias - (time * averagedark)
        # Now we need to normalize the flat. What normalization is this?
        normflat= cleanflat / np.median(cleanflat)

        # Append this to our list
        cleanflatlist.append(normflat)

    # Turn our list into a cube
    cleancube = np.array(cleanflatlist)
    # What kind of average should we take?
    cleancube = sigma_clip(cleancube, 5)
    averageflat = np.mean(cleancube, axis=0)
    return averageflat


def science_exposure(rawscidata,rawsciheader,averagebias,averagedark,
    averageflat):
    """ Produce an science frame by combining all the calibration steps.

    Args:
        rawscidata (np.array): The raw science data
        rawsciheader (fitsio.header.FITSHDR): A fits header file for the
            raw science image.
        averagebias (np.array): The average bias array returned from
            average_bias.
        averagedark (np.array): The normalized, average dark array returned
            from average_dark.
        averageflat (np.array): The normalized, average flat array returned
            from average_flat.

    Returns:
        (np.array): The science image.
    """

    # Grab the exposure time from the header
    expotime = rawsciheader['EXPTIME']

    # What statistic do we need from the average flat?
    m = np.median(averageflat)
    # Fill in the formula from image calibration from the parts you have
    scienceimage = m * (rawscidata - averagebias - (averagedark * expotime)) / averageflat
    return scienceimage
