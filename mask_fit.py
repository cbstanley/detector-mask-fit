#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Performs a custom sigmoidal fit to find the lower edge position (in pixels)
of a motorized mask in front of a neutron detector. The 2D detector array is
1 x 1 m and digitized into 192 x 256 pixels.

The custom fit function is a sigmoidal with sloping baseline.
'''

import os
import numpy as np
import pylab
from scipy.special import expit
from scipy.optimize import curve_fit


def fit_func(x, base, slope, max_val, k0, x0, k1, x1):
    '''
    Fit function using two sigmoidals and slope baseline
    '''
    y = base + slope*x + max_val*(expit(k0*(x-x0)) - expit(k1*(x-x1)))

    return y


def fit_range(mask_pos):
    '''
    Sets the fit range
    '''
    fitL = max(11, int(mask_pos*0.23 - 18))
    fitR = min(244, int(mask_pos*0.25 + 25))

    return fitL, fitR


def mask_fit(mask_pos, run, tube, popt):
    '''
    Fits over detector mask using fit_func()

    '''
    # Get fit range
    fitL, fitR = fit_range(mask_pos)

    # Always set x0, x1 to these initial guesses
    popt[4] = 0.23*mask_pos
    popt[6] = 0.23*mask_pos + 10 + 0.06*tube

    # Load mask
    fname = 'run_' + str(run) + '_tube' + str(tube) \
            + '_mask' + str(mask_pos) + '.txt'

    dtype = np.dtype([('pixel', 'int'), ('counts', 'f8'), ('error', 'f8')])
    mask_scan = np.loadtxt(fname, dtype=dtype, skiprows=2, usecols=(0, 1, 2))

    xdata = np.array(mask_scan['pixel'][fitL:fitR])
    ydata = np.array(mask_scan['counts'][fitL:fitR])
    edata = np.array(mask_scan['error'][fitL:fitR])
    popt, pcov = curve_fit(fit_func, xdata, ydata, p0=popt, sigma=edata)

    return xdata, ydata, edata, popt


def make_plot(mask_pos, xdata, ydata, edata, popt):
    '''
    Generate a plot with fit for a given mask_pos
    '''
    # Get fit range
    fitL, fitR = fit_range(mask_pos)

    x = np.linspace(fitL, fitR, num=300)
    y = fit_func(x, *popt)
    pylab.plot(xdata, ydata, 'o', label='data')
    pylab.errorbar(xdata, ydata, edata)
    pylab.plot(x, y, label='fit', zorder=3)
    pylab.ylim(-10, 950)
    pylab.legend(loc='best')
    pylab.show()


# Go to directory where data files are located
os.chdir('mask_data')

tube = 0

# Loop over all tubes
while tube < 192:
    try:
        # Start a file with this header
        mask_pixel_list = 'tube ' + str(tube) + '\n' \
            + '''mask_pos mask_edge (pixel)    mask_width (pixels)'''

        run_start = 70424

        # First fit mask position = 700 to get initial guesses
        mask_pos = 700

        run = run_start + (1000 - mask_pos) // 10

        popt_guess = [816, 0, -800, 1.2, 0.23*mask_pos, 1.6,
                      0.23*mask_pos + 10 + 0.06*tube]

        xdata, ydata, edata, popt = mask_fit(mask_pos, run, tube, popt_guess)

        # Plot a given tube and mask_pos to assess fit quality
        if tube == 10 and mask_pos == 700:
            make_plot(mask_pos, xdata, ydata, edata, popt)
        else:
            pass

        # Reset mask position to 1000
        mask_pos = 1000

        # Loop over all mask positions
        while mask_pos >= 0:
            try:
                run = run_start + (1000 - mask_pos) // 10

                # Perform fit
                xdata, ydata, edata, popt = mask_fit(mask_pos, run, tube, popt_guess)

                # Calculate mask edge and width
                pixel_edge = min(popt[4], popt[6])
                pixel_edge = '%.1f' % pixel_edge
                pixel_width = abs(popt[6] - popt[4])
                pixel_width = '%.1f' % pixel_width
                mask_pixel_list += '\n' + str(mask_pos) + '\t' \
                    + str(pixel_edge) + '\t' + str(pixel_width)

#                # Plot a given tube and mask_pos to assess fit quality
#                if tube == 10 and mask_pos == 100:
#                    make_plot(mask_pos, xdata, ydata, edata, popt)
#                else:
#                    pass

            except (IOError, RuntimeError):
                pass

            mask_pos -= 10

        # Write file with mask fit results
        f = open('tube' + str(tube) + '_mask_fit.txt', 'w')
        f.write(mask_pixel_list)
        f.close()

    except (IOError, RuntimeError):
        pass

    tube += 1
