import sys, os

# Pre-pending everything. Shhhh nobody will know.
sys.path.insert(0, '/user/phstf/md0046/tools/xpa/build/lib')
sys.path.insert(0, '/user/phstf/md0046/tools/astropy/build/lib.linux-x86_64-2.7')
sys.path.insert(0, '/user/phstf/md0046/tools/imexam/build/lib/python2.7/site-packages/imexam-0.5.2-py2.7-linux-x86_64.egg')
sys.path.insert(0, '/user/phstf/md0046/tools/photutils/build/lib.linux-x86_64-2.7/')

# Importing basic utils
import pyds9
import imexam
import time
import pyfits as fits
import numpy as np
from scipy import optimize
import math
import matplotlib.pyplot as plt
from ctypes import c_longlong as ll
from imexam import math_helper

# Parameter dictionnary for our imexam routine
fit2d_pars = {
    'title': ['', 'Title of the plot'],
    'xlabel': ['Radius', 'X-axis label'], 
    'ylabel': ['Pixel value', 'Y-axis label'],
    'fitplot': [True, 'Overplot profile fit?'],
    'fittype': ['gaussian', 'Profile type to fit'],
    'center': [True, 'Center object in aperture?'],
    'background': [True, 'Fit and subtract background?'],
    'radius': [17., 'Object radius'],
    'buffer': [5., 'Background buffer width'],
    'width': [17., 'Background width'],
    'iterations': [5, 'Number of radius adjustment iterations'],
    'xorder': [0, 'Background x order'],
    'yorder': [0, 'Background y order'],
    'magzero': [25., 'Magnitude zero point'],
    'beta': ['INDEF', 'Moffat beta parameter'],
    'rplot': [20., 'Plotting radius'],
    'x1': ['INDEF', 'X-axis window limit'],
    'x2': ['INDEF', 'X-axis window limit'],
    'y1': ['INDEF', 'Y-axis window limit'],
    'y2': ['INDEF', 'Y-axis window limit'],
    'marker': ['+', 'Marker character?'],
    'szmarker': [1., 'Marker size'],
    'logx': [False, 'log scale x-axis'],
    'logy': [False, 'log scale y-axis'],
    'ticklab': [True, 'Label tick marks'],
    'majrx': ['5', 'Number of major divisions along x grid'],
    'minrx': ['5', 'Number of minor divisions along x grid'],
    'majry': ['5', 'Number of major divisions along y grid'],
    'minry': ['5', 'Number of minor divisions along y grid']}

# Functions not found in any package that are required for the lab
def combine(file_list, ctype='median'):
    """
    Combines the files listed in file_list using a specific method
    
    Arguments:
    file_list -- a string indicating a filename. 
    ctype -- the method to use to combine the images listed in file_list. For the moment, only 'median' is implemented
    """
    f = open(file_list)
    fits_files = []
    for filename in f:
        filename = filename.strip()
        print('Adding {0} to the list of combination'.format(filename))
        fits_files += [fits.getdata(filename)]
    f.close()

    sys.stdout.write('Combining images ... ')
    sys.stdout.flush()
    nr, nc = fits_files[0].shape
    output = np.zeros((nr, nc))

    # SUPER INEFFICIENT ! Must be improved
    ok = False
    if ctype == 'median':
        med = len(fits_files)//2
        for i in xrange(output.size):
            x, y = i%nc, i//nc
            vals = []
            for f in fits_files:
                vals += [f[y, x]]
            vals.sort()
            #print(vals)
            output[y, x] = vals[med]
        ok = True
    else:
        print('ERROR : Other ctypes than \'median\' have not been implemented yet')

    if ok:
        sys.stdout.write('DONE\n')
            
    return output

def export_headers(file_out, folder):
    """
    Exports all of the observation and exposure times present in folder to a file (file_out). 
    """
    fout = open(file_out, 'w')
    for f in os.listdir(folder):
        if f[-5:] == '.fits' and f[0] != '.':
            header = fits.getheader(folder+'/'+f)
            if 'DATE-OBS' not in header or 'EXPTIME' not in header:
                print('Warning : Date of observation or exposition time not found in meta-data of file : ' + folder + '/' + f)
            else:
                fout.write('{0} {1}\n'.format(header['DATE-OBS'], header['EXPTIME']))

    fout.close()


def gauss_p(x, I0, sigma, offset):
    """ Gaussian function on a vector of points (x,y)"""
    return I0 * np.exp(-0.5*(x[:,0]**2+x[:,1]**2)/sigma**2) + offset

def gauss_r(r, I0, sigma, offset):
    """ Gaussian function on a vector of distances (r)"""
    return I0 * np.exp(-0.5*(r/sigma)**2) + offset

def gauss_center_p(x, I0, sigma, offset, x0, y0):
    """ Gaussian function used to fit the center of the point """
    return I0 * np.exp(-0.5*((x0-x[:,0])**2 + (y0-x[:,1])**2)/sigma**2) + offset


def fit_2d_new(self, x, y, data, form=None, subsample=4, fig=None):
    params = self.fit2d_pars

    if not form:
        form = params["fittype"][0]

    if form not in ('gaussian', 'moffat'):
        warnings.warn('Error : The fitting function must be either \'gaussian\' or \'moffat\'')
        return

    center = params["center"][0]
    background = params["background"][0]
    radius = params["radius"][0]
    buffer = params["buffer"][0]
    width = params["width"][0]
    xorder = params["xorder"][0]
    yorder = params["yorder"][0]
    medsky = xorder <= 0 or yorder <= 0
    loc_beta = params["beta"][0]

    magzero = params["magzero"][0]
    rplot = params["rplot"][0]
    fitplot = params["fitplot"][0]
        
    centerx, centery = x, y

    chunk = data[centery-rplot:centery+rplot, centerx-rplot:centerx+rplot]
    xr = np.arange(centerx-rplot, centerx+rplot) - centerx
    yr = np.arange(centery-rplot, centery+rplot) - centery
    X, Y = np.meshgrid(xr, yr)
    p = np.dstack((X.ravel(), Y.ravel()))[0]

    if center:
        if form == 'gaussian':
            func = gauss_center_p
            parnames = ['I0', 'sigma', 'x0', 'y0', 'offset']
        else:
            func = 0 # TODO : moffat_center
    else:
        if form == 'gaussian':
            func = gauss_p
            parnames = ['I0', 'sigma', 'offset']
        else:
            func = 0 # TODO moffat_p

    popt, pcov = optimize.curve_fit(func, p, chunk.ravel())

    print('2D fit, results : ')
    print('\t'.join(parnames))
    print('\t'.join(str(x) for x in popt))

    # If we have enabled centering then we retrieve the new values :
    if center:
        # TODO : retrieve center in image coordinates
        centerx, centery = popt[3:]
        popt = popt[0:3]

        

    # No background subtraction yet
    max_dist = math.sqrt(rplot**2+rplot**2)
    nfitpts = 100
    fx = np.linspace(0.0, max_dist, nfitpts)
    if form == 'gaussian':
        fy = gauss_r(fx, *popt)
    else:
        fy = 0 # TODO moffat_r

    dist = np.sqrt(p[:,0]**2+p[:,1]**2)

    if fig is None:
        fig = plt.figure(self._figure_name)

    fig.clf()
    fig.add_subplot(111)
    ax = fig.gca()
    ax.set_xlabel(params['xlabel'][0])
    ax.set_ylabel(params['ylabel'][0])
    if params["logx"][0]:
        ax.set_xscale("log")
    if params["logy"][0]:
        ax.set_yscale("log")

    ax.plot(dist, chunk.ravel(), params['marker'][0], label='data')
    ax.plot(fx, fy, c='r', label= form + " fit")
    plt.legend()
    plt.draw()
    plt.show(block=False)    
    time.sleep(self.sleep_time)

def astro_lab_register(viewer) :
    viewer.exam.register({'f': (fit_2d_new, 'Compute the 2d fit to the sample of data using the specified form')})
    viewer.exam.fit2d_pars = fit2d_pars
