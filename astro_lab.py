import sys, os

# Pre-pending everything. Shhhh nobody will know.
sys.path.insert(0, '/user/phstf/md0046/tools/xpa/build/lib')
sys.path.insert(0, '/user/phstf/md0046/tools/astropy/build/lib.linux-x86_64-2.7')
sys.path.insert(0, '/user/phstf/md0046/tools/imexam/build/lib/python2.7/site-packages/imexam-0.5.2-py2.7-linux-x86_64.egg')
sys.path.insert(0, '/user/phstf/md0046/tools/photutils/build/lib.linux-x86_64-2.7/')

# Importing basic utils
import pyds9
import imexam
import pyfits as fits
import numpy as np
from scipy import optimize
import math
from ctypes import c_longlong as ll

fit2d_pars = {
    'title': ['', 'Title of the plot'],
    'xlabel': ['Radius', 'X-axis label'], 
    'ylabel': ['Pixel value', 'Y-axis label'],
    'fitplot': [True, 'Overplot profile fit?'],
    'fittype': ['moffat', 'Profile type to fit'],
    'center': [False, 'Center object in aperture?'],
    'background': [True, 'Fit and subtract background?'],
    'radius': [17., 'Object radius'],
    'buffer': [5., 'Background buffer width'],
    'width': [17., 'Background width'],
    'iterations': [5, 'Number of radius adjustment iterations'],
    'xorder': [0, 'Background x order'],
    'yorder': [0, 'Background y order'],
    'magzero': [25., 'Magnitude zero point'],
    'beta': ['INDEF', 'Moffat beta parameter'],
    'rplot': [30., 'Plotting radius'],
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

# Functions used for the lab that are not present (or weren't found in astropy)
def combine(file_list, ctype='median'):
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
    fout = open(file_out, 'w')
    for f in os.listdir(folder):
        if f[-5:] == '.fits' and f[0] != '.':
            header = fits.getheader(folder+'/'+f)
            fout.write('{0} {1}\n'.format(header['DATE-OBS'], header['EXPTIME']))

    fout.close()


def ie_center(im, radius, xcntr, ycntr):
    """
    Find the center of a star image given approximate coords.  Uses
    Mountain Photometry Code Algorithm as outlined in Stellar Magnitudes
    from Digital Images.
    """
    for k in range(3):
        # Extract region around center
        xlast = xcntr
        ylast = ycntr
        x1 = int(xcntr - radius + 0.5)
        x2 = int(xcntr + radius + 0.5)
        y1 = int(ycntr - radius + 0.5)
        y2 = int(ycntr + radius + 0.5)
        data = im[y1:y2+1, x1:x2+1]

        ny, nx = data.shape
        npts = nx * ny

        # Find center of gravity for marginal distributions above mean.
        S = data.sum()
        mean = S / nx
        sum1 = 0
        sum2 = 0

        for i in range(nx):
            sum3 = data[:,i].sum() - mean
            if sum3 > 0:
                sum1 += i * sum3
                sum2 += sum3

        if abs(sum2) > 1e-5:
            xcntr = sum1 / sum2

        mean = S / ny
        sum1 = 0.
        sum2 = 0.

        for i in range(ny):
            sum3 = data[i,:].sum() - mean
            if sum3 > 0:
                sum1 += i * sum3
                sum2 += sum3

        if abs(sum2) > 1e-5:
            ycntr = sum1 / sum2

        if abs(xcntr - xlast) < 1e-5 and abs(ycntr - ylast) < 1e-5:
            break

    return xcntr, ycntr


def gs_eval(x, y, sigma):
    return np.exp(-0.5*(x**2 + y**2) / sigma**2)

def gs_eval_r(r, i0, sigma):
    return i0 * np.exp(-0.5*r**2 / sigma**2)

def moffat_eval_r(r, i0, sigma, b):
    y = 1.0 + (r / sigma) ** 2
    return i0 * y**b

def ie_gauss(x, i0, sigma):
    r = np.sqrt(x[:,0]**2 + x[:,1]**2)
    r2 = r**2 / (2.0*sigma)
    return i0 * np.exp(-r2)

beta = 0.0
def ie_moffat_fixed_b(x, i0, sigma):
    r = np.sqrt(x[:,0]**2 + x[:,1]**2)
    y = 1.0 + (r / sigma) ** 2
    return i0 * y**beta

def ie_moffat_b(x, i0, sigma, b):
    r = np.sqrt(x[:,0]**2 + x[:,1]**2)
    y = 1.0 + (r / sigma) ** 2
    return i0 * y**b

def fit_2d(self, x, y, im, form=None, subsample=4, fig=None):
        """Compute the 2d fit to the sample of data using the specified form


        Parameters
        ----------
        form: string
            This is the functional form specified line fit parameters
            Currently 'gaussian' or 'moffat'

        """
        cnt = 0
        print('Dbg : ' + str(cnt))
        cnt+=1

        amp = 0
        sigma = 0
        params = self.fit2d_pars
        print(params)

        if not form:
            form = params["fittype"][0]

        if form not in ('gaussian', 'moffat'):
            warnings.warn('Error : The fitting function must be either \'gaussian\' or \'moffat\'')
            return

        print('Dbg : ' + str(cnt))
        cnt+=1  

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

	# Center
	if center:
            xcntr, ycntr = ie_center(im, radius, x, y) # Does not work yet
        else:
            xcntr, ycntr = x, y

        print('Dbg : ' + str(cnt))
        cnt+=1

        """ No idea the purpose of this part ...
        # Do the enclosed flux and direct FWHM measurments using the
	# PSFMEASURE routines.

	call stf_measure (im, xcntr, ycntr, beta, 0.5, radius, nit, buffer,
	    width, INDEF, NULL, NULL, dbkg, r, dfwhm, gfwhm, efwhm)
	if (fittype == FITGAUSS)
	    efwhm = gfwhm             
        """
        dfwhm = 1.0 # TODO : compute the correct fhwm

        # Get data including a buffer and background annulus.
        if not background:
            buffer = 0.0
            width = 0.0

        print('Dbg : ' + str(cnt))
        cnt+=1

        import numpy as np 

        r = max(rplot, radius + buffer + width)
        print(rplot, radius + buffer + width, r)
	x1 = int(xcntr - r)
	x2 = int(xcntr + r)
	y1 = int(ycntr - r)
	y2 = int(ycntr + r)
        print(x1, x2, y1, y2)
        data = im[y1:y2+1, x1:x2+1]
        lindat = np.ravel(data)
        ny, nx = data.shape
        npts = nx*ny

        print(npts)
        xs = np.zeros((npts,), dtype=np.float32)
        ys = np.zeros((npts,), dtype=np.float32)
        zs = np.zeros((npts,), dtype=np.float32)

        print('Dbg : ' + str(cnt))
        cnt+=1

        # Extract the background data if background subtracting.
	ns = 0
        lid = 0
	if background and width > 0:
	    r1 = radius**2
	    r2 = (radius + buffer)**2
	    r3 = (radius + buffer + width)**2

            for j in range(y1, y2+1):
	        dy = (ycntr - j)**2
                for i in range(x1, x2+1):
		    r = (xcntr - i)**2 + dy
		    if r <= r1:
                        continue
		    elif r >= r2 and r <= r3:
		        xs[ns] = i
		        ys[ns] = j
		        zs[ns] = lindat[lid]
		        ns += 1

		    lid += 1

        
        print('Dbg : ' + str(cnt))
        cnt+=1

        # Accumulate the various sums for the moments and the gaussian fit.
	no = 0
	nps = 0
	zcntr = 0.
	sumo = 0. 
        sums = 0.
        sumxx = 0.
        sumyy = 0.
        sumxy = 0.
	lid = 0

        # Background subtraction :
        if ns > 0:
	    # If background points are defined fit a surface and subtract
	    # the fitted background from within the object aperture. 

	    if medsky:
                print(zs[:ns])
                tmp = sorted(list(zs[:ns]))
		bkg = tmp[len(tmp)//2] # np.median(zs[:ns], axis=1)
	    else:
                done = False
                ## Fitting a gaussian curve
                
                ## Init
                p = np.dstack(xs, ys).T
                popt, = optimize.curve_fit(gs_eval, p, zs)
                gs = popt[0]
		bkg = gseval(gs, real(x1), real(y1))
	    
            for j in range(y1, y2+1):
	        dy = j - ycntr

                for i in range(x1, x2+1):
		    dx = i - xcntr
		    r = math.sqrt(dx ** 2 + dy ** 2)
		    r3 = max(0., min(5., 2 * r / dfwhm - 1.))

		    if medsky:
			r2 = bkg
		    else:
			r2 = gseval(gs, i, j)
			bkg = min(bkg, r2)
		    r1 = lindat[lid] - r2

		    if r <= radius:
			sumo += r1
			sums += r2
			sumxx += dx * dx * r1
			sumyy += dy * dy * r1
			sumxy += dx * dy * r1
			zcntr = max(r1, zcntr)
			if r <= rplot:
                            xs[no] = r
                            ys[no] = r1
			    zs[no] = math.exp(-r3**2) / max (.1, r**2)
			    no += 1
			else:
			    nps += 1
                            xs[npts-nps] = r
			    ys[npts-nps] = r1
			    zs[npts-nps] = math.exp (-r3**2) / max (.1, r**2)
			
		    elif r <= rplot:
		        nps += 1
		        xs[npts-nps] = r
		        ys[npts-nps] = r1
		    
		    lid += 1

	else:
            bkg = 0.
            for j in range(y1, y2-1):
	        dy = j - ycntr
                for i in range(x1, x2-1):
		    dx = i - xcntr
		    r = np.sqrt (dx ** 2 + dy ** 2)
		    r3 = max (0., min (5., 2 * r / dfwhm - 1.))
		    r1 = lindat[lid]

		    if r <= radius:
			sumo += r1
			sumxx += dx * dx * r1
			sumyy += dy * dy * r1
			sumxy += dx * dy * r1
			zcntr = max (r1, zcntr)
			if r <= rplot:
			    xs[no] = r
			    ys[no] = r1
			    zs[no] = math.exp(-r3**2) / max (.1, r**2)
			    no += 1
			else:
			    nps += 1
			    xs[npts-nps] = r
			    ys[npts-nps] = r1
			    zs[npts-nps] = math.exp(-r3**2) / max (.1, r**2)
			
		    elif r <= rplot:
		        nps += 1
		        xs[npts-nps] = r
		        ys[npts-nps] = r1

                    lid += 1

        print('Dbg : ' + str(cnt))
        cnt+=1
        

        # What's that ?
        if nps > 0:
            xs[npts-nps:npts] = xs[no:no+nps]
            ys[npts-nps:npts] = xs[no:no+nps]
            zs[npts-nps:npts] = xs[no:no+nps]
        
	if rplot <= radius:
	    no += nps
	    nps = no - nps 
	else:
	    nps += no

        print(xs.shape, ys.shape)
        p = np.dstack((xs, ys))
        # Compute the photometry and profile fit parameters.
        if form == 'gaussian':
            par = [zcntr, dfwhm**2 / (8 * log(2.))]
            popt, = optimize.curve_fit(ie_gauss, p, zs, p0=par)
            zcntr = popt[0]
            fwhm = math.sqrt(8.0*log(2.0)*popt[1])
        elif form == 'moffat':
            global beta
            if loc_beta == 'INDEF':
                beta = -3.0
                par = [zcntr, dfwhm / 2. / math.sqrt (2.**(-1./beta) - 1.)]
                func = ie_moffat_fixed_b
            else:
                beta = -loc_beta
                par = [zcntr, dfwhm / 2. / math.sqrt (2.**(-1./beta) - 1.), beta]
                func = ie_moffat_b
            
            popt, = optimize.curve_fit(ie_gauss, p, zs, p0=par)
            if popt[1] < 0:
                zcntr = None
                fhwm = None
                beta = None
            else:
                zcntr = popt[0]
                beta  = -popt[2]
                fwhm = abs(popt[1])*2.0*sqrt(2.0**(-1.0/popt[2])-1.0)

        if sumo > 0:
            mag = magzero - 2.5 * math.log(sumo, 10)
	    r2 = sumxx + sumyy
	    if r2 > 0:
                if form == 'gaussian':
		    r = 2 * math.sqrt (log (2.) * r2 / sumo)
                elif form == 'moffat':
		    if beta > 2:
			r = 2 * math.sqrt ((beta-2.)*(2.**(1./beta)-1) * r2 / sumo)

	        r1 =(sumxx-sumyy)**2+(2*sumxy)**2
		if r1 > 0:
		    e = math.sqrt (r1) / r2
	        else:
		    e = 0.
	    
	    if e < 0.01:
		e = 0.
	    else:
		pa = 180.0/math.pi * (0.5 * atan2 (2*sumxy, sumxx-sumyy))

        '''
        call ie_mwctran (ie, xcntr, ycntr, wxcntr, wycntr)
	if (xcntr == wxcntr && ycntr == wycntr)
	    call strcpy ("%.2f %.2f", Memc[title], IE_SZTITLE)
	else {
	    call sprintf (Memc[title], IE_SZTITLE, "%s %s")
		if (IE_XFORMAT(ie) == '%')
		    call pargstr (IE_XFORMAT(ie))
		else
		    call pargstr ("%g")
		if (IE_YFORMAT(ie) == '%')
		    call pargstr (IE_YFORMAT(ie))
		else
		    call pargstr ("%g")
	}
	call sprintf (Memc[coords], IE_SZTITLE, Memc[title])
	    call pargr (wxcntr)
	    call pargr (wycntr)
        '''
        if fig is None:
            fig = plt.figure(self._figure_name)

        dist = np.zeros((npts,))
        val  = np.zeros((npts,))
        
        lid = 0
        for i in xrange(x1, x2+1):
            for j in xrange(y1, y2+1):
                dist[j, i] = math.sqrt((xcntr - i)**2 + (ycntr - j)**2)
                val = data[j, i]
        
        ncurve = 100
        gx = np.linspace(dist.min(), dist.max(), num=ncurve)
        if form == 'gaussian':
            gy = gs_eval_r(dist, zcntr, fhwm)
        else:
            gy = gs_eval_r(dist, zcntr, fhwm, beta)

        fig.clf()
        fig.add_subplot(111)
        ax = fig.gca()
        ax.set_xlabel(params['xlabel'][0])
        ax.set_ylabel(params['ylabel'][0])
        if params["logx"][0]:
            ax.set_xscale("log")
        if params["logy"][0]:
            ax.set_yscale("log")

        ax.plot(dist, val, chunk, params['marker'], label='data')
        ax.plot(gx, gy, c='r', label= form + " fit")
        plt.legend()
        plt.draw()
        plt.show(block=False)

        # ------ RESUSE SOME OF THIS
        '''
        delta = self.line_fit_pars["rplot"][0]

        if self.line_fit_pars["center"][0]:
            popt = self.gauss_center(x, y, data, delta=delta)
            if popt.count(0) > 1:  # an error occurred in centering
                centerx = x
                centery = y
                warnings.warn(
                    "Problem fitting center, using original coordinates")
            else:
                amp, x, y, sigma, offset = popt

        

        line = data[y, :]
        chunk = line[x - delta:x + delta]

        # use x location as the first estimate for the mean, use 20 pixel
        # distance to guess center
        argmap = {'a': amp, 'mu': len(chunk) / 2, 'sigma': sigma, 'b': 0}
        args = inspect.getargspec(form)[0][1:]  # get rid of "self"

        xline = np.arange(len(chunk))
        popt, pcov = curve_fit(
            form, xline, chunk, [
                argmap.get(
                    arg, 1) for arg in args])
        # do these so that it fits the real pixel range of the data
        fitx = np.arange(len(xline), step=1. / subsample)
        fity = form(fitx, *popt)

        # calculate the std about the mean
        # make a plot
        if fig is None:
            fig = plt.figure(self._figure_name)
        fig.clf()
        fig.add_subplot(111)
        ax = fig.gca()
        ax.set_xlabel(self.line_fit_pars["xlabel"][0])
        ax.set_ylabel(self.line_fit_pars["ylabel"][0])
        if self.line_fit_pars["logx"][0]:
            ax.set_xscale("log")
        if self.line_fit_pars["logy"][0]:
            ax.set_yscale("log")

        if self.line_fit_pars["pointmode"][0]:
            ax.plot(xline + x - delta, chunk, 'o', label="data")
        else:
            ax.plot(xline + x - delta, chunk, label="data", linestyle='-')

        if self.line_fit_pars["func"][0] == "gaussian":
            sigma = np.abs(popt[2])
            fwhm = math_helper.gfwhm(sigma)
            fitmean = popt[1] + x - delta
        elif self.line_fit_pars["func"][0] == "moffat":
            alpha = popt[0]
            beta = popt[1]
            fwhm = math_helper.mfwhm(alpha, beta)
            fitmean = popt[2] + x - delta
        else:
            warnings.warn("Unsupported functional form used in line_fit")
            raise ValueError

        ax.set_title("{0:s}  mean = {1:s}, fwhm= {2:s}".format(
            self.line_fit_pars["title"][0], str(fitmean), str(fwhm)))
        ax.plot(
            fitx +
            x -
            delta,
            fity,
            c='r',
            label=str(
                form.__name__) +
            " fit")
        plt.legend()
        plt.draw()
        plt.show(block=False)
        time.sleep(self.sleep_time)
        pstr = "({0:d},{1:d}) mean={2:0.3f}, fwhm={3:0.3f}".format(
            int(x + 1), int(y + 1), fitmean, fwhm)
        print(pstr)
        logging.info(pstr)
        '''

def astro_lab_register(viewer) :
    viewer.exam.register({'f': (fit_2d, 'Compute the 2d fit to the sample of data using the specified form')})
    viewer.exam.fit2d_pars = fit2d_pars

# Remove this later on
ds9 = pyds9.DS9('lab')
view = imexam.connect('lab')
astro_lab_register(view)
data = fits.getdata('east-14/2014_02_07/d0004.fits')
view.view(data)
view.imexam()
