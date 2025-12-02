import numpy as np
import click
import astropy.io.fits as fitsio
import logging
import pylab
import json

# create logger
log = logging.getLogger('simple_example')
log.setLevel(logging.DEBUG)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
log.addHandler(ch)

pi = np.pi
def mod2pi(theta):
    return theta%(2*pi)

def quadrant(theta):
    """
    Find quandrant of angle in 2D cartesian plane. w.r.t +ve axis in acw direction
    """
    theta = mod2pi(theta)
    if np.logical_and(theta > 0, theta < pi/2):
        return 1
    elif np.logical_and(theta > pi/2, theta < pi):
        return 2
    elif np.logical_and(theta > pi, theta < 0.75*pi):
        return 3
    else:
        return 4

def ends(theta, npix):
    """
    Find ends of a line in a square grid going through the centre and
    rotated by an angle theta w.r.t the horizontal (acw)

    Parameters:
    -----------

    theta (float): Angle in degrees
    npix (int): number of pixels in square grid (image)

    Returns:
    -----------

    ndarray: numpy.array([[x1,y1], [x2,y2]], dtype=int)
    """
    theta = mod2pi(np.deg2rad(theta))

    x0 = npix//2
    y0 = npix//2

    slope = -np.tan(theta) # slope of line with no y-intercept
    def fx(x):
        return slope*(x - x0) + y0
    def fy(y):
        return (y-y0)/slope + x0

    # avoid boundary
    def avoid(x):
        if x==npix:
            return x-1
        else:
            return x
    # Top/Bottom boundary
    if abs(slope) > 1:
        y1 = 0
        x1 = fy(y1)
        x1 = avoid(x1)

        y2 = npix - 1
        x2 = fy(y2)
        x2 = avoid(x2)
    # Left/Right boundary
    elif abs(slope) < 1:
        x1 = 0
        y1 = fx(x1)
        y1 = avoid(y1)

        x2 = npix - 1
        y2 = fx(x2)
        y2 = avoid(y2)

    if slope>0:
        return np.array([[x1,y1], [x2,y2]], dtype=int)
    else:
        return np.array([[x2,y2], [x1,y1]], dtype=int)
        
# nicked from https://stackoverflow.com/questions/47704008/fastest-way-to-get-all-the-points-between-two-x-y-coordinates-in-python
def crosseccros(ends):
    d0, d1 = np.diff(ends, axis=0)[0]
    if np.abs(d0) > np.abs(d1): 
        return np.c_[np.arange(ends[0, 0], ends[1,0] + np.sign(d0), np.sign(d0), dtype=np.int32),
                     np.arange(ends[0, 1] * np.abs(d0) + np.abs(d0)//2,
                               ends[0, 1] * np.abs(d0) + np.abs(d0)//2 + (np.abs(d0)+1) * d1, d1, dtype=np.int32) // np.abs(d0)]
    else:
        return np.c_[np.arange(ends[0, 0] * np.abs(d1) + np.abs(d1)//2,
                               ends[0, 0] * np.abs(d1) + np.abs(d1)//2 + (np.abs(d1)+1) * d0, d0, dtype=np.int32) // np.abs(d1),
                     np.arange(ends[0, 1], ends[1,1] + np.sign(d1), np.sign(d1), dtype=np.int32)]

def turning_points(profile, n):
    """
    Find first n turning points in a 1D array (profile)
    """

    size = len(profile)
    tp = []
    N = 0
    for i in range(1, size):
        if i+1 > size-1:
            log.debug("Current index larger than size of array [turning_points()]")
            break
        elif N == n:
            break
        elif (profile[i] < profile[i-1] and profile[i] < profile[i+1]) or\
            (profile[i] > profile[i-1] and profile[i] > profile[i+1]):

            tp.append(profile[i])
            N += 1
    log.info(f"Profile has {len(tp)} turning points")

    return tp


class PSFObject(object):
    def __init__(self, name, channel=0):
        self.name = name
        self.channel = channel
        self.angles = []
        self.profiles = []
        self.turning_points = []
        self.sidelobe_area = []

        with fitsio.open(self.name) as hdul:
            self.hdu = hdul[0]
            self.hdr = self.hdu.header

            naxis = self.hdr['NAXIS']
            # find frequency axis
            for i in range(1, naxis + 1):
                if self.hdr['CTYPE%d' % i].lower().startswith("freq"):
                    freq = naxis - i  # fits to numpy indexing

            ra, dec = naxis - 2, naxis - 1
            slc = [0]*naxis
            slc[ra] = slice(None)
            slc[dec] = slice(None)
            slc[freq] = self.channel


            print(slc, self.hdu.data.shape)
            self.data = self.hdu.data[tuple(slc)]

        self.npix = self.data.shape[0]
        self.scale = abs(self.hdr['CDELT1'])

    def set_profiles(self, angles=None):
        self.angles = angles or self.angles
        ends_ = None
        for theta in self.angles:
            # if too close (or equal) to zero deg
            if theta < self.scale/2:
                idxs = np.array([range(self.npix), [self.npix//2]*self.npix], dtype=int).T
            # if too close (or equal) to 90 deg
            elif theta < 91 and theta > 89:
                idxs = np.array([[self.npix//2]*self.npix, range(self.npix)], dtype=int).T
            else:
                ends_ = ends(theta, self.npix)
                idxs = crosseccros(ends_)
            profile = self.data[(idxs[:,0], idxs[:,1])]
            self.profiles.append(profile)

    def any_first_tp_above_tol(self, turning_points=None, tol=0.05):
        self.tuning_points = turning_points or self.turning_points
        for tps in self.turning_points:
            if abs(tps[0]) > tol:
                return True
        return False

    
@click.command()
@click.argument("psfs", nargs=-1, type=click.types.Path(writable=False))


@click.option("--channel-index", "-ci", "channel", type=int, default=0,
              help="Channel index if providing a multi-channel PSF. Default is 0")

@click.option("--ignore-psf-fit", "-ipf", "nopsf", is_flag=True, 
              help="Ignore fitted PSF positiona angle.")

@click.option("--angles", "-a", multiple=True, type=float,
              help="Angles (w.r.t horizontal in acw direction) which to take cross-sections." \
                   "If --ignore-psf-fit is not set, then the reference will be set to the fitted PSF position in the header")

@click.option("--max-tp", "--m","maxtp", default=6, type=int,
              help="Maximimm number of turning points to consider when comparing PSFs")
@click.option("--tolerance", "-tol","tol", default=0.05, type=float,
              help="Tolerance for how far from zero the first turning point of PSF can be along any angle")

@click.option("--save", "-s", 
              help="Filename for comparison results (JSON format). If not given a system default will be used.")
    
@click.option("--savefig", "-sf","savefig", is_flag=True,
              help="Save plot of comparison of qualifying psfs (cross-section along major axis will be used)")

def cli(psfs, nopsf, angles, maxtp, tol, save, savefig, channel):
    """
    Driver function
    """
    method = "turning_points"
    tps_psf = []
    angles = list(angles)
    save = save or "psftune_output.json"
    for name in psfs:
        psf = PSFObject(name, channel)
        log.info(f"Loading PSF: {psf.name}")

        if not angles:
            psf.angles = [-45, 0.0, 45, 90]

        if not nopsf:
            # use BPA as reference
            psf.angles = [a + psf.hdr["BPA"] for a in psf.angles]

        psf.set_profiles()
        for profile in psf.profiles:
            N = len(profile)
            if method == "turning_points":
                NN = N//2+2
                tp = turning_points(profile[NN:], maxtp)
                psf.turning_points.append(tp)

        # compute turning points for different cross-sections of a PSF
        if psf.any_first_tp_above_tol(tol=tol):
            log.info("PSF has at least one cross-section that has a 1st sidelobe that is too far from zero. This PSF does not make the cut")
        else:
            psf.metric = np.absolute(psf.turning_points).sum(axis=(0,1))
            tps_psf.append(psf)
    if len(tps_psf) == 0:
        raise RuntimeError("None of the input PSFs made cut.")
        
    neo = np.argmin([m.metric for m in tps_psf])
    names = [x.name for x in tps_psf]
    x = np.linspace(-psf.npix//2, psf.npix//2, psf.npix)*psf.scale
    L = psf.npix * psf.scale
    savejson = {
        "info": "PSF comparison output from psftune. The chosen one has psf.neo=True.",
        }

    
    for i,psf in enumerate(tps_psf):
        name = names[i]
        if i==neo:
            savejson[name] = dict(neo=True)
            log.info(f"The best PSF is: {name}")
            pt = "k-x"
        else:
            savejson[name] = dict(neo=False)
            pt = "-"
        if savefig:
            pylab.plot(x, psf.profiles[1], pt, label=name)
    
    with open(save, "w") as stdw:
        json.dump(savejson, stdw)

    if savefig:
        pylab.grid()
        pylab.ylim(-tol*3,tol*3)
        pylab.xlim(-L/5, L/5)
        pylab.title(f"Neo={names[neo]}")
        pylab.legend()
        pylab.savefig(f"{save}.png")
