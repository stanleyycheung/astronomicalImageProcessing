from astropy.io import fits
import numpy as np
from pprint import pprint

hdulist = fits.open('A1_mosaic.fits')
header = hdulist[0].header
pprint(header['MAGZPT'])
pprint(header['MAGZRR'])
pixelData = hdulist[0].data
exit()
np.save('pixelData', pixelData)
