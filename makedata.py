from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# hdulist = fits.open('A1_mosaic.fits')
# pixelData = hdulist[0].data
#
# testData = pixelData[1878:2168, 1598:1961]  # small slice
# # testData = pixelData[1963:2446, 1591:2244]
# testData = pixelData[2240:2631, 1574:1790]  # bigger slice

mData = fits.open('masked1.fits')[0].data
orgData = fits.open('A1_mosaic.fits')[0].data
print(mData.shape)
exit()
np.save('realmaskedData', mData)
np.save('realOrgData', orgData)

fig, ax = plt.subplots()
plt.imshow(orgData, norm=LogNorm(), origin='lower')
plt.colorbar()
plt.show()
