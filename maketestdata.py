from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

hdulist = fits.open('A1_mosaic.fits')
pixelData = hdulist[0].data

testData = pixelData[1878:2168, 1598:1961]  # small slice
# testData = pixelData[1963:2446, 1591:2244]
testData = pixelData[2240:2631, 1574:1790]  # bigger slice
np.save('testData_bigger', testData)

fig, ax = plt.subplots()
plt.imshow(testData, norm=LogNorm())
plt.colorbar()
plt.show()
