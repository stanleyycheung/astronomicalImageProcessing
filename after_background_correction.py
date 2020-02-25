import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
plt.rcParams.update({'font.size': 14})


def gaussian(x, mean, amplitude, standard_deviation):
    return amplitude * 1/(standard_deviation*2*np.sqrt(np.pi)) * np.exp(-0.5*((x - mean) / standard_deviation) ** 2)


# data = fits.open('masked1.fits')[0].data
# xlower = 3340
# xupper = 3550
# # xlower, xupper = 3400, 3440
# background = [i for i in data.flatten() if xlower <= i <= xupper]
# np.save('background_whole_correction', background)


x = np.load('background_whole_correction.npy')
bin_heights, bin_borders, _ = plt.hist(x, bins=3550-3340+1, density=True)
bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
popt, _ = curve_fit(gaussian, bin_centers, bin_heights, p0=[3420, 0.03, 50])
print(popt)
print(_)
# plt.plot(bin_centers, bin_heights, 'o')
x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)
plt.plot(x_interval_for_fit, gaussian(x_interval_for_fit, *popt), label='fit')
plt.xlabel(r'Count per pixel ($n_p$)')
plt.ylabel(r'P($n_p$)')
plt.tight_layout()
plt.savefig('fig/after_correction.pdf')
plt.show()
