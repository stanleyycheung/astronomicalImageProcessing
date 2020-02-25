import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14})
x = np.load('background_whole.npy')
plt.hist(x, density=True, bins=3550-3340+1)
plt.xlabel(r'Count per pixel ($n_p$)')
plt.ylabel(r'P($n_p$)')
plt.tight_layout()
plt.savefig('fig/before_correction.pdf')
plt.show()
