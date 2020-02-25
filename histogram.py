import numpy as np
import matplotlib.pyplot as plt

differences = np.load('differences.npy')

# lowrangeplot = []
# for i in range(len(differences)):
#     if differences[i] <200:
#         lowrangeplot.append(differences[i])

plt.hist([i for i in differences if i <= 250000], bins=200)
#plt.hist(differences, bins = 100)
plt.show()
