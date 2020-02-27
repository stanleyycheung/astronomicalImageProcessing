import operator
import json
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14})
with open(r"galaxies.json", "r") as read_file:
    galaxies = json.load(read_file)
m_array = []
m_err_array = []

for galaxy in galaxies:
    m = galaxy['m']
    m_err = galaxy['m_err']
    m_array.append(m)
    m_err_array.append(m_err)

L = sorted(zip(m_array, m_err_array), key=operator.itemgetter(0))
sorted_m, m_err_array_sorted = zip(*L)


m_step = 0.2
m_plot = np.arange(int(sorted_m[0])+1, sorted_m[-1]+m_step, step=m_step)
count = []
for i in range(len(m_plot)):
    c = 0
    for j in sorted_m:
        if m_plot[i] > j:
            c += 1
    count.append(c)

count = np.array(count)
count_err = np.sqrt(count)
count_err_new = [[], []]
count_err = 1/count * count_err
count = np.log10(count)
for i in range(len(count_err)):
    if count[i] - count_err[i] > 0:
        count_err_new[0].append(count_err[i])
        count_err_new[1].append(count_err[i])
    else:
        count_err_new[0].append(count[i])
        count_err_new[1].append(count_err[i])

# plt.plot(m_plot, np.log10(count), 'x')
plt.errorbar(m_plot, count, fmt='x', yerr=count_err_new, capsize=3, zorder=1)
lower = 8
upper = 26
fit_m = m_plot[lower:-upper]
fit_count = count[lower:-upper]
z, cov = np.polyfit(fit_m, fit_count, 1, cov=True, w=1/count_err[lower:-upper])
p = np.poly1d(z)

print(z, cov[0][0])
plt.errorbar(m_plot[-upper], count[-upper], fmt='x',
             yerr=count_err[-upper], capsize=3, color='red', zorder=2)
plt.errorbar(m_plot[lower], count[lower], fmt='x',
             yerr=count_err[lower], capsize=3, color='red', zorder=2)
plt.plot(m_plot[:45], p(m_plot[:45]), zorder=3)
plt.xlabel('m')
plt.ylabel(r"$\log N(m)$")
plt.grid()
plt.savefig('fig/distribution.pdf')
plt.tight_layout()
plt.show()
