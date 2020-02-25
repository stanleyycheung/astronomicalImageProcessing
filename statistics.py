import json
import numpy as np
import matplotlib.pyplot as plt
import operator

with open(r"galaxies.json", "r") as read_file:
    galaxies = json.load(read_file)
m_array = []
m_err_array = []

# m_min = 100
# m_min_pos = [0, 0]
for galaxy in galaxies:
    m = galaxy['m']
    m_err = galaxy['m_err']
    # real_count = galaxy['real_count']
    # if m < m_min:
    #     m_min = m
    #     m_min_pos = galaxy['pos']
    m_array.append(m)
    m_err_array.append(m_err)
    # if galaxy['avg_background'] > 3450:
    #     print(galaxy)

# img = np.load('realmaskedData.npy')
# offset = 200
# plt.imshow(img[int(m_min_pos[0]-offset):int(m_min_pos[0]+offset),
#                int(m_min_pos[1]-offset):int(m_min_pos[1]+offset)])
# plt.show()
# exit()

lower = 20
upper = 1150

L = sorted(zip(m_array, m_err_array), key=operator.itemgetter(0))
sorted_m, m_err_array_sorted = zip(*L)

sorted_m = sorted_m[3:]
m_err_array_sorted = m_err_array_sorted[3:]

# specific_m = sorted_m[8]
# for galaxy in galaxies:
#     if galaxy['m'] == specific_m:
#         x, y = galaxy['pos']
# print(specific_m)
# img = np.load('realmaskedData.npy')
# offset = 20
# plt.imshow(img[int(x-offset):int(x+offset),
#                int(y-offset):int(y+offset)])
# plt.show()
# exit()
count = np.arange(1, len(sorted_m)+1)
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


plt.errorbar(sorted_m, count, fmt='x', yerr=count_err_new, capsize=3)
# plt.plot(sorted_m, count, 'x')

fit_m = sorted_m[lower:-upper]
fit_count = count[lower:-upper]

z, cov = np.polyfit(fit_m, fit_count, 1, cov=True, w=1/count_err[lower:-upper])
p = np.poly1d(z)

print(z, cov[0][0])

plt.plot(sorted_m, p(sorted_m))
# plt.plot(sorted_m[lower], count[lower], 'o')
# plt.plot(sorted_m[-upper], count[-upper], 'o')
# plt.yscale('log')

plt.xlabel('m')
plt.ylabel(r"$\log N(m'<m)$")
plt.savefig('fig/distribution.pdf')
plt.tight_layout()
plt.show()
