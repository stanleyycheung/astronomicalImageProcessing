import json
import numpy as np
import matplotlib.pyplot as plt

with open(r"galaxies.json", "r") as read_file:
    galaxies = json.load(read_file)
m_array = []

# m_min = 100
# m_min_pos = [0, 0]
for galaxy in galaxies:
    m = galaxy['m']
#     if m < m_min:
#         m_min = m
#         m_min_pos = galaxy['pos']
    m_array.append(m)

# img = np.load('realmaskedData.npy')
# offset = 23.5
# plt.imshow(img[int(m_min_pos[0]-offset):int(m_min_pos[0]+offset),
#                int(m_min_pos[1]-offset):int(m_min_pos[1]+offset)])
# plt.show()
# exit()

lower = 20
upper = 1700

sorted_m = np.sort(m_array)
count = np.arange(1, len(sorted_m) + 1)
count = np.log10(count)

print(sorted_m)

exit()

plt.plot(sorted_m, count, 'x')
plt.plot(sorted_m[lower], count[lower], 'o')
plt.plot(sorted_m[-upper], count[-upper], 'o')

fit_m = sorted_m[lower:-upper]
fit_count = count[lower:-upper]

z, cov = np.polyfit(fit_m, fit_count, 1, cov=True)
p = np.poly1d(z)


print(z)
plt.plot(sorted_m, p(sorted_m))


# plt.yscale('log')

plt.xlabel('m')
plt.ylabel(r"$\log N(m'<m)$")
plt.tight_layout()
plt.show()
