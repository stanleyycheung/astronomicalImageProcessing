import json
import numpy as np
import matplotlib.pyplot as plt

with open(r"galaxies.json", "r") as read_file:
    galaxies = json.load(read_file)
real_count = []

# m_min = 100
# m_min_pos = [0, 0]
for galaxy in galaxies:
    real_count.append(galaxy['real_count'])

plt.hist(real_count, bins=40)
plt.show()
