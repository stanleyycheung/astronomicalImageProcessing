import json
import numpy as np
import matplotlib.pyplot as plt

with open(r"galaxies.json", "r") as read_file:
    galaxies = json.load(read_file)
m_array = []
m_step = 0.1
for galaxy in galaxies:
    m_array.append(galaxy['m'])

sorted_m = np.sort(m_array)
count = np.arange(1, len(sorted_m) + 1)
print(sorted_m)
print(count)
