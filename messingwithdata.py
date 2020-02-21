import json
import numpy as np
import matplotlib.pyplot as plt

with open(r"galaxies.json", "r") as read_file:
    galaxies = json.load(read_file)

ZP = 25.3
m = []
for galaxy in galaxies:
    total_count = galaxy['total_count']
    m.append(ZP - 2.5*np.log10(total_count-3421))

print(m)
