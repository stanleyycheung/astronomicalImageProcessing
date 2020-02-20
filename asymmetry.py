import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

galaxies_points = np.load('galaxies_points.npy',allow_pickle=True)
img = np.load('testData_noisy.npy')
#print(galaxies_points)

x_resid = []
y_resid = []

for galaxy in galaxies_points:
    brightestValue = 0
    brightestPoint = [0, 0]
    x_min, x_max = 1e9, 0
    y_min, y_max = 1e9, 0
    for point in galaxy:
        brightness = img[point[0], point[1]]
        if brightness > brightestValue:
            brightestValue = brightness
            brightestPoint = point
        if point[0] < x_min:
            x_min = point[0]
        elif point[0] > x_max:
            x_max = point[0]
        if point[1] < y_min:
            y_min = point[1]
        elif point[1] > y_max:
            y_max = point[1]
    x_mid, y_mid = (x_max + x_min)/2, (y_max + y_min)/2
    # radius = max((x_max-x_mid), (x_mid-x_min), (y_max-y_mid), (y_mid-y_min))

    x_diff = abs((x_max-brightestPoint[0])-(brightestPoint[0]-x_min)) #difference in the distance to brightest point from xmin and xmax\
    y_diff = abs((y_max-brightestPoint[1])-(brightestPoint[1]-y_min))

    # #for later in the calc galaxies
    # if x_diff > 3 or y_diff > 3:
    #     return False

    x_resid.append(x_diff)
    y_resid.append(y_diff)

    # if x_diff > 5:
    #     print(x_min, x_max)
    #     print(y_min, y_max)
    #     print(x_mid, y_mid)
    #     fig,ax = plt.subplots()
    #     plt.imshow(img[x_min:x_max, y_min:y_max], norm=LogNorm())

    if y_diff > 5 or x_diff > 5:
        print(x_min, x_max)
        print(y_min, y_max)
        print(x_mid, y_mid)
        fig,ax = plt.subplots()
        plt.imshow(img[x_min:x_max, y_min:y_max], norm=LogNorm())

# print('x diff are', x_resid)
# print('y diff are', y_resid)

fig,ax = plt.subplots()
plt.hist(y_resid, bins ='auto')
plt.title("difference in y distances")
plt.xlabel("difference in pixel length")
plt.ylabel("number of counts")
plt.show()
