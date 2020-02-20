import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def radial_plot(galaxies_points, img):
    for galaxy in [galaxies_points[0]]:
        brightestValue = 0
        brightestPoint = [0, 0]
        for point in galaxy:
            brightness = img[point[0], point[1]]
            if brightness > brightestValue:
                brightestValue = brightness
                brightestPoint = point
        print(brightestPoint, brightestValue)
        print(midPointAlgorithm(img, brightestPoint, 1))


def midPointAlgorithm(data, center, radius):
    pixelCounts = []
    f = 1 - radius
    ddf_x = 1
    ddf_y = -2 * radius
    x = 0
    y = radius
    pixelCounts.append(data[center[0], center[1] + radius])
    pixelCounts.append(data[center[0], center[1] - radius])
    pixelCounts.append(data[center[0] + radius, center[1]])
    pixelCounts.append(data[center[0] - radius, center[1]])

    while x < y:
        if f >= 0:
            y -= 1
            ddf_y += 2
            f += ddf_y
        x += 1
        ddf_x += 2
        f += ddf_x
        pixelCounts.append(data[center[0] + x, center[1] + y])
        pixelCounts.append(data[center[0] - x, center[1] + y])
        pixelCounts.append(data[center[0] + x, center[1] - y])
        pixelCounts.append(data[center[0] - x, center[1] - y])
        pixelCounts.append(data[center[0] + y, center[1] + x])
        pixelCounts.append(data[center[0] - y, center[1] + x])
        pixelCounts.append(data[center[0] + y, center[1] - x])
        pixelCounts.append(data[center[0] - y, center[1] - x])
    return pixelCounts


if __name__ == '__main__':
    img = np.load('testData_noisy.npy')
    galaxies_points = np.load('galaxies_points.npy', allow_pickle=True)
    # print(img.shape)
    # print(len(galaxies_points))
    # radial_plot(galaxies_points, img)
