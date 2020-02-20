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


def midPointAlgorithm(data, center, radius, value):
    f = 1 - radius
    ddf_x = 1
    ddf_y = -2 * radius
    x = 0
    y = radius
    data[center[0], center[1] + radius] = value
    data[center[0], center[1] - radius] = value
    data[center[0] + radius, center[1]] = value
    data[center[0] - radius, center[1]] = value

    while x < y:
        if f >= 0:
            y -= 1
            ddf_y += 2
            f += ddf_y
        x += 1
        ddf_x += 2
        f += ddf_x
        data[center[0] + x, center[1] + y] = value
        data[center[0] - x, center[1] + y] = value
        data[center[0] + x, center[1] - y] = value
        data[center[0] - x, center[1] - y] = value
        data[center[0] + y, center[1] + x] = value
        data[center[0] - y, center[1] + x] = value
        data[center[0] + y, center[1] - x] = value
        data[center[0] - y, center[1] - x] = value
    return data


if __name__ == '__main__':
    img = np.load('testData_noisy.npy')
    galaxies_points = np.load('galaxies_points.npy', allow_pickle=True)
    # print(img.shape)
    # print(len(galaxies_points))
    # radial_plot(galaxies_points, img)

    data = np.zeros((40, 40))
    data = midPointAlgorithm(data, [20, 20], 10, 1)
    data = midPointAlgorithm(data, [20, 20], 15, 2)
    plt.imshow(data, origin='lower')
    plt.show()
