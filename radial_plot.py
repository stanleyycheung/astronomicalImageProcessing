import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter


def radial_plot(galaxies_points, img):
    for i in range(100, 200):
        galaxy = galaxies_points[i]
        brightestValue = 0
        brightestPoint = [0, 0]
        x_min, x_max = 1e9, 0
        y_min, y_max = 1e9, 0
        for point in galaxy:
            brightness = img[point[0], point[1]]
            if point[0] < x_min:
                x_min = point[0]
            elif point[0] > x_max:
                x_max = point[0]
            if point[1] < y_min:
                y_min = point[1]
            elif point[1] > y_max:
                y_max = point[1]
            if brightness > brightestValue:
                brightestValue = brightness
                brightestPoint = point
        # print(brightestPoint, brightestValue)
        # print(midPointAlgorithm(img, brightestPoint, 1))
        # print(midPointAlgorithm(img, brightestPoint, 2))
        radial_brightness = [brightestValue]
        for i in range(1, 11):
            radial_brightness.append(midPointAlgorithm(img, brightestPoint, i))
        r_plot = np.arange(0, 11)
        smoothed_radial_brightness = savgol_filter(radial_brightness, 7, 3)

        if sorted(np.diff(smoothed_radial_brightness))[-1] > 5:
            fig, ax = plt.subplots()
            plt.plot(r_plot, radial_brightness, '-')
            plt.plot(r_plot, smoothed_radial_brightness, 'o-')
            fig, ax = plt.subplots()
            plt.imshow(img[x_min:x_max+1, y_min:y_max+1], origin='upper', norm=LogNorm())


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
    return np.average(pixelCounts)


if __name__ == '__main__':
    img = np.load('realmaskedData.npy')
    galaxies_points = np.load('filtered_galaxies_points.npy', allow_pickle=True)
    # print(img.shape)
    # print(len(galaxies_points))
    radial_plot(galaxies_points, img)
    plt.show()
