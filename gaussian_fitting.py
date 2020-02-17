import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def gauss2d(xy, amp, x0, y0, a, b, c):
    x, y = xy
    inner = a * (x - x0)**2
    inner += 2 * b * (x - x0)**2 * (y - y0)**2
    inner += c * (y - y0)**2
    return amp * np.exp(-inner)


def gaussian_fitting(galaxies_points, img):
    for counter in range(1):
        galaxy = galaxies_points[counter]
        x_min, x_max = 1e9, 0
        y_min, y_max = 1e9, 0
        for point in galaxy:
            if point[0] < x_min:
                x_min = point[0]
            elif point[0] > x_max:
                x_max = point[0]
            if point[1] < y_min:
                y_min = point[1]
            elif point[1] > y_max:
                y_max = point[1]
        buffer = 3
        x_min = max(int(x_min-buffer), 0)
        # x_max = min(int(x_max+buffer), img.shape[1] - 1)
        y_min = max(int(y_min-buffer), 0)
        # y_max = min(int(y_max+buffer), img.shape[0] - 1)
        # print(x_min, x_max, y_min, y_max)
        if x_min == x_max or y_min == y_max:
            print(counter)

        x = np.arange(x_min, x_max, 1)
        y = np.arange(y_min, y_max, 1)
        xy = np.meshgrid(y, x)

        img_slice = img[x_min:x_max, y_min:y_max]
        print(len(xy[0]))
        print(img_slice.shape)
        popt, pcov = curve_fit(gauss2d, xy, img_slice)
        fig, ax = plt.subplots()
        plt.imshow(img[x_min:x_max, y_min:y_max], origin='upper')


if __name__ == '__main__':
    img = np.load('testData_noisy.npy')
    galaxies_points = np.load('galaxies_points.npy', allow_pickle=True)
    # print(img.shape)
    # print(len(galaxies_points))
    gaussian_fitting(galaxies_points, img)
    plt.show()
