import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# unfinished 18/02
# hi


def twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    (x, y) = xdata_tuple
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp(- (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                                     + c*((y-yo)**2)))
    return g.ravel()


def gaussian_fitting(galaxies_points, img):
    for counter in range(10):
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
            print(f'weird stuff {counter}')
            continue

        x = np.arange(x_min, x_max, 1)
        y = np.arange(y_min, y_max, 1)
        x, y = np.meshgrid(y, x)
        img_slice = img[x_min:x_max, y_min:y_max]
        popt, pcov = curve_fit(twoD_Gaussian, (x, y), img_slice.ravel(),
                               p0=[100, 200, 6, 1, 1, 0, 3000])
        print(popt)
        fig, ax = plt.subplots()
        plt.imshow(img_slice, origin='upper')
        # ax.contour(x, y, twoD_Gaussian((x, y), *popt).reshape(7, 12), 1, colors='w')


if __name__ == '__main__':
    img = np.load('testData_noisy.npy')
    galaxies_points = np.load('galaxies_points.npy', allow_pickle=True)
    # print(img.shape)
    # print(len(galaxies_points))
    gaussian_fitting(galaxies_points, img)
    plt.show()
