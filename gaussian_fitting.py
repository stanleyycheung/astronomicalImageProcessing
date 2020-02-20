import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
# unfinished 18/02


def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x, y: height*np.exp(
        -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)


def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y


def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    def errorfunction(p): return np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                          data)
    p, success = optimize.leastsq(errorfunction, params)
    return p


def gaussian_fitting(galaxies_points, img):
    for counter in range(5):
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
        buffer = 0
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
        params = fitgaussian(img_slice)
        fit = gaussian(*params)
        fig, ax = plt.subplots()
        plt.imshow(img_slice, origin='upper')
        plt.contour(fit(*np.indices(img_slice.shape)))
        # ax.contour(x, y, twoD_Gaussian((x, y), *popt).reshape(7, 12), 1, colors='w')


if __name__ == '__main__':
    img = np.load('testData_noisy.npy')
    galaxies_points = np.load('galaxies_points.npy', allow_pickle=True)
    # print(img.shape)
    # print(len(galaxies_points))
    gaussian_fitting(galaxies_points, img)
    plt.show()
