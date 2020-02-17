import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_prominences


def galaxyDistribution(galaxies_points, img):
    for counter in range(len(galaxies_points)):
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
        if x_min == x_max or y_min == y_max:
            np.delete(galaxies_points, counter)
            # print(x_min, x_max, y_min, y_max)
            continue

        y_counts = []
        for i in range(x_min, x_max):
            counts = img[i, y_min:y_max+1]
            y_counts.append(counts)
        y_counts = np.sum(y_counts, axis=0)
        pos_y = np.arange(y_min, y_max+1)

        x_counts = []
        for i in range(y_min, y_max):
            counts = img[x_min:x_max+1, i]
            x_counts.append(counts)
        x_counts = np.sum(x_counts, axis=0)
        pos_x = np.arange(x_min, x_max+1)

        # if len(y_counts) == 0 or len(x_counts) == 0:
        #     raise ValueError(f'{x_min, x_max, y_min, y_max}')

        y_diff = np.sort(np.diff(y_counts))
        x_diff = np.sort(np.diff(x_counts))
        # y_peak, _ = find_peaks(y_counts, prominence=(200,))
        # x_peak, _ = find_peaks(x_counts, prominence=(200,))
        if x_diff[-1] > 1e4 or y_diff[-1] > 1e4 or x_diff[0] < -1e4 or y_diff[0] < -1e4:
            # print('exclude', i)
            # fig, ax = plt.subplots()
            # plt.plot(pos_y, y_counts, '.-')
            # plt.plot(pos_y[y_peak], y_counts[y_peak], 'ro')
            #
            # fig, ax = plt.subplots()
            # plt.plot(pos_x, x_counts, '.-')
            # plt.plot(pos_x[x_peak], x_counts[x_peak], 'ro')
            #
            # fig, ax = plt.subplots()
            # img_plot = img[x_min:x_max, y_min:y_max]
            # img_plot[img_plot == 0] = None
            # plt.imshow(img_plot, origin='upper')
            np.delete(galaxies_points, counter)
            continue
        y_peak, _ = find_peaks(y_counts, prominence=(600,))
        x_peak, _ = find_peaks(x_counts, prominence=(600,))

        if len(x_peak) > 1 or len(y_peak) > 1:
            print(i)
            print(y_diff[-1] > 1e4)
            print(x_diff[-1] > 1e4)
            py = peak_prominences(y_counts, y_peak)[0]
            px = peak_prominences(x_counts, x_peak)[0]

            fig, ax = plt.subplots()
            plt.plot(pos_y, y_counts, '.-')
            plt.title(f'y:{py}')
            plt.plot(pos_y[y_peak], y_counts[y_peak], 'ro')

            fig, ax = plt.subplots()
            plt.plot(pos_x, x_counts, '.-')
            plt.title(f'x:{px}')
            plt.plot(pos_x[x_peak], x_counts[x_peak], 'ro')

            fig, ax = plt.subplots()
            plt.imshow(img[x_min:x_max, y_min:y_max], origin='upper')


if __name__ == '__main__':
    img = np.load('testData_noisy.npy')
    galaxies_points = np.load('galaxies_points.npy', allow_pickle=True)
    test_galaxy = galaxies_points[0]
    print(img.shape)
    print(len(galaxies_points))
    galaxyDistribution(galaxies_points, img)
    plt.show()
