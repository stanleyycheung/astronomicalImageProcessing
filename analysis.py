import matplotlib.pyplot as plt
from processing import Processor
from matplotlib import colors
import numpy as np
from matplotlib.colors import LogNorm
from scipy import stats

# https://en.wikipedia.org/wiki/Sersic_profile
# test

class Analyzer:
    def __init__(self):
        self.mask = np.load('mask.npy')
        self.testImg = np.load('testData_bigger.npy')
        self.imgMasked = np.load('imgMasked.npy')

        self.background = 3412
        self.sigma = 20
        self.ZP = 25.3
        self.threshold = self.background + 5 * self.sigma

        self.galaxies_points = []
        self.galaxies = []

    def run(self):
        """Run code"""
        self.labelGalaxies(self.testImg)
        self.calcGalaxies(self.testImg, self.digitalMap)
        print(self.galaxies)
        self.test()
        self.plotTestImg(self.testImg)
        # print(self.galaxies_points[-2])

    def labelGalaxies(self, data):
        """Creates a digital map and clusters galaxies together
        Puts the galaxies into galaxy_points"""
        self.digitalMap = data > self.threshold
        self.digitalMap = np.asarray(self.digitalMap, dtype=np.int32)
        print(f'threshold = {self.threshold}')

        for i in range(len(self.digitalMap)):
            for j in range(len(self.digitalMap[i])):
                if self.digitalMap[i, j] == 1:
                    self.floodFill(i, j)

    def calcGalaxies(self, data, digitalMap):
        """Calculates the fluxes from galaxy_points"""
        for galaxy in np.array(self.galaxies_points):
            pixel_count = 0
            x_min, x_max = 1e9, 0
            y_min, y_max = 1e9, 0
            size = 0
            for point in galaxy:
                # print(point, data[point[0], point[1]])
                pixel_count += data[point[0], point[1]]
                if point[0] < x_min:
                    x_min = point[0]
                elif point[0] > x_max:
                    x_max = point[0]
                if point[1] < y_min:
                    y_min = point[1]
                elif point[1] > y_max:
                    y_max = point[1]
                size += 1
            x_mid, y_mid = (x_max + x_min)/2, (y_max + y_min)/2
            # print(x_mid, y_mid)
            # self.findBackground(x_mid, y_mid, (x_max - x_min)/2 * 10, data)
            background = self.findBackground(x_mid, y_mid, 70, data, digitalMap)
            real_count = pixel_count - background * size
            mag_i = -2.5 * np.log10(real_count)
            m = self.ZP + mag_i
            galaxy = {'pos': (x_mid, y_mid), 'm': m, 'size': size, 'real_count': real_count,
                      'total_count': pixel_count, 'background_count': background * size}
            self.galaxies.append(galaxy)

    def findBackground(self, x_mid, y_mid, radius, data, digitalMap, mode=0):
        """Finds the background value by drawing a big circle around star
        Values:
        0 - no star - background
        1 - star but not grouped
        2 - star but grouped
        3 - no star - background but chosen in find background
        """
        xlower = max(int(x_mid - radius), 0)
        xhigher = max(int(x_mid + radius), 0)
        ylower = max(int(y_mid - radius), 0)
        yhigher = max(int(y_mid + radius), 0)
        background_arr = []
        for row in range(xlower, xhigher + 1, 1):
            for column in range(ylower, yhigher + 1, 1):
                try:
                    if radius ** 2 > (row-x_mid)**2 + (column-y_mid)**2 and (digitalMap[row][column] == 0 or digitalMap[row][column] == 3):
                        background_arr.append(data[row][column])
                        digitalMap[row][column] = 3
                        # print(row, column)
                except IndexError:
                    pass
        if len(background_arr) == 0:
            raise ValueError('Something is wrong')
        if mode == 0:
            background = np.mean(background_arr)
            return background
        elif mode == 1:
            fig, ax = plt.subplots()
            plt.hist(background_arr, bins='auto')
            xt = plt.xticks()[0]
            xmin, xmax = min(xt), max(xt)
            lnspc = np.linspace(xmin, xmax, len(background_arr))
            m, s = stats.norm.fit(background_arr)  # get mean and standard deviation
            print(m, s)
            pdf_g = stats.norm.pdf(lnspc, m, s)  # now get theoretical values in our interval
            plt.plot(lnspc, pdf_g, label=f"Norm, mean = {m:.2f}")
            return m

    def floodFill(self, x, y):
        """Calculates points that are clustered together"""
        size_threshold = 2
        object = []
        size = 0
        toFill = set()
        toFill.add((x, y))
        while toFill:
            (x, y) = toFill.pop()
            object.append([x, y])
            self.digitalMap[x, y] = 2
            size += 1
            try:
                if self.digitalMap[x, y+1] == 1:
                    toFill.add((x, y+1))
            except IndexError:
                pass
            try:
                if self.digitalMap[x, y-1] == 1:
                    toFill.add((x, y-1))
            except IndexError:
                pass
            try:
                if self.digitalMap[x+1, y] == 1:
                    toFill.add((x+1, y))
            except IndexError:
                pass
            try:
                if self.digitalMap[x-1, y] == 1:
                    toFill.add((x-1, y))
            except IndexError:
                pass
        if size > size_threshold:
            self.galaxies_points.append(object)
        else:
            for x, y in object:
                self.digitalMap[x, y] = 0

    def test(self):
        fig, ax = plt.subplots()
        cmap = colors.ListedColormap(['black', 'white', 'red', 'yellow'])
        plt.imshow(self.digitalMap, cmap=cmap)

        # plt.hist(data)
        # plt.show()

    def plotTestImg(self, data):
        fig, ax = plt.subplots()
        plt.imshow(data, norm=LogNorm())
        plt.colorbar()


if __name__ == '__main__':
    a = Analyzer()
    a.run()
    plt.show()
    # a.plotTestImg()
