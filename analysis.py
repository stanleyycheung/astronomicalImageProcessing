import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
from matplotlib.colors import LogNorm
from scipy import stats
from pprint import pprint
import time
import json
# https://en.wikipedia.org/wiki/Sersic_profile


class Analyzer:
    def __init__(self):
        self.background = 3412  # mean of gaussian
        self.sigma = 20  # spread
        self.ZP = 25.3  # zeropoint ZP
        self.sigma_num = 3
        self.threshold = self.background + self.sigma_num * self.sigma

        self.galaxies_points = []
        self.galaxies = []

    def load(self, testImgFile, orgImgFile):
        # self.mask = np.load('mask.npy')
        self.testImg = testImgFile  # from maketestdata
        self.orgImg = orgImgFile
        # self.imgMasked = np.load('imgMasked.npy')
        # load blooming data

    def run(self):
        """Run code"""
        t1 = time.process_time()
        self.labelGalaxies(self.testImg)
        t2 = time.process_time()
        print(f'Took {t2-t1:.2f}s to label galaxies')
        self.calcGalaxies(self.testImg, self.digitalMap)  # digital map: creating catalog
        t3 = time.process_time()
        print(f'Took {t3-t2:.2f}s to calculate galaxies')
        # pprint(self.galaxies)
        print(f'No. of galaxies = {len(self.galaxies)}')
        self.write()
        self.plotDigital()  # plot digital graph
        self.plotData(self.testImg)
        np.save('galaxies_points', self.galaxies_points)

    def labelGalaxies(self, data):
        """Creates a digital map and clusters galaxies together
        Puts the galaxies into galaxy_points"""  # above threshold = 1, rest =0
        self.digitalMap = data > self.threshold  # in True or False
        self.digitalMap = np.asarray(self.digitalMap, dtype=np.int32)  # Convert into 0 or 1
        print(f'threshold = {self.threshold}')

        for i in range(len(self.digitalMap)):
            for j in range(len(self.digitalMap[i])):
                if self.digitalMap[i, j] == 1:
                    self.floodFill(i, j)

    def calcGalaxies(self, data, digitalMap):
        """Calculates the fluxes from galaxy_points"""
        interval = len(self.galaxies_points)//100
        differences = []
        self.clip_counter = 0
        for i, galaxy in enumerate(np.array(self.galaxies_points)):
            pixel_count = 0
            x_min, x_max = 1e9, 0
            y_min, y_max = 1e9, 0
            size = 0
            for point in galaxy:
                # print(point, data[point[0], point[1]])
                pixel_count += data[point[0], point[1]]  # each are x and y coordinate of the point
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
            radius = max((x_max-x_mid), (x_mid-x_min), (y_max-y_mid), (y_mid-y_min))

            # print(x_mid, y_mid)
            # self.findBackground(x_mid, y_mid, (x_max - x_min)/2 * 10, data)
            if self.clipper(x_mid, y_mid, radius, differences, data):
                background = self.findBackground(x_mid, y_mid, 70, data, digitalMap)
                real_count = pixel_count - background * size
                mag_i = -2.5 * np.log10(real_count)
                m = self.ZP + mag_i
                galaxy_dict = {'pos': (x_mid, y_mid), 'm': m, 'size': size, 'real_count': real_count,
                               'total_count': pixel_count, 'background_count': background * size}
                self.galaxies.append(galaxy_dict)
            if i % interval == 0:
                print(f'Galaxy {i} out of {len(self.galaxies_points)}')
        # access differences here
        fig, ax = plt.subplots()
        plt.hist(differences, bins='auto')
        plt.xlabel('difference in number of counts')
        plt.ylabel('number of entries')
        plt.title('difference plot - find threshold')

        print("Number of clipped galaxies = ", self.clip_counter)

    def clipper(self, x_mid, y_mid, radius, differences, data):
        # patches = np.load('clipper.npy')
        # circle = patches[0]
        # square = patches[1]

        # method 1: compare distance to blooming region and radius
        # radius = np.sqrt(size/np.pi)
        xlower = max(int(x_mid - radius), 0)
        ylower = max(int(y_mid - radius), 0)
        xhigher = int(x_mid + radius)
        yhigher = int(y_mid + radius)

        orig_count = []
        mask_count = []  # try to eliminate these two later on
        for row in range(xlower, xhigher + 1, 1):
            for column in range(ylower, yhigher + 1, 1):
                try:  # make circular aperture out of it
                    if radius ** 2 > (row-x_mid)**2 + (column-y_mid)**2:
                        orig_count.append(self.orgImg[row][column])
                        mask_count.append(data[row][column])
                        # print(row, column)
                except IndexError:
                    pass

        orig_tot = np.sum(orig_count)
        mask_tot = np.sum(mask_count)
        differences.append(abs(orig_tot-mask_tot))
        if abs(orig_tot-mask_tot) == 0:
            return True
        else:
            print('cut galaxy has midpoints', x_mid, ',', y_mid, 'radius is', radius)
            self.clip_counter += 1
            return False

    def findBackground(self, x_mid, y_mid, radius, data, digitalMap, mode=0):
        """Finds the background value by drawing a big circle around star"""
        xlower = max(int(x_mid - radius), 0)
        xhigher = max(int(x_mid + radius), 0)
        ylower = max(int(y_mid - radius), 0)
        yhigher = max(int(y_mid + radius), 0)  # make square
        background_arr = []
        for row in range(xlower, xhigher + 1, 1):
            for column in range(ylower, yhigher + 1, 1):
                try:  # make circular aperture out of it
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
        object = []  # list of points of one galaxy detected
        size = 0
        toFill = set()
        toFill.add((x, y))
        while toFill:
            (x, y) = toFill.pop()
            if x < 0 or y < 0:
                continue
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
            self.galaxies_points.append(object)  # arrays in array: list of objects
        else:
            for x, y in object:
                self.digitalMap[x, y] = 0

    def write(self):
        with open('galaxies.json', 'w') as fout:
            json.dump(self.galaxies, fout)

    def plotDigital(self):
        fig, ax = plt.subplots()
        cmap = colors.ListedColormap(['black', 'white', 'red', 'yellow'])
        plt.imshow(self.digitalMap, cmap=cmap, origin='lower')

        # plt.hist(data)
        # plt.show()

    def plotData(self, data):
        fig, ax = plt.subplots()
        plt.imshow(data, norm=LogNorm(), origin='lower')
        plt.colorbar()


if __name__ == '__main__':
    a = Analyzer()
    a.load(np.load('testData_noisy.npy'), np.load('orgData_noisy.npy'))
    a.run()
    plt.show()
    # a.plotData()
