from astropy.io import fits
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
from pprint import pprint
import seaborn as sns

'''
To do:

Stanley
- exclude grey part of data

Hanna
- exclude bright spot
- write function to create fits file
'''

pixelData = np.load('pixelData.npy')
hdulist = fits.open('A1_mosaic.fits')


class Processor:
    def __init__(self):
        self.img = pixelData
        self.mask = np.full(self.img.shape, True)

    def plotBackground(self, data):
        xlower = 3340
        xupper = 3550
        # xlower, xupper = 3400, 3440
        background = [i for i in data.flatten() if xlower <= i <= xupper]
        plt.hist(background, bins=3550-3340+1, density=True)

        print(np.mean(background))
        print(np.std(background, ddof=1))

        xt = plt.xticks()[0]
        xmin, xmax = min(xt), max(xt)
        lnspc = np.linspace(xmin, xmax, len(background))
        m, s = stats.norm.fit(background)  # get mean and standard deviation
        print(m, s)
        pdf_g = stats.norm.pdf(lnspc, m, s)  # now get theoretical values in our interval
        plt.plot(lnspc, pdf_g, label=f"Norm, mean = {m:.2f}")  # plot it
        plt.legend(loc=1)
        plt.savefig('background_aftercorrection.png')

    def fillMask(self):
        self.edgeExclude()
        self.starsExclude()
        for i in range(len(self.img)):
            for j in range(len(self.img[0])):
                self.greyAreaExclude(i, j)
                if self.img[i][j] >= 48000:
                    self.mask[i][j] = False
                else:
                    continue
        np.save('mask', self.mask)

    def greyAreaExclude_old(self, img):
        def helper(xstart, ystart, xdir, ydir):
            running = True
            currentX = xstart
            currentY = ystart
            while running:
                if self.img[currentX][currentY] == self.img[currentX][currentY + ydir]:
                    self.mask[currentX][currentY] = False
                    currentY = currentY + ydir
                else:
                    currentY = ystart
                    if self.img[currentX][currentY] == self.img[currentX + xdir][currentY]:
                        currentX = currentX + xdir
                    else:
                        break
                print(currentX, currentY)
        helper(1, 1, 1, 1)
        # helper(img.shape[0] - 2, 2, -1, 1)
        # helper(2, img.shape[1] - 2, 1, -1)
        # helper(img.shape[0] - 2, img.shape[1] - 2, -1, -1)
        return self.mask

    def greyAreaExclude(self, i, j):
        # need to modify
        if self.img[i][j] == 3421:
            try:
                if self.img[i][j + 1] == 3421 and self.img[i + 1][j] == 3421 and self.img[i - 1][j] == 3421 and self.img[i][j - 1] == 3421:
                    # self.mask[i + 1][j + 1] = False
                    # self.mask[i - 1][j - 1] = False
                    # self.mask[i + 1][j - 1] = False
                    # self.mask[i - 1][j + 1] = False
                    self.mask[i][j] = False
                    self.mask[i][j + 1] = False
                    self.mask[i + 1][j] = False
                    self.mask[i - 1][j] = False
                    self.mask[i][j-1] = False
            except IndexError:
                pass

    def edgeExclude(self):
        self.mask[0] = np.full(self.mask[0].shape, False)
        # self.mask[1] = np.full(self.mask[0].shape, False)
        self.mask[:, 0] = np.full(len(self.mask), False)
        self.mask[:, len(self.mask[0])-1] = np.full(len(self.mask), False)
        self.mask[:, len(self.mask[0])-2] = np.full(len(self.mask), False)
        self.mask[len(self.mask) - 1] = np.full(self.mask[0].shape, False)
        self.mask[len(self.mask) - 2] = np.full(self.mask[0].shape, False)

    def bloomingExclude(self, i, j):
        # i, j left pointer
        # loop throught image only once -- call from fillmask
        threshold = 30000
        try:
            if self.img[i][j+1]-self.img[i][j] > threshold:
                # print(i, j)

                pass
        except IndexError:
            pass

    def starsExclude(self):
        def excluder(xlower, ylower, xhigher, yhigher):
            for row in range(xlower, xhigher + 1, 1):
                for column in range(ylower, yhigher + 1, 1):
                    # print(row, column)
                    self.mask[column][row] = False

        def excluder_circle(center, radius):
            xlower = int(center[0] - radius)
            xhigher = int(center[0] + radius)
            ylower = int(center[1] - radius)
            yhigher = int(center[1] + radius)
            for row in range(xlower, xhigher + 1, 1):
                for column in range(ylower, yhigher + 1, 1):
                    if radius ** 2 > (row-center[0])**2 + (column-center[1])**2:
                        self.mask[column][row] = False
        # excluder(735, 3284, 820, 3367)  # left up star
        # excluder(1248, 3046, 1623, 3410)  # center big star
        # excluder(2112, 3736, 2159, 3787)  # right up star
        # excluder(946, 2741, 1008, 2814)  # left center star
        # excluder(876, 2259, 942, 2321)  # left center star

        excluder_circle([(820+735)/2, (3367+3284)/2], (820-735)/2 + 7)
        excluder_circle([(1248+1623)/2, (3046+3410)/2], (3410-3046)/2 + 10)
        excluder_circle([(2112+2159)/2, (3787+3736)/2], (2159-2112)/2 + 7)
        excluder_circle([(946+1008)/2, (2814+2741)/2], (1008-946)/2 + 7)
        excluder_circle([(876+942)/2, (2321+2259)/2], (942-876)/2 + 7)

        excluder(1420, 5, 1450, 4604)  # large blooming line

        excluder(1027, 425, 1045, 452)  # small cross
        excluder(1639, 333, 1650, 356)  # small L

        excluder(1100, 426, 1649, 233)
        excluder(1295, 434, 1544, 459)
        # 2nd bloom from top
        excluder(1020, 314, 1701, 330)
        excluder(1382, 330, 1511, 370)
        # 3rd bloom
        excluder(1389, 217, 1475, 264)

    def createFile(self):
        self.imgMasked = np.zeros(self.img.shape)
        self.imgMasked[self.mask] = self.img[self.mask]
        print(self.imgMasked.shape)
        hdu = fits.PrimaryHDU(self.imgMasked)
        hdul = fits.HDUList([hdu])
        hdul.writeto('masked1.fits')

    def createMask(self):
        # p.plotBackground()
        p.fillMask()
        p.createFile()

    def readMask(self):
        self.mask = np.load('mask.npy')
        self.imgMasked = np.zeros(self.img.shape)
        self.imgMasked[self.mask] = self.img[self.mask]
        np.save('imgMasked', self.imgMasked)


if __name__ == '__main__':
    p = Processor()
    p.createMask()
    # p.readMask()
    # p.fitBackground(p.imgMasked)
    # p.plotBackground(p.imgMasked)
    # p.plotBackground(p.img)
    plt.show()

    # for i in np.linspace(0, len(p.img)-1, num=10):
    #     row = p.img[int(i//1)]
    #     plt.plot(np.arange(0, len(row)), row, '.-')
    # plt.show()
