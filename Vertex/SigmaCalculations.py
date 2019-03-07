__author__ = "Navratan Bagwan"
import numpy as np
from scipy.optimize import leastsq
import pdb
import os
import shutil
import matplotlib
matplotlib.use('agg')
import all_stats
from pylab import *
separations = "/"

def SigmaCalculation(deltaPPMlist, processingFileList):
    outfile = []
    bins = np.arange(-20, 20 + .5, .5)
    ydata, xdata = np.histogram(deltaPPMlist, bins)
    xdata = xdata.astype(float).tolist()
    ydata = ydata.astype(float).tolist()
    xdata.pop()
    fitfunc = lambda p, x: p[0] * exp(-0.5 * ((x - p[1]) / p[2]) ** 2) + p[3]
    errfunc = lambda p, x, y: (y - fitfunc(p, x))
    init = [1.0, 0.5, 0.5, 0.5]
    out = leastsq(errfunc, init, args=(xdata, ydata))
    # out = leastsq(errfunc, init, args=(aa,bb))
    c = out[0]

    outfile.append(
        "A exp[-0.5((x-mu)/sigma)^2] + k" + "\n" + "Parent Coefficients:" + "\n" + "1.000, 0.200, 0.300, 0.625" + "\n" + "Fit Coefficients:" + str(
            c[0]) + "\t" + str(c[1]) + "\t" + str(abs(c[2])) + "\t" + str(c[3]) + "\n" + "threeSigma: " + str(
            abs(c[2] * 3)))
    firstfilepath = os.path.dirname(processingFileList[0])
    fileName = firstfilepath + "/sigmaCalculations.txt"
    w = open(fileName, "w")
    for i in outfile:
        w.writelines(i)
    w.close()

    for eachfile in processingFileList[1:len(processingFileList)]:
        foldername = os.path.dirname(eachfile)
        shutil.copy(fileName, foldername)
