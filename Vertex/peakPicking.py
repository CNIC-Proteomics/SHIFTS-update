__author__ = "Navratan Bagwan"
import os
import sys
import pdb
from scipy.stats.stats import linregress
import shutil
import numpy as np
import all_stats
import pandas
import time
import warnings
import concurrent.futures
from functools import partial
import itertools

separations = "/"


def target_slope(filelist, target_dic, bintoHist, window, rerunTF, extensionName ,target=True):
    firstFiller, startPoint, LastFiller = all_stats.listFiller(window=window)


    outputlist = [
        ["Bin", "\t", "Frequency", "\t", "Slope1", "\t", "Slope2", "\t", "peak-Width", "\t", "peak-Apex", "\t",
         "intercept_mass", "\n"]]
    separations = "/"

    firstfilepath = os.path.dirname(filelist[0])
    if rerunTF!= "0":
        intialFilename = "/" + str(extensionName)+ "_" + "target_Peak_identification_histogram.txt"
        if target:
            fileName = firstfilepath + intialFilename
        else:
            fileName = firstfilepath + "/Decoy_Peak_identification_histogram.txt"
    else:
        intialFilename = "/target_Peak_identification_histogram.txt"
        if target:
            fileName = firstfilepath + intialFilename
    #


    w1 = open(fileName, "w")
    dictlist = []

    start_time = time.time()
    for key, value in target_dic.items():
        temp = [key, value]
        dictlist.append(temp)

    sorted_dictlist = sorted(dictlist, key=lambda x: float(x[0]))
    df = pandas.DataFrame.from_records(sorted_dictlist)
    binlist = df[0].astype(float).tolist()
    frelist = df[1].astype(int).tolist()

    slope1 = [0] * firstFiller
    deltamassItrator = all_stats.window(binlist, int(window))
    FrequencyItrator = all_stats.window(frelist, int(window))
    for Dm, freqeuncey in zip(deltamassItrator, FrequencyItrator):
        s, intercept, r, p, std_error = linregress(Dm, freqeuncey)
        slope1.append(s)

    slope1 = slope1 + [0] * LastFiller

    slope2 = [0] * firstFiller
    deltamassItrator = all_stats.window(binlist, int(window))
    slope2Itraotor = all_stats.window(slope1, int(window))
    for Dm1, gaussFirst in zip(deltamassItrator, slope2Itraotor):
        s1, intercept1, r1, p1, std_error1 = linregress(Dm1, gaussFirst)
        slope2.append(s1)
    slope2 = slope2 + [0] * LastFiller

    ApexInterCept = [0]
    deltamassItrator = all_stats.window(binlist, 2)
    slope2Itraotor = all_stats.window(slope1, 2)

    for x, y in zip(deltamassItrator, slope2Itraotor):

        if sum(y) == 0.0:
            intercept2 = 0
            ApexInterCept.append(intercept2)
        else:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                s2, intercept2, r2, p2, std_error2 = linregress(y, x)
                ApexInterCept.append(intercept2)

    Apex = [0]
    for x1, y1 in zip(slope1[::], slope1[1::]):
        if x1 >= 0.0 and y1 < 0.0:
            Apex.append(1)
        else:
            Apex.append(0)

    width = []
    [width.append(0) if x2 > 0.0 else width.append(1) for x2 in slope2]

    ######## part for second derivative filer input ##########
    # Slope2Threshold = [0]
    # deltamassItrator = all_stats.window(binlist, 2)
    # slopeDeterminator = all_stats.window(slope2, 2)
    # for x2, y2 in zip(deltamassItrator, slopeDeterminator):
    #     if sum(y2) == 0.0:
    #         intercept3 = 0
    #         Slope2Threshold.append(intercept3)
    #     else:
    #         with warnings.catch_warnings():
    #             warnings.simplefilter("ignore", category=RuntimeWarning)
    #             s3, intercept3, r3, p3, std_error3 = linregress(y2, x2)
    #             Slope2Threshold.append(intercept3)

    # print len(slope1)
    # print len(slope2)
    # print len(ApexInterCept)
    # print len(Apex)
    # print len(width)

    outTest = pandas.DataFrame(np.column_stack([binlist, frelist, slope1, slope2, width, Apex, ApexInterCept]),
                               columns=["bin", "freq", "first derivative", "second derivative", "gauss width of peak",
                                        "apex", "intercept(final Apex)"])

    outTest.to_csv(fileName, index=False, sep="\t")

    for eachfile in filelist[1:len(filelist)]:
        foldername = os.path.dirname(eachfile)
        shutil.copy(fileName, foldername)

    print("---%s seconds ---" % (time.time() - start_time))
    return fileName


def SlopeCalculation(filelist, bintoHist):
    target_histDic = {}
    decoy_histDic = {}
    dictlist = []
    #for file in filelist:
    with open(filelist) as PTMfile:
        next(PTMfile)
        for line3 in PTMfile:
            if line3 != "\n":
                splits3 = line3.split("\t")
                charge = int(splits3[2].strip())
                corXcor = float(splits3[16].strip())
                Delta_bin = float(bintoHist)
                calibrated_Delta_MH = float(splits3[20].strip())
                label = str(splits3[18].strip())

                if calibrated_Delta_MH > 0:
                    intformula = int(
                        calibrated_Delta_MH / Delta_bin) * Delta_bin + Delta_bin / 2  #### intformula makes the masses to center to 0.5 values using bins
                    # truncateDmass = float("%.3f" % intformula)
                    truncateDmass = all_stats.check(intformula, 4)
                    truncateDmass1 = truncateDmass

                else:
                    intformula = int(calibrated_Delta_MH / Delta_bin) * Delta_bin - Delta_bin / 2
                    # truncateDmass = float("%.3f" % intformula)
                    truncateDmass = all_stats.check(intformula, 4)
                    truncateDmass1 = truncateDmass

                if label == "Target":
                    if truncateDmass1 not in target_histDic:
                        target_histDic[truncateDmass1] = 1
                    else:
                        target_histDic[truncateDmass1] += 1
                # else:
                #     if truncateDmass1 not in decoy_histDic:
                #         decoy_histDic[truncateDmass1] = 1
                #     else:
                #         decoy_histDic[truncateDmass1] += 1

    return target_histDic


def callSlopeMethod(filelist, bintoHist, windowSize, rerun, extensionName):
    # if __name__ == '__main__':
    targetDiclist = []
    func1 = partial(SlopeCalculation, bintoHist=bintoHist)
    listfile = []
    if rerun != "0":
        for line in open(filelist):
            listfile.append(line.strip())
        filelist = listfile
    else:
        filelist = filelist

    with concurrent.futures.ProcessPoolExecutor() as executor:
        #targetDic.update(executor.map(func1, filelist))
        targetDiclist += executor.map(func1, filelist)

    #### finally two lists are created, one for target and one decoys and both are than send to target_slope method for guassian modelling ####
        useDic = {}
        for i in targetDiclist:
            for dm in i:
                if dm not in useDic:

                    useDic[dm] = i[dm]
                else:
                    useDic[dm] = useDic[dm] + i[dm]

        targetFileName = ["", ""]
        targetFileName[0] = target_slope(filelist, target_dic=useDic, bintoHist=bintoHist, target=True, window=windowSize, rerunTF=rerun, extensionName=extensionName)
    # targetFileName[1] = target_slope(filelist, target_dic=decoy_histDic, bintoHist=bintoHist, target=False,
    #                                  window=windowSize)
    # return targetFileName


