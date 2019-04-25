__author__ = "Navratan Bagwan"
from scipy.stats.stats import linregress
import os
import itertools
import pdb
import pandas as pd
#alpha = all_stats.getAlpha(alphafile)
def iDawindwoEstimator(processingFileList):
    alpha = 1
    slope_intercept = []
    with open(processingFileList) as PTMfile:
        next(PTMfile)
        for line3 in PTMfile:
            if line3 != "\n":
                splits3 = line3.split("\t")
                calibrated_Delta_MH = float(splits3[20].strip())

                if calibrated_Delta_MH > 0:
                    intOFmass = int(calibrated_Delta_MH + 0.5)
                    slope_intercept.append(intOFmass)
                else:
                    intOFmass = int(calibrated_Delta_MH - 0.5)
                    slope_intercept.append(intOFmass)

        return slope_intercept

def iDawindwoEstimator1(processingFileList):
    slopeMassList = []
    with open(processingFileList) as PTMfile:
        next(PTMfile)
        for line3 in PTMfile:
            if line3 != "\n":
                splits3 = line3.split("\t")
                calibrated_Delta_MH = float(splits3[20].strip())
                slopeMassList.append(calibrated_Delta_MH)
    return slopeMassList



    # minOFslope = min(slope_intercept)
    # maxOFslope = max(slope_intercept)
    # listofPsuedoNumber = range(minOFslope, maxOFslope)
    #
    # slopeValue, intercept, r, p, std_erro = linregress(slope_intercept, slopeMassList)
    #
    # for everyPsuedo in listofPsuedoNumber:
    #     medianPoint = everyPsuedo * slopeValue + intercept
    #     slopeMedian.append(medianPoint)

def meadianSlope(intDelatamass, originalDeltaMas):
    localCenter = []
    flat1int = itertools.chain.from_iterable(intDelatamass)
    flat2float = itertools.chain.from_iterable(originalDeltaMas)
    newListint = list(flat1int)
    newListfloat = list(flat2float)
    minOFslope = min(newListint)
    maxOFslope = max(newListint)

    listofPsuedoNumber = range(minOFslope, maxOFslope)
    slopeValue, intercept, r, p, std_erro = linregress(newListint, newListfloat)

    for everyPsuedo in listofPsuedoNumber:
        medianPoint = everyPsuedo * slopeValue + intercept
        localCenter.append(medianPoint)


    return localCenter

def ClosestDic(processingFileList, listofOneDaltonCenter, ApexList):
    takeClosest = lambda num, collection: min(collection, key=lambda x: abs(x - num))
    mod2Xcor = {}

    with open(processingFileList) as PTMfile:
        next(PTMfile)
        for newline in PTMfile:
            if newline != "\n":
                splitsnew = newline.split("\t")
                label = splitsnew[18].strip()
                delta_modification = float(splitsnew[20].strip())
                xscore = float(splitsnew[16].strip())

                if delta_modification > -502.00 and delta_modification < 502.00:

                    closestSlope = takeClosest(delta_modification, listofOneDaltonCenter)
                    apexClosest = takeClosest(delta_modification,ApexList)

                    if closestSlope not in mod2Xcor:
                        mod2Xcor[closestSlope] = [[xscore, str(label).strip(), newline.strip(), delta_modification, apexClosest, str(processingFileList)]]
                    else:
                        mod2Xcor[closestSlope].append([xscore, str(label).strip(), newline.strip(), delta_modification, apexClosest, str(processingFileList)])

    return mod2Xcor


def calcualteFDR(mod2Xcor1,outPath):
    alpha = 1
    firstfilepath = os.path.dirname(outPath)
    fileName = firstfilepath + "/SlopeFDRfile.txt"
    w = open(fileName, "w")
    mainList = ["Scan\tSearchengineRank\tcharge\texpMH\ttheoMH\tExpmz\tXcor\tSeq\tRetentionTime\tProtAcc\tDeltaMod"
                "\tB-series\tY.series\tIsotpicJump\tdeltaPeptide\tfilename\tCorXcor\tnewExpMH\tLabel\tMedian"
                "\tCal_Delta_MH\tCal_Delta_M/Z\tCalExp_MZ\t1DawindowCenter\tdecoyCount\ttargetCount\tlocalFDR\t1Da_start\t1Da_end\tPeakApex\tExperiment Name\n"]
    for mod in mod2Xcor1:
        decoy, target = 0.0, 0.0
        mod2Xcor1[mod].sort(key=lambda row: (row[0]), reverse=True)
        MinMaxList = [item[3] for item in mod2Xcor1[mod]]
        startingPoint = min(MinMaxList)
        endingPoint = max(MinMaxList)
        for eachline in mod2Xcor1[mod]:
            if eachline[1] == "Decoy":
                decoy += 1.0
                eachline.append(decoy * float(alpha))
                eachline.append(target)
                try:
                    localFDR = float((decoy * float(alpha)) / target)
                except ZeroDivisionError:
                    localFDR = "1"
                eachline.append(str(localFDR).strip())
                eachline.append(str(startingPoint))
                eachline.append(str(endingPoint))

            elif eachline[1] == "Target":
                target += 1.0
                eachline.append(decoy * float(alpha))
                eachline.append(target)
                try:
                    localFDR = float((decoy * float(alpha)) / target)
                except ZeroDivisionError:
                    localFDR = "1"
                eachline.append(str(localFDR))
                eachline.append(str(startingPoint))
                eachline.append(str(endingPoint))

    for fdrline in mod2Xcor1:
        for everyLine in mod2Xcor1[fdrline]:

            mainList.append("\t".join(everyLine[2].split("\t")) + "\t" + str(fdrline) + "\t" + str(
                everyLine[6]).strip() + "\t" + str(
                everyLine[7]).strip() + "\t" + str(everyLine[8]).strip() + "\t" + str(
                everyLine[9]).strip() + "\t" + str(everyLine[10]).strip() + "\t" + str(everyLine[4]).strip() + "\t" + str(everyLine[5]) + "\n")

    for newline in mainList:
        w.write(str(newline))
    w.close()


def peakClosestDic(Slope_peak_Dict,mad_fwhm,filepath):

    DiffExpData = {}
    outputFilelist = [
        "Scan\tSearchengineRank\tcharge\texpMH\ttheoMH\tExpmz\tXcor\tSeq\tRetentionTime\tProtAcc\tDeltaMod"
        "\tB-series\tY.series\tIsotpicJump\tdeltaPeptide\tfilename\tCorXcor\tnewExpMH\tLabel\tMedian"
        "\tCal_Delta_MH\tCal_Delta_M/Z\tCalExp_MZ\tPeakApex\tpeakDecoyCount\tpeakTargetCount\tPeakFDR\tLocalFDR\tExperimentName\n"]
    # outputFilelist = []
    outputfilename = os.path.dirname(filepath) + "/Peak_and_Slope_FDRfile.txt"
    w = open(outputfilename, "w")
    slopePeak_list = Slope_peak_Dict.keys()
    ApexList = []
    outputFilelist_dic = {}

    takeClosest = lambda num, collection: min(collection, key=lambda x: abs(x - num))
    for closeinSlopePeak_List in Slope_peak_Dict:
        for element in Slope_peak_Dict[closeinSlopePeak_List]:

            everyApex = element[4]
            theo_mass = float(everyApex)
            experimental_mass = float(element[2].split("\t")[20])
            corrXcor = element[2].split("\t")[16]
            label = element[2].split("\t")[18]
            charge = element[2].split("\t")[2]

#            ppm_error = abs((experimental_mass - theo_mass) * 1000000 / theo_mass)
            massDiff = abs(experimental_mass - theo_mass)
            # if ppm_error <= 10:
            if massDiff <= float(charge) * float(mad_fwhm[element[5]]):  ##### checking if a scan falls within a PTM peak by using MAd and Charge
                # if massDiff <= float(mad_fwhm):
                # outputFilelist.append(element)

                if everyApex not in outputFilelist_dic:
                    outputFilelist_dic[everyApex] = [[float(corrXcor), label, element[2], str(everyApex), element[8], str(element[5])]]
                else:
                    outputFilelist_dic[everyApex].append([float(corrXcor), label, element[2], str(everyApex), element[8], str(element[5])])

                Slope_peak_Dict[closeinSlopePeak_List].remove(element)

    #return outputFilelist_dic
    for apex_mass in outputFilelist_dic:
        outputFilelist_dic[apex_mass].sort(key=lambda row: (row[0]), reverse=True)
        decoy, target = 0.0, 0.0
        for everyline in outputFilelist_dic[apex_mass]:
            expName = everyline[-1]
            if everyline[1] == "Decoy":
                decoy += 1
                everyline.append(decoy)
                everyline.append(target)
                try:
                    PeakFDR = float(decoy / target)
                except ZeroDivisionError:
                    PeakFDR = "1"
                everyline.append(PeakFDR)
                outputFilelist.append("\t".join(everyline[2].split("\t")) + "\t" + str(everyline[3]) + "\t"+ str(decoy) + "\t"+ str(target) + "\t" + str(PeakFDR) + "\t" + str(everyline[4]) + "\t" + str(everyline[5]) +"\n" )

            elif everyline[1] == "Target":
                target += 1
                everyline.append(decoy)
                everyline.append(target)
                try:
                    PeakFDR = float(decoy / target)
                except ZeroDivisionError:
                    PeakFDR = "1"
                everyline.append(PeakFDR)


                outputFilelist.append("\t".join(everyline[2].split("\t")) + "\t" + str(everyline[3]) + "\t"+ str(decoy) + "\t" + str(target) + "\t" + str(PeakFDR) + "\t" + str(everyline[4]) + "\t" + str(everyline[5]) + "\n")

                ####### creating outpu return for for final assignation

                # if expName not in DiffExpData:
                #     DiffExpData[expName] = ["\t".join(everyline[2].split("\t")) + "\t" + str(everyline[3]) + "\t"+ str(decoy) + "\t"+ str(target) + "\t" + str(PeakFDR) + "\t" + str(everyline[4]) + "\t" + str(everyline[5]) +"\n" ]
                # else:
                #     DiffExpData[expName].append("\t".join(everyline[2].split("\t")) + "\t" + str(everyline[3]) + "\t"+ str(decoy) + "\t"+ str(target) + "\t" + str(PeakFDR) + "\t" + str(everyline[4]) + "\t" + str(everyline[5]) +"\n" )

    for restMass in Slope_peak_Dict:
        for mass in Slope_peak_Dict[restMass]:
            outputFilelist.append("\t".join(mass[2].split("\t")) + "\t" + "Orphan" + "\t" + "1"+ "\t" + "1" + "\t" + "1" + "\t" + str(mass[8]) + "\t" + str(mass[5]) + "\n")
            # if mass[5] not in DiffExpData:
            #     DiffExpData[mass[5]] = ["\t".join(mass[2].split("\t")) + "\t" + "Orphan" + "\t" + "1" + "\t" + "1" + "\t" + "1" + "\t" + str(mass[8]) + "\t" + str(mass[5]) + "\n"]
            # else:
            #     DiffExpData[mass[5]].append("\t".join(mass[2].split("\t")) + "\t" + "Orphan" + "\t" + "1" + "\t" + "1" + "\t" + "1" + "\t" + str(mass[8]) + "\t" + str(mass[5]) + "\n")


    # df = pd.DataFrame(outputFilelist)

    # for test in DiffExpData:
    #     print(len(DiffExpData[test]))


    for towrite in outputFilelist:
        w.write(str(towrite))
    w.close()

    return outputfilename
    # return DiffExpData

def ProduceFileforEveryExperiment(listofpath,bigDataFile):
    outputFilelist = [
        "Scan\tSearchengineRank\tcharge\texpMH\ttheoMH\tExpmz\tXcor\tSeq\tRetentionTime\tProtAcc\tDeltaMod"
        "\tB-series\tY.series\tIsotpicJump\tdeltaPeptide\tfilename\tCorXcor\tnewExpMH\tLabel\tMedian"
        "\tCal_Delta_MH\tCal_Delta_M/Z\tCalExp_MZ\tPeakApex\tpeakDecoyCount\tpeakTargetCount\tPeakFDR\tLocalFDR\tExperimentName\n"]
    outputfilename = os.path.dirname(listofpath) + "/Peak_and_Slope_FDRfile-Smallchunk.txt"
    w = open(outputfilename, "w")
    with open(bigDataFile) as bigF:
        next(bigF)
        for line in bigF:
            ExpName = line.split("\t")[-1].strip()
            if listofpath == ExpName:
                w.write(line)


    w.close()
