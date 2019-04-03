__author__ = "Navratan Bagwan"
"version = 0.2, Date: 2019-04"
#navratan.bagwan@sund.ku.dk
# solved the MAD issue for different experiment.
# should not fail with the larger data but, have really tested with huge-huge data.
import pdb
import os
import localFDR_median
import all_stats
import psmAssignment
import time
import concurrent.futures
from functools import partial
from optparse import OptionParser

usage = "usage:-e -a -f  -l -p \n Use -h for detailed help for parameter\nversion = 0.1, Date: 2019-02" \
        "\nPTMdesignator is a part (module 2) of SHIFTS \n(Systematic, Hypothesis-free Identification of PTMs with controlled \nFDR " \
        "based on ultra-Tolerant database search)is a program made in the Jesus \nVazquez Cardiovascular Proteomics Lab at Centro Nacional de " \
        "Investigaciones \nCardiovasculares,for highthroughput PTM( Post translation modifications )processing.\n" \
        "PTMdesignator, can be used for FDR filter at local and peak lever\nIn addition this program will also assign the PSMs to peak and tag ORPHANS "

parser = OptionParser(usage=usage)

parser.add_option("-e", "--PathtoExperiment", action="store", type="string", dest="Dir", help="Enter the path to experiment log file")

parser.add_option("-a", "--ApexFile", action="store", type="string", dest="apexfilePath", help="Enter the path to PTM apex file")

parser.add_option("-f", "--fastaFile", action="store", type="string", dest="fasta", help="Enter the fasta file used in CometPTM Search")

parser.add_option("-l", "--Localthreshold", action="store", type="string", dest="LocalFDRs", help="Enter the local FDR threshold; 0.01 for 1% FDR")

parser.add_option("-p", "--Peakthreshold", action="store", type="string", dest="PeakFDRs", help="Enter the Peak FDR threshold; 0.01 for 1% FDR")

(options, args) = parser.parse_args()


if options.Dir:
    Dir = options.Dir
else:
    parser.error("Please enter the path to Experiment LOG file -h")

if options.apexfilePath:
    apexfilePath = options.apexfilePath
else:
    parser.error("Please enter the path to PTm apex file -h")

if options.fasta:
    fasta = options.fasta
else:
    parser.error("Please enter the path to fasta file -h")

if options.LocalFDRs:
    LocalFDRs = options.LocalFDRs
else:
    parser.error("please enter the local fdr threshold -h")

if options.PeakFDRs:
    PeakFDRs = options.PeakFDRs
else:
    parser.error("please enter the Peak fdr threshold -h")


def LocalandPeakFDR(experirmetlist, apexFile, fastafile, localFDRthreshold, peakFDRthreshold):

    start_time = time.time()
    listofFile = []
    fwhmListDic = {}
    for line in open(experirmetlist):
        if line != "\n":
            listofFile.append(line.strip())

            tempPath = os.path.dirname(line.strip()) + "/" + "MAD_and_FWHM_calculations.txt"

            fwhmValue = all_stats.fullWidthHalfMaximum(tempPath)
            if line.strip() not in fwhmListDic:
                fwhmListDic[line.strip()] = fwhmValue


    ApexList = all_stats.createApexList(apexfile=apexFile)


    if __name__ == '__main__':
        #if onlyFirstPart == "1" and OrphanClassTF == "0" and IsotopTF == "0":

        floatList = []
        intList = []
        with concurrent.futures.ProcessPoolExecutor() as executor:
            intList += executor.map(localFDR_median.iDawindwoEstimator, listofFile)
            floatList += executor.map(localFDR_median.iDawindwoEstimator1, listofFile)

            OnDaWindowCenteList = localFDR_median.meadianSlope(intDelatamass=intList, originalDeltaMas=floatList)

        listofDic = []
        func1 = partial(localFDR_median.ClosestDic, listofOneDaltonCenter=OnDaWindowCenteList, ApexList=ApexList)
        with concurrent.futures.ProcessPoolExecutor() as executor:
            listofDic += executor.map(func1, listofFile)


        updatedList = {}
        #print len(listofDic)
        for evryDic in listofDic:
            for firstElement in evryDic:
                if firstElement not in updatedList:
                    updatedList[firstElement]=evryDic[firstElement]
                else:
                    for everySmallList in evryDic[firstElement]:
                        updatedList[firstElement].append(everySmallList)

        updatedList1 = updatedList
        updatedList2 = updatedList

        print("1da window center is calculated")

        localFDR_median.calcualteFDR(mod2Xcor1=updatedList1, outPath=listofFile[0])
        print("local fdr is done")
        AllDataDicforAssignation = localFDR_median.peakClosestDic(Slope_peak_Dict=updatedList2, mad_fwhm =fwhmListDic, filepath=listofFile[0])
        # print len(AllDataList)
        print("peak fdr is done")
        func2 = partial(psmAssignment.taggingMethodForCOMET_PTM, DictofallData=AllDataDicforAssignation, fasta=fastafile, localFDRfilter= localFDRthreshold, peakFDRfilter=peakFDRthreshold)
        with concurrent.futures.ProcessPoolExecutor() as executor:
            executor.map(func2, listofFile)
            #psmAssignment.taggingMethodForCOMET_PTM(DictofallData=AllDataDicforAssignation, calibrationFile=i, fasta=fastafile, localFDRfilter= 0.05, peakFDRfilter=0.05)

            print("PTMesignator has finished")
    print("---%s seconds ---" % (time.time() - start_time))



if __name__ == "__main__":
    LocalandPeakFDR(experirmetlist=Dir, apexFile=apexfilePath, fastafile=fasta, localFDRthreshold=LocalFDRs, peakFDRthreshold=PeakFDRs)
    # LocalandPeakFDR(experirmetlist=r"E:\Users\nbagwan\Desktop\SHIFTS_lastChange\heart\TestEXE1\Peakpicking_Log.txt",
    #             apexFile=r"E:\Users\nbagwan\Desktop\SHIFTS_lastChange\heart\TestEXE1\fakeApex.txt",
    #             fastafile=r"E:\Users\nbagwan\Desktop\SHIFTS_lastChange\ID_uniprot_MusMusculus-concatanated.fasta",
    #             localFDRthreshold= 0.05, peakFDRthreshold= 0.05)
