__author__ = "Navratan Bagwan"
"version = 0.1, Date: 2019-02"

import pdb
import all_stats
import time
import OrphanClassifier
import Isotop_CorrectionAdvance
import concurrent.futures
from functools import partial
from optparse import OptionParser

usage = "usage:-e -s -m  -O -I -o -i \n Use -h for detailed help for parameter\nversion = 0.1, Date: 2019-02" \
        "\nPTMcorrector is a part (module 4) of SHIFTS \n(Systematic, Hypothesis-free Identification of PTMs with controlled \nFDR " \
        "based on ultra-Tolerant database search)is a program made in the Jesus \nVazquez Cardiovascular Proteomics Lab at Centro Nacional de " \
        "Investigaciones \nCardiovasculares,for highthroughput PTM( Post translation modifications )processing.\n" \
        "PTMcorrector, can be used for Orphan calissifiction and isotop correction of PTM data"

parser = OptionParser(usage=usage)

parser.add_option("-e", "--PathtoExperiment", action="store", type="string", dest="Dir", help="Enter the path to experiment log file")

parser.add_option("-s", "--sigmafile", action="store", type="string", dest="sigmafilename", help="Enter the name of Sigma log file")

parser.add_option("-m", "--madfile", action="store", type="string", dest="MADfilename", help="Enter the name of MAD log file")

parser.add_option("-O", "--OrphanOUT", action="store", type="string", dest="orphanOUTname", help="Enter the output file name for Orphan classificaition")

parser.add_option("-I", "--targetOfiso", action="store", type="string", dest="Isotargetfile", help="enter the name of target file for isotop correction")

parser.add_option("-i", "--isotopTF", action = "store", type= "string", dest= "isoTF", help="1 to perform istop correction, 0 for False")

parser.add_option("-o", "--orphansTF", action = "store", type= "string", dest= "orphanTF", help="1 to perform orphans classification, 0 for False")

(options, args) = parser.parse_args()


if options.Dir:
    Dir = options.Dir
else:
    parser.error("Please enter the path to Experiment LOG file -h")

if options.sigmafilename:
    sigmafilename = options.sigmafilename
else:
    parser.error("Enter the name of Sigma log file -h")

if options.MADfilename:
    MADfilename = options.MADfilename
else:
    parser.error("Enter the name of  MAD log file -h")

if options.orphanOUTname:
    orphanOUTname = options.orphanOUTname
else:
    parser.error("Enter the output file name for Orphan classificaition -h")

if options.Isotargetfile:
    Isotargetfile = options.Isotargetfile
else:
    parser.error("enter the name of target file for isotop correction -h")

if options.isoTF:
    isoTF = options.isoTF
else:
    parser.error("1 to perform istop correction, 0 for False -h")

if options.orphanTF:
    orphanTF = options.orphanTF
else:
    parser.error("1 to perform orphans classification, 0 for False -h")


def creatList(listoffiles):
    listofFile = []
    for line in open(listoffiles):
        if line != "\n":
            listofFile.append(line.strip())
    return listofFile

def orphanCalssificationANDisotopCorrection(experirmntlist, sigmaFiles, MADfiles, orphansOutName, isotopTargetFile, OrphanTF, IsotopTF):
    start_time = time.time()
    if __name__ == '__main__':
        listofFile = creatList(listoffiles=experirmntlist)

        if OrphanTF == "1" and IsotopTF == "1":

            func = partial(OrphanClassifier.ptmClassifier, sigmafIle = sigmaFiles, Outfilename =orphansOutName)
            with concurrent.futures.ProcessPoolExecutor() as executor:
                executor.map(func, listofFile)

            func1 = partial(Isotop_CorrectionAdvance.isotop, targetFile=isotopTargetFile, madFIle=MADfiles)
            with concurrent.futures.ProcessPoolExecutor() as executor:
                executor.map(func1, listofFile)

        elif OrphanTF == "1" and IsotopTF == "0":
            func = partial(OrphanClassifier.ptmClassifier, sigmafIle=sigmaFiles, Outfilename=orphansOutName)
            with concurrent.futures.ProcessPoolExecutor() as executor:
                executor.map(func, listofFile)

        elif OrphanTF == "0" and IsotopTF == "1":
            func1 = partial(Isotop_CorrectionAdvance.isotop, targetFile=isotopTargetFile, madFIle=MADfiles)
            with concurrent.futures.ProcessPoolExecutor() as executor:
                executor.map(func1, listofFile)


if __name__ == "__main__":
    orphanCalssificationANDisotopCorrection(experirmntlist=Dir, sigmaFiles=sigmafilename, MADfiles=MADfilename, orphansOutName=orphanOUTname,
                                            isotopTargetFile=Isotargetfile, OrphanTF=orphanTF, IsotopTF=isoTF)
    # orphanCalssificationANDisotopCorrection(experirmntlist=r"E:\Users\nbagwan\Desktop\SHIFTS_lastChange\heart\TestEXE1\Peakpicking_Log.txt",
    #                                     sigmaFiles=r"E:\Users\nbagwan\Desktop\SHIFTS_lastChange\heart\TestEXE1\sigmaCalculations.txt",
    #                                     MADfiles=r"E:\Users\nbagwan\Desktop\SHIFTS_lastChange\heart\TestEXE1\MAD_and_FWHM_calculations.txt",
    #                                     orphansOutName="OrphanWithOtherData.txt", isotopTargetFile="AllWithSequence-massTag.txt", OrphanTF="0", IsotopTF="1")
