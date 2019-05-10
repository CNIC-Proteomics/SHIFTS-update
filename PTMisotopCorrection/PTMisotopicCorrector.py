__author__ = "Navratan Bagwan"
"version = 0.1, Date: 2019-02"

import pdb
import all_stats
import time

import IsotopCorrectionDM
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

parser.add_option("-m", "--madfile", action="store", type="string", dest="MADfilename", help="Enter the name of MAD log file")

parser.add_option("-I", "--targetOfiso", action="store", type="string", dest="Isotargetfile", help="enter the name of target file for isotop correction")

parser.add_option("-a", "--Apexfile", action="store", type="string", dest="ListDM", help="enter the name of file, containing peak apexes")

parser.add_option("-c", "--DeltpeptideColnumber", action="store", type="string", dest="DpepCol", help="enter the coloumn number, which has deltapeptide")

(options, args) = parser.parse_args()


if options.Dir:
    Dir = options.Dir
else:
    parser.error("Please enter the path to Experiment LOG file -h")


if options.MADfilename:
    MADfilename = options.MADfilename
else:
    parser.error("Enter the name of  MAD log file -h")

if options.Isotargetfile:
    Isotargetfile = options.Isotargetfile
else:
    parser.error("enter the name of target file for isotop correction -h")

if options.ListDM:
    ListDM = options.ListDM
else:
    parser.error("enter the name of file, containing peak apexes")

if options.DpepCol:
    DpepCol = options.DpepCol
else:
    parser.error("enter the coloumn number, which has deltapeptide")


def creatList(listoffiles):
    listofFile = []
    for line in open(listoffiles):
        if line != "\n":
            listofFile.append(line.strip())
    return listofFile

def orphanCalssificationANDisotopCorrection(experirmntlist, MADfiles, isotopTargetFile, peakFIle, colnumber):
    start_time = time.time()
    if __name__ == '__main__':
        listofFile = creatList(listoffiles=experirmntlist)

        # for i in listofFile:
        #     IsotopCorrectionDM.isotop(folder=i, targetFile=isotopTargetFile, madFIle=MADfiles,apexList=r"D:\CNIC\SHIFTS0.2\data\testTissue\bin002_w7_SL18000Target_ApexList_mass.txt")
        func1 = partial(IsotopCorrectionDM.isotop, targetFile=isotopTargetFile, madFIle=MADfiles, apexList=peakFIle, col= colnumber)
        with concurrent.futures.ProcessPoolExecutor() as executor:
            executor.map(func1, listofFile)


if __name__ == "__main__":
    orphanCalssificationANDisotopCorrection(experirmntlist=Dir, MADfiles=MADfilename,isotopTargetFile=Isotargetfile, peakFIle=ListDM, colnumber=DpepCol)
    # orphanCalssificationANDisotopCorrection(experirmntlist=r"E:\Users\nbagwan\Desktop\SHIFTS_lastChange\heart\TestEXE1\Peakpicking_Log.txt",
    #                                     sigmaFiles=r"E:\Users\nbagwan\Desktop\SHIFTS_lastChange\heart\TestEXE1\sigmaCalculations.txt",
    #                                     MADfiles=r"E:\Users\nbagwan\Desktop\SHIFTS_lastChange\heart\TestEXE1\MAD_and_FWHM_calculations.txt",
    #                                     orphansOutName="OrphanWithOtherData.txt", isotopTargetFile="AllWithSequence-massTag.txt", OrphanTF="0", IsotopTF="1")
