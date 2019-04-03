__author__ = "Navratan Bagwan"
"version = 0.1, Date: 2019-04"
#navratan.bagwan@sund.ku.dk
#change from the last version: this version uses average frequency based on the user input##

import glob
import os
import concurrent.futures
from functools import partial
import itertools
import time
import shutil
from optparse import OptionParser
import multiprocessing
import Comet_PTM_processingScript
import all_stats
import SigmaAndFDR_process
import SigmaCalculations
import peakPicking
import matplotlib
matplotlib.use('agg')

usage = "usage:-P -B -X  -f -O -w -p -e -s\n Use -h for detailed help for parameter\nversion = 0.1, Date: 2019-02" \
        "\nvertex.exe is a part (module 1) of SHIFTS \n(Systematic, Hypothesis-free Identification of PTMs with controlled \nFDR " \
        "based on ultra-Tolerant database search)is a program made in the Jesus \nVazquez Cardiovascular Proteomics Lab at Centro Nacional de " \
        "Investigaciones \nCardiovasculares,for highthroughput PTM( Post translation modifications )processing.\n" \
        "vertex.exe, can be used to preprocess and calibrate the Comet output\nIn addition you can also choose to do peak picking modelling with or without first part "

parser = OptionParser(usage=usage)

parser.add_option("-P", "--Path2master", action="store", type="string", dest="Dir", help="Enter the path to masterFile")

parser.add_option("-B", "--BinSize", action="store", type="string", dest="binSize", help="Enter the bin size for peak modelling")

parser.add_option("-X", "--Xcor", action="store", type="string", dest="CorrXcor", help="Enter the Xcor thershold for Calibration")

parser.add_option("-f", "--FDRthreshold", action="store", type="string", dest="FDRs", help="Enter the global FDR threshold; 0.01 for 1% FDR")

parser.add_option("-O", "--OutPutname", action="store", type="string", dest="Outname", help="Enter the output folder name")

parser.add_option("-w", "--WindowRange", action="store", type="string", dest="window", help="Enter the sliding window for peak picking; Ex 7 or 9")

parser.add_option("-p", "--PeakReRun", action="store", type="string", dest="rerunPeak", help="chose to rerun peakpicking independently; 1 for T, 0 for F")

parser.add_option("-e", "--experimentLog", action="store", type="string", dest="ExpLog", help="Enter the path to peakpickingLog file if rerun (-p=1)")

parser.add_option("-s", "--suffix", action="store", type="string", dest="suffixname", help="Enter the suffix name for output file if rerun (-p=1)")

parser.add_option("-n", "--normalizeFreq", action="store", type="string", dest="norm", help="choose if you want to normalize frequency; 1 for T, 0 for F")

(options, args) = parser.parse_args()


if options.Dir:
    Dir = options.Dir
else:
    parser.error("Please enter the path to mater file -h")

if options.binSize:
    binSize = options.binSize
else:
    parser.error("Please enter the Bin size you want use -h")

if options.CorrXcor:
    CorrXcor = options.CorrXcor
else:
    parser.error("Please enter Xcor thershold you want use for Calibration -h")

if options.FDRs:
    FDRs = options.FDRs
else:
    parser.error("please enter global-FDR you want to use; 0.01 for 1% FDR")

if options.Outname:
    Outname = options.Outname
else:
    parser.error("Please enter name for output folder -h")

if options.window:
    window = options.window
else:
    parser.error("Please the window size for peakpicking model -h")

if options.rerunPeak:
    rerunPeak = options.rerunPeak
else:
    parser.error("peakRerun true = 1 false = 0 -h")

if options.ExpLog:
    ExpLog = options.ExpLog
else:
    parser.error("path to peakpicking log file, if -p is 1, -h")

if options.suffixname:
    suffixname = options.suffixname
else:
    parser.error("Please enter the suffix name for rerun output file name -h")

if options.norm:
    norm = options.norm
else:
    parser.error("choose if you want to normalize frequency; 1 for T, 0 for F, -h")


def calibrate_every_path(pathTomasterFIle, bins, Xcorthershold, outputfolder, fdrFilter, window, PeakRerun,expLog,rexentension, smoothTF ):
    LogString = "Location of master file: " + str(pathTomasterFIle) + "\n" + "Bins used for peak picking model: " + str(Xcorthershold)\
                + "\n" + "Output folder name: " + str(outputfolder) + "\n" + "Global FDR: " + str(fdrFilter) + "\n" + \
                "Sliding window for guassian modelling (peak picking): " + str(window) + "\n" + "Peak picking part rerun TRUR/FALSE:  0 denotes False and 1 True): " \
                + str(PeakRerun) + "\n" + "Experiment Log file (if you have 1 in Rerun, you need to use the log file made from calibration step : )" + str(expLog) + "\n"\
                + "Suffix name for outfile if you Rerun peapicking: " + str(rexentension)

    if PeakRerun == "0":
        ##### the first method "Comet_PTM_processingScript.ProcessingFile" takes all the raw comet output from master file and calculates a machine error and creates an output folder###
        #### for more details check the method Comet_PTM_processingScript.py#####
        start_time = time.time()
        if __name__ == '__main__':
            multiprocessing.freeze_support()

            filelist = []
            master_list, child_list = all_stats.master_table(masterFile=pathTomasterFIle)
            func = partial(Comet_PTM_processingScript.ProcessingFile, outfoldername=outputfolder, scoreInput=Xcorthershold)
            with concurrent.futures.ProcessPoolExecutor() as executor:
                filelist += executor.map(func, child_list)

            listofDM =[]
            listofDic = {}
            listofFiles = filelist
            func1 = partial(SigmaAndFDR_process.filemakerforSigmaFinder, scoreInput= Xcorthershold,fdrthershol=fdrFilter)
            with concurrent.futures.ProcessPoolExecutor() as executor:
                listofDM += executor.map(func1, filelist)
            listofDM_ppm = list(itertools.chain.from_iterable(listofDM))
            SigmaCalculations.SigmaCalculation(deltaPPMlist=listofDM_ppm, processingFileList=filelist)


            funcMAD = partial(all_stats.madCalculation, scoreInput = Xcorthershold)
            with concurrent.futures.ProcessPoolExecutor() as executor:
                executor.map(funcMAD, filelist)

            #histogramFilenameList = peakPicking.SlopeCalculation(filelist=filelist, bintoHist=bins, windowSize=window)
            peakPickingLogfile = "Peakpicking_Log.txt"
            w = open(peakPickingLogfile, "w")
            for logfile in filelist:
                w.write(logfile + "\n")
            w.close()
            vertexLogfile = "vertextParameter_Log.txt"

            w1 = open(vertexLogfile, "w")
            for i in LogString:
                w1.writelines(i)
            w1.close()

            for eachfile in filelist:
                foldername = os.path.dirname(eachfile)
                shutil.copy(peakPickingLogfile, foldername)
                shutil.copy(vertexLogfile, foldername)
            os.remove(peakPickingLogfile)
            os.remove(vertexLogfile)
            print ("starting the peakPIcking model")
            histogramFilenameList = peakPicking.callSlopeMethod(filelist=filelist, bintoHist=bins, windowSize=int(window), rerun=PeakRerun,extensionName=rexentension, smotheeingTF=int(smoothTF))
            print("---%s seconds ---" % (time.time() - start_time))
    else:
        if __name__ == '__main__':
            multiprocessing.freeze_support()
            vertexLogfile = str(rexentension)+"_"+"vertextParameter_Log.txt"
            w1 = open(vertexLogfile, "w")
            for i in LogString:
                w1.writelines(i)
            w1.close()

            for eachfile in open(expLog):
                foldername = os.path.dirname(eachfile)
                shutil.copy(vertexLogfile, foldername)
            os.remove(vertexLogfile)

            print ("starting the peakPicking model with different parameter than previous")
            histogramFilenameList = peakPicking.callSlopeMethod(filelist=expLog, bintoHist=bins, windowSize=int(window), rerun=PeakRerun,extensionName=rexentension, smotheeingTF=int(smoothTF))

if __name__ == "__main__":
    calibrate_every_path(pathTomasterFIle=Dir, bins=binSize,Xcorthershold=CorrXcor, outputfolder=Outname, fdrFilter=FDRs,
                         window=window,PeakRerun=rerunPeak,expLog=ExpLog, rexentension=suffixname, smoothTF=norm)

# calibrate_every_path(pathTomasterFIle=r"E:\Users\nbagwan\Desktop\SHIFTS_lastChange\MasterFile.txt", Xcorthershold=0.22,
#                      outputfolder = "ModulesTest", bins = 0.001, fdrFilter=0.05, window= 7, PeakRerun="0",
#                      expLog=r"E:\Users\nbagwan\Desktop\SHIFTS_lastChange\liver\Module1OUT\Peakpicking_Log.txt", rexentension="test")

