__author__ = "Navratan Bagwan"
"version = 0.1, Date: 2019-02"

import os
import pdb
from optparse import OptionParser

usage = "usage:-P -B -X  -f -O -w -p -e -s\n Use -h for detailed help for parameter\nversion = 0.1, Date: 2019-02" \
        "\npeakFiltering is a part (module 2) of SHIFTS \n(Systematic, Hypothesis-free Identification of PTMs with controlled \nFDR " \
        "based on ultra-Tolerant database search)is a program made in the Jesus \nVazquez Cardiovascular Proteomics Lab at Centro Nacional de " \
        "Investigaciones \nCardiovasculares,for highthroughput PTM( Post translation modifications )processing.\n" \
        "peakFiltering, can be used to reprocess and filter slope and peak obtained from previous step"

parser = OptionParser(usage=usage)

#parser.add_option("")
parser.add_option("-i", "--histofile", action="store", type="string", dest="histofilename", help="Enter the path to histogramfile")

parser.add_option("-a", "--apexthreshold", action="store", type="string", dest="apexfreq", help="Enter the apex freq filter, else use 0")

parser.add_option("-g", "--guass", action="store", type="string", dest="guassFilter", help="Enter the firstGaussian filter, else use 0")

parser.add_option("-s", "--analysisname", action="store", type="string", dest="preFixname", help="Enter the name for outputfile")

(options, args) = parser.parse_args()

if options.histofilename:
    histofilename = options.histofilename
else:
    parser.error("Enter the path to histogramfile -h")

if options.apexfreq:
    apexfreq = options.apexfreq
else:
    parser.error("Enter the apex freq filter, else use 0, -h")

if options.guassFilter:
    guassFilter = options.guassFilter
else:
    parser.error("Enter the firstGaussian filter, else use 0, -h")

if options.preFixname:
    preFixname = options.preFixname
else:
    parser.error("Enter the name for outputfile -h")

def MassTagging2Seq_T(histogramFile, apex_thershold, firstGauss, prefixname):
    filepath = os.path.dirname(histogramFile)
    width = []
    mass = 0
    apexList = []
    #if target:
    filename = [filepath + str("/" + prefixname) + "Target_Peak_for_sequence_Tagging.txt", filepath + str("/" + prefixname) + "Target_ApexList_mass.txt"]
    # else:
    #     filename = [filepath + "/decoy_Peak_for_sequence_Tagging.txt", filepath + "/decoy_ApexList_mass.txt"]
    w2 = open(filename[0], "w")
    ApexList = open(filename[1], "w")
    Apex = []
    outputlist1 = [
        ["Bin", "\t", "Frequency", "\t", "Slope1", "\t", "Slope2", "\t", "peak-Width", "\t", "peak-Apex", "\t",
         "Intercept", "\t", "Tag",
         "\n"]]
    with open(histogramFile) as slopeFile:
        next(slopeFile)
        firstDerivative = []
        width = []
        mass = 0
        for line4 in slopeFile:
            if line4 != "\n":
                splits4 = line4.split("\t")
                if splits4[4].strip() == "1.0":  ### width coloumn
                    width.append(line4.strip())
                    firstDerivative.append(float(splits4[2]))
                    if splits4[5].strip() == "1.0":  ### apex coloumn
                        mass = splits4[6].strip()
                        freqeuncy = splits4[1].strip()
                        # apexList.append(mass)

                elif splits4[4].strip() == "0.0":
                    if len(width) != 0 and mass!=0:
                        firstelement = abs(float(firstDerivative[0]))
                        lastelement = abs(float(firstDerivative[-1]))
                        if firstelement >= float(firstGauss) and lastelement >= float(firstGauss) and float(freqeuncy) >= float(apex_thershold):
                            apexList.append(mass)
                            #outputlist1.append(" " + "\t" + " " + "\t"+ " " + "\t" + " " + "\t" + " " + "\t" + " " + "\t" + " "+ "\n")
                            outputlist1.append("\n")
                            for each in width:
                                # print each.strip(), "\t", mass
                                outputlist1.append([str(each), "\t", str(mass), "\n"])
                    width = []
                    firstDerivative = []
                    mass = 0
                    freqeuncy = 0

    for mas in outputlist1:
        w2.writelines(mas)
    w2.close()

    for apex in apexList:
        ApexList.writelines(apex + "\n")
    ApexList.close()



    return filename[1]

if __name__ == "__main__":
    MassTagging2Seq_T(histogramFile=histofilename, apex_thershold=apexfreq, firstGauss=guassFilter, prefixname=preFixname)
# MassTagging2Seq_T(histogramFile=r"E:\Users\nbagwan\Desktop\SHIFTS_lastChange\heart\TestEXE1\testhistogramFile.txt",
#                   apex_thershold=150, firstGauss = 300, prefixname = "test")