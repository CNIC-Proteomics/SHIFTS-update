__author__ = 'Navratan Bagwan'
import all_stats
import pdb
import os
import glob
import re
from collections import Counter
import numpy
def creatPeakApexfromList(peakapexfile):
    listofApex = []
    for line in open(peakapexfile):
        if line!= "\n":
            apex = line.split("\t")[0].strip()
            curatedApex = all_stats.check(float("%.6f" % float(apex)), 6)

            listofApex.append(str(curatedApex))

    return listofApex

# creatPeakApexfromList(peakapexfile=r"D:\CNIC\SHIFTS0.2\data\heart\OUTPUT-test\testApexTarget_ApexList_mass.txt")

def CountDic_and_WorkingDic(fileName, folder, listofApexs, column):
    Tdic = {}
    apexCountlist = []
    filetouse = folder + "/" + fileName

    with open(filetouse) as filename:
        next(filename)
        for line in filename:

            splits = line.split("\t")
            # label = splits[18].strip()
            tag = splits[23].strip()
            deltapeptide = splits[int(column)-1].strip()
            plainSeq = re.sub("[^a-zA-Z]+", "", deltapeptide)
            apexinDeltapeptide = "".join(re.findall("[^a-zA-Z]+", deltapeptide))
            apexinDeltapeptideDM = str(apexinDeltapeptide.split("[")[1].split("]")[0])
            # if apexinDeltapeptideDM == "0.000043":
            #     print(apexinDeltapeptideDM)
            #     pdb.set_trace()
            if apexinDeltapeptideDM in listofApexs and tag != "Orphan":
                apexCountlist.append(apexinDeltapeptideDM)
                if plainSeq not in Tdic:
                    Tdic[plainSeq] = [[deltapeptide, apexinDeltapeptide, apexinDeltapeptideDM, line.strip()]]
                else:
                    Tdic[plainSeq].append([deltapeptide, apexinDeltapeptide, apexinDeltapeptideDM, line.strip()])


    countDic = Counter(apexCountlist)
    # print(len(apexCountlist))

    return countDic, Tdic



    #         if label == "Target" and tag!= "Orphan":
    #             Tdic[line.strip()] = 1
    #         else:
    #             Ddic[line.strip()] = 1
    #
    # return Tdic, Ddic


def isotop(folder, targetFile, madFIle, apexList, col):

    folder1 = os.path.dirname(folder)
    writefile = folder1 + "/IsotopCorrection_TargetData_withSequence-massTag.txt"
    headerList = all_stats.getHeader(folder1 + "/" + targetFile)
    header = "\t".join(headerList) + "\t" + "Corr_Seq-mass" + "\t" + "Corr_mass" + "\t" +"Monoisotop_T/F" + "\n"

    apexilstfromInput = creatPeakApexfromList(peakapexfile=apexList)


    wfile = open(writefile, "w")
    wfile.write(header)

    seconslastElement = all_stats.getIndex(folder1 + "/" + targetFile)

    print (madFIle)
    madfilePath = glob.glob(os.path.join(folder1, str(madFIle)))
    print (madfilePath)
    FWHM = float(all_stats.fullWidthHalfMaximum(MADfile=madfilePath[0]))
    ApexFequency, targetFileDic = CountDic_and_WorkingDic(targetFile, folder=folder1, listofApexs=apexilstfromInput, column=col)


    c13Mass = 1.003354
    c13end = c13Mass + float(FWHM)
    c13start = c13Mass - float(FWHM)

    for everySeq in targetFileDic:
        AllMassList = targetFileDic[everySeq]
        AllMassList.sort(key=lambda x: x[2])

        massList = [item[2] for item in AllMassList]

        uniqueMasses = numpy.unique(massList)
        DicforTrueIsotop = {}
        c = 0
        for firstE in range(len(uniqueMasses)-1):
            for SecondE in range(firstE+1, len(uniqueMasses)):
                diff = float(uniqueMasses[SecondE]) - float(uniqueMasses[firstE])

                if c13end >= abs(float(diff)) >= c13start:
                    # tempList = [(ApexFequency[str(uniqueMasses[firstE])]) ,int(ApexFequency[str(uniqueMasses[SecondE])])]

                    if int(ApexFequency[str(uniqueMasses[firstE])]) < int(ApexFequency[str(uniqueMasses[SecondE])]):
                        to_beReplaced = str(uniqueMasses[firstE])
                        replace_by = str(uniqueMasses[SecondE])
                        c = c + 1

                        DicforTrueIsotop[to_beReplaced] = replace_by
                    elif int(ApexFequency[str(uniqueMasses[firstE])]) > int(ApexFequency[str(uniqueMasses[SecondE])]):

                        to_beReplaced = str(uniqueMasses[SecondE])
                        replace_by = str(uniqueMasses[firstE])

                        c = c + 1

                        DicforTrueIsotop[to_beReplaced] = replace_by
                        break

                # elif c13end >= abs(float(diff)) >= c13start and int(ApexFequency[str(uniqueMasses[firstE])]) > int(ApexFequency[str(uniqueMasses[SecondE])]):
                #     c = c + 1
                #     #DicforTrueIsotop[str(uniqueMasses[SecondE])] = str(uniqueMasses[firstE])
                #     DicforTrueIsotop[str(uniqueMasses[firstE])] = str(uniqueMasses[SecondE])
                #     break



        # if c == 0:
        #     print(everySeq)
        #     pdb.set_trace()

        if c > 0:
            for test in AllMassList:

                if str(test[2]) not in DicforTrueIsotop:
                    wfile.write(test[3]+ "\t" + str(test[0]) + "\t" + str(test[2]) + "\t" + "no"+ "\n")

                else:

                    newmass = DicforTrueIsotop[str(test[2])]
                    deltaP = test[0]
                    if "[" in deltaP:
                        for index in range(0, len(deltaP)):
                            if deltaP[index] == "[":
                                startofMOD = index
                            if deltaP[index] == "]":
                                endofMOD = index
                        tempPepe = deltaP[startofMOD + 1:endofMOD]
                        # newclose = all_stats.check(list2[0], 2)
                        finalSeqMassTag = deltaP.replace(tempPepe, str(newmass))

                        wfile.write(test[3] + "\t" + str(finalSeqMassTag) + "\t" + str(newmass) + "\t" + "yes"+ "\n")

        else:
            for test in AllMassList:
                wfile.write(test[3]+ "\t" + str(test[0]) + "\t" + str(test[2]) + "\t" + "no"+ "\n")

    wfile.close()







