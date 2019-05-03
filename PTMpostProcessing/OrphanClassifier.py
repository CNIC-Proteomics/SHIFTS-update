__author__ = "Navratan Bagwan"

import all_stats
import pdb
import numpy
import os
import glob


def ptmClassifier(listofpaths, sigmafIle, Outfilename):
    print (sigmafIle)
    #sigmaValue = float(all_stats.sigmaFinder(sigmafIle))
    outpath = os.path.dirname(listofpaths)
    assignationFile = glob.glob(os.path.join(outpath, "AllWithSequence-massTag.txt"))
    sigmafilePath = glob.glob(os.path.join(outpath, str(sigmafIle)))
    print (sigmafilePath)
    sigmaValue = float(all_stats.sigmaFinder(sigmafilePath[0]))
    tag = outpath + "/" + Outfilename
    tagFile = open(tag, "a")
    header = ["Scan", "\t", "SearchEngineRank", "\t", "Charge", "\t", "Exp Mh", "\t", "Theo mh+", "\t", "Exp MZ", "\t", "Xcor", "\t",
         "Seq", "\t", "RetentionTIme", "\t", "Protein", "\t", "Delta_MOD" ,"\t", "B_series","\t", "Y_series", "\t", "Jumps", "\t",
         "DeltaPeptide", "\t", "FileName", "\t", "CorXcor" "\t", "New_expMH", "\t", "label", "\t",
         "median", "\t", "Calibrated Delta MH", "\t", "Calibrated Delta MZ","\t", "Calibrated EXP MZ","\t", "Tags", "\t",
                     "fastaDescription", "\t", "seq_mass", "\t", "PeakApex/Dmass", "\t", "DeltaPeptide", "\t", "PPMerror Orphancalss",
              "\t", "ClassNumber", "\t", "Apex/OrphanClassMass", "\t", "FinalSeq", "\n"]

    tagFile.writelines("".join(header))

    mainList = []
    finalList = []
    nonOrphanList = []
    # for everypath in listofpaths:
    #     filename = everypath + "/" + "NotassignedSequences.txt"

    with open(assignationFile[0]) as unassignedFile:
        next(unassignedFile)
        for line in unassignedFile:
            if line != "\n":
                splits = line.split("\t")
                calibratedDelta_MH = splits[20].strip()
                if splits[23].strip() == "Orphan":
                    mainList.append([outpath + "\t" + line.strip() + "\t", float(calibratedDelta_MH.strip())])
                else:
                    nonOrphanList.append([outpath + "\t" + line.strip()])

    # print len(mainList)
    # print len(nonOrphanList)
    for notOrphan in nonOrphanList:

        lineToadd = "\t".join(notOrphan[0].split("\t")[1:]) + "\t" + "NA" + "\t" + "NA" + "\t" + str(notOrphan[0].split("\t")[27]) + "\t" + str(notOrphan[0].split("\t")[28]) + "\n"
        tagFile.writelines(lineToadd)


    mainList.sort(key=lambda x: x[1])
    for index in range(len(mainList)):
        m1 = mainList[index][1]
        if index == 0:
            m0 = 0
        else:
            m0 = mainList[index - 1][1]
        try:
            ppmError = abs((m0 - m1) / m0) * 1000000
        except ZeroDivisionError:
            ppmError = 0.0

        finalList.append([mainList[index][0], str(ppmError)])

    classList = []
    for index1 in range(len(finalList)):
        if index1 == 0:
            classificationTerm = 1
            classList.append([finalList[index1][0] + finalList[index1][1], str(classificationTerm)])
        else:
            if float(finalList[index1][1]) <= sigmaValue:
                classificationTerm = classificationTerm
                classList.append([finalList[index1][0] + finalList[index1][1], str(classificationTerm)])
            else:
                classificationTerm = classificationTerm + 1
                classList.append([finalList[index1][0] + finalList[index1][1], str(classificationTerm)])

    classDic = {}
    for classes in classList:
        if classes[1].strip() not in classDic:
            classDic[classes[1].strip()] = [float(classes[0].split("\t")[21])]
        else:
            classDic[classes[1].strip()].append(float(classes[0].split("\t")[21]))

    for dikey in classDic:
        medianValuelist = classDic[dikey]
        if len(medianValuelist) < 2:
            massOfOrphan = medianValuelist[0]
            # intformula = int(massOfOrphan / float(binvalue)) * float(binvalue) + float(binvalue) / 2

            orphanPTM = all_stats.check(massOfOrphan, 6)
            classDic[dikey].append(orphanPTM)
        else:
            massOfOrphan = numpy.median(medianValuelist)
            # intformula = int(massOfOrphan / float(binvalue)) * float(binvalue) + float(binvalue) / 2
            orphanPTM = all_stats.check(massOfOrphan, 6)
            classDic[dikey].append(orphanPTM)

    for appendingPTMinList in classList:
        if appendingPTMinList[1] in classDic:
            ptmINT = classDic[appendingPTMinList[1]][1:][0]
            appendingPTMinList.append(all_stats.check(float(ptmINT),6))
            appendingPTMinList.append("Orphan")

    DicTomergeFiles = {}
    for ii in classList:

        PathTO_ptmFile = ii[0].split("\t")[0]
        if PathTO_ptmFile not in DicTomergeFiles:
            DicTomergeFiles[PathTO_ptmFile] = [ii]
        else:
            DicTomergeFiles[PathTO_ptmFile].append(ii)



    for pathInDic in DicTomergeFiles:

        if pathInDic in outpath:

            # tag = pathInDic + "/" + "ClassTest_AllWithSequence-massTag.txt"
            # tagFile = open(tag, "a")
            for linestring in DicTomergeFiles[pathInDic]:

                deltaPeptide = linestring[0].split("\t")[15]
                if deltaPeptide.isalpha():
                    finalSeqMassTag = deltaPeptide
                else:
                    for index in range(0, len(deltaPeptide)):
                        if deltaPeptide[index] == "[":
                            startofMOD = index
                        if deltaPeptide[index] == "]":
                            endofMOD = index
                    tempPepe = deltaPeptide[startofMOD + 1:endofMOD]
                    newclose = all_stats.check(float(linestring[2]), 6)
                    finalSeqMassTag = deltaPeptide.replace(tempPepe, str("{:.6f}".format(float(newclose))))


                lineToadd = "\t".join(linestring[0].split("\t")[1:]) + "\t" + str(linestring[1]) + "\t" + str(linestring[2]) + "\t" + finalSeqMassTag + "\n"
                tagFile.writelines(lineToadd)

    tagFile.close()

# def callCallisifer(listoffiles, sigmafiles):
#     pdb.set_trace()
#     calssifierFucntion = partial(ptmClassifier, sigmafile=sigmafiles)
#     with concurrent.futures.ProcessPoolExecutor() as executor:
#         executor.map(calssifierFucntion, listoffiles)
