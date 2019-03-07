__author__ = "Navratan Bagwan"

import os
import shutil
import all_stats

separations = "/"


#dir_list(dir_name="E:/Users/nbagwan/Desktop/OpenSearchPaper/Improved_isPTM", subdir= True, args="[txt]", fn="Carot_ND_Fr1_PSMs.txt")
###### the method takes three input, scoreInput for CorrXcor threshold, file produced from Comet_PTM_processingScript.py, and FDr treshold provided by user as input #####

def filemakerforSigmaFinder(inputFile,scoreInput,fdrthershol):
    deltaPPMdic = {}
    deltaPPMlist = []
    outfile = []
    globalFDRdic = {}
    #for everyProcssedFile in processingFileList:
    with open(inputFile) as f:
        next(f)
        for line in f:
            line = line.strip("\n")
            if line != "\n":
                splits = line.strip().split("\t")
                experimental_MZ = float(splits[22].strip("\" "))
                theoritical_MZ = float(splits[4].strip("\" ")) / float(splits[2].strip("\" "))
                rawFilename = str(splits[15]).strip("\" ")
                Xcorthershold = float(scoreInput)
                newXcor = float(splits[16].strip("\" "))
                delta_mz = float(splits[21].strip("\" "))
                labels = str(splits[18]).strip("\" ")
                calibrated_deltaMH = float(splits[20].strip("\" "))

                KeyID = os.path.dirname(inputFile)
                if newXcor >= Xcorthershold and delta_mz > -0.1 and delta_mz < 0.1:  ##### new charge condition. search engine rank fixed
                    delta_ppm = (experimental_MZ - theoritical_MZ) / theoritical_MZ * 1000000
                    if delta_ppm > -20 and delta_ppm < 20:
                        deltaPPMlist.append(delta_ppm)
                        if rawFilename not in deltaPPMdic:
                            deltaPPMdic[rawFilename] = [delta_ppm]
                        else:
                            deltaPPMdic[rawFilename].append(delta_ppm)

                        ################################## making a dictionary for global FDR  ########################
                if KeyID not in globalFDRdic:
                    if calibrated_deltaMH >= -56.0 and calibrated_deltaMH <=500: #### we used -56 as start for global
                        globalFDRdic[KeyID] = [[newXcor, labels]]
                else:
                    if calibrated_deltaMH >= -56.0 and calibrated_deltaMH <=500:
                        globalFDRdic[KeyID].append([newXcor, labels])

    alpha = 1
    for mod in globalFDRdic:
        decoy, target = 0.0, 0.0
        globalFDRdic[mod].sort(key=lambda row: (row[0]), reverse=True)
        for eachline in globalFDRdic[mod]:
            if eachline[1] == "Decoy":
                decoy += 1.0
                eachline.append(decoy * float(alpha))
                eachline.append(target)
                try:
                    localFDR = float(decoy * float(alpha) / target)
                except ZeroDivisionError:
                    localFDR = "1"
                eachline.append(str(localFDR).strip())

            elif eachline[1] == "Target":
                target += 1.0
                eachline.append(decoy * float(alpha))
                eachline.append(target)
                try:
                    localFDR = float(decoy * float(alpha) / target)
                except ZeroDivisionError:
                    localFDR = "1"
                eachline.append(str(localFDR))

    ##### method writes a file (FDR log file) mentioning the Global FDR trreshold in CorrXcor based on user input #####
    for thershold in globalFDRdic:
        outfile1 = []
        fileName1 = thershold + "/globalFDRforSmallDB.txt"
        for each in globalFDRdic[thershold]:
            if float(each[4]) < float(fdrthershol):
                globalFDR_xcor = float(each[0])
                numberOFtarget = each[3]
                numberOFdecoy = each[2]
            else:
                break

        outfile1.append("globalFDR: " + str(globalFDR_xcor) + "\n" + "number of Target: " + str(
            numberOFtarget) + "\n" + "number of decoys: " + str(numberOFdecoy))
        w1 = open(fileName1, "w")
        for ii in outfile1:
            w1.writelines(ii)
        w1.close()

    return deltaPPMlist
