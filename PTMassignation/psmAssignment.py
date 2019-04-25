__author__ = "Navratan Bagwan"
import pdb
import glob
import os
import all_stats
def createFasta(fasta):
    proteinFast = {}
    for line6 in open(fasta):
        if line6 != "\n":
            if line6.startswith(">"):
                if "DECOY_" not in line6 and "INV_" not in line6 and "XXX_" not in line6 and "_INV_" not in line6:
                    splits6 = line6.split("|")
                    if splits6[1].strip() not in proteinFast:
                        proteinFast[splits6[1].strip()] = line6.strip()

    return proteinFast

def modfiyPeptide(plainSeq, deltapaptide, peakapex):
    for index in range(0, len(deltapaptide)):
        if deltapaptide[index] == "[":
            startofMOD = index
        if deltapaptide[index] == "]":
            endofMOD = index
    tempPepe = deltapaptide[startofMOD + 1:endofMOD]

    newclose = all_stats.check(float(peakapex), 6)
    finalSeqMassTag = deltapaptide.replace(tempPepe, str(newclose))
    seq_mass = plainSeq + "_" + str(newclose)
    return seq_mass, finalSeqMassTag

def taggingMethodForCOMET_PTM(calibrationFile, fasta, localFDRfilter, peakFDRfilter):
    outpath = os.path.dirname(calibrationFile)
    NonOrphanList = ["Scan", "\t", "SearchEngineRank", "\t", "Charge", "\t", "Exp Mh", "\t", "Theo mh+", "\t", "Exp MZ", "\t", "Xcor", "\t",
         "Seq", "\t", "RetentionTIme", "\t", "Protein", "\t", "Delta_MOD" ,"\t", "B_series","\t", "Y_series", "\t", "Jumps", "\t",
         "DeltaPeptide", "\t", "FileName", "\t", "CorXcor" "\t", "New_expMH", "\t", "label", "\t",
         "median", "\t", "Calibrated Delta MH", "\t", "Calibrated Delta MZ","\t", "Calibrated EXP MZ","\t", "Tags", "\t",
                     "fastaDescription", "\t", "seq_mass", "\t", "PeakApex/Dmass", "\t", "DeltaPeptide", "\n"]


    OrphanList = []
    outfileName = outpath + "/AllWithSequence-massTag.txt"
    w = open(outfileName, "w")
    outlist = []
    globalFDRFile = glob.glob(os.path.join(outpath, "globalFDRforSmallDB.txt"))

    GlobalXcorrThreshold = all_stats.globalfdr(globalFDRFile[0])
    fastaHeaderDic = createFasta(fasta=fasta)

    DictofallData = outpath + "/Peak_and_Slope_FDRfile-Smallchunk.txt"

    for line in open(DictofallData):#[calibrationFile]:
        splits = line.split("\t")
        proteinComet = splits[9].strip()
        Label = splits[18].strip()
        corrXcorr = splits[16].strip()
        localFDRinput = splits[27].strip()
        PeakFDRinput = splits[26].strip()
        peakApex = splits[23].strip()
        deltaPepetide = splits[14].strip()
        plainseq = splits[7].strip()
        calDeltaMH = splits[20].strip()
        if Label == "Target" and peakApex != "Orphan":
            if "|" in proteinComet:
                proteinAccestion = proteinComet.split("|")[1]
            else:
                proteinAccestion = proteinComet.split("|")[0]

            fastaDescription = fastaHeaderDic[proteinAccestion]

            if float(localFDRinput) < float(localFDRfilter): ##### fall in local fdr
                if float(corrXcorr) > float(GlobalXcorrThreshold): ### should also fall in global fdr
                    plainseq_mass, NewdeltaPpeptide = modfiyPeptide(plainSeq = plainseq, deltapaptide = deltaPepetide, peakapex = float(peakApex))
                    NonOrphanList.append(["\t".join(splits[:23]), "\t", "NA", "\t", str(fastaDescription), "\t", str(plainseq_mass),"\t", str(peakApex), "\t", str(NewdeltaPpeptide), "\n"])


            elif float(localFDRinput) > float(localFDRfilter) and float(PeakFDRinput) < float(peakFDRfilter):
                if float(corrXcorr) > float(GlobalXcorrThreshold):
                    plainseq_mass, NewdeltaPpeptide = modfiyPeptide(plainSeq = plainseq, deltapaptide = deltaPepetide, peakapex = float(peakApex))
                    NonOrphanList.append(["\t".join(splits[:23]), "\t", "NA", "\t", str(fastaDescription), "\t", str(plainseq_mass), "\t", str(peakApex), "\t", str(NewdeltaPpeptide), "\n"])

        if Label == "Target" and peakApex == "Orphan":

            if "|" in proteinComet:
                proteinAccestion = proteinComet.split("|")[1]
            else:
                proteinAccestion = proteinComet.split("|")[0]

            fastaDescription = fastaHeaderDic[proteinAccestion]
            if float(localFDRinput) < float(localFDRfilter): ##### fall in local fdr
                if float(corrXcorr) > float(GlobalXcorrThreshold): ### should also fall in global fdr
                    plainseq_mass, NewdeltaPpeptide = modfiyPeptide(plainSeq = plainseq, deltapaptide = deltaPepetide, peakapex = float(calDeltaMH))
                    OrphanList.append(["\t".join(splits[:23]), "\t", "Orphan", "\t", str(fastaDescription), "\t", str(plainseq_mass), "\t", str(calDeltaMH), "\t", str(NewdeltaPpeptide), "\n"])

            elif float(localFDRinput) > float(localFDRfilter) and float(PeakFDRinput) < float(peakFDRfilter):
                if float(corrXcorr) > float(GlobalXcorrThreshold):
                    plainseq_mass, NewdeltaPpeptide = modfiyPeptide(plainSeq = plainseq, deltapaptide = deltaPepetide, peakapex = float(calDeltaMH))
                    OrphanList.append(["\t".join(splits[:23]), "\t", "Orphan", "\t", str(fastaDescription), "\t", str(plainseq_mass), "\t", str(calDeltaMH), "\t", str(NewdeltaPpeptide), "\n"])




    for nonorphan in NonOrphanList:
        w.writelines(nonorphan)

    for orphans in OrphanList:
        w.writelines(orphans)
    w.close()