from os import listdir
from os.path import isfile, join

def test_read_remote(txt_file = "fileList_v9/GluGluToH_HToJPsiG_JPsiToMuMu_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL18.txt"):
    fileList = []
    with open(txt_file) as file:
        for line in file:
            fileList.append("root://cms-xrd-global.cern.ch/" + line.rstrip("\n"))
    return fileList

def test_read_local(dir_file = "/Users/javi/Documents/Padova/Research.nosync/ResearchActivities/data.nosync/OutputTest/"):
    fileList = [f for f in listdir(dir_file) if (isfile(join(dir_file, f)) and ".root" in f)]
    return fileList

