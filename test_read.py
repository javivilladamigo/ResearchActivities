txt_file = "fileList_v9/GluGluToH_HToJPsiG_JPsiToMuMu_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL18.txt"
fileList = []
with open(txt_file) as file:
    for line in file:
        fileList.append("root://cms-xrd-global.cern.ch/" + line.rstrip("\n"))