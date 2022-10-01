#!/bin/env python

import os, sys
from global_paths import *
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from samesign import *
from hzz4l import *
from higgs import *
from triggerJPsi import *
from triggerMuPhoton import *
from triggerPhoton import *
from triggerMu import *

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)

parser.add_option("-l", "--list", action="store", type="string", dest="list", default="")
parser.add_option("-i", "--input", action="store", type="string", dest="input", default="")
parser.add_option("-o", "--output", action="store", type="string", dest="output", default="OutputTest/")
parser.add_option("-f", "--filter", action="store", type="string", dest="filter", default=None)
parser.add_option("-y", '--year', action='store', type=int, dest="year", default=0)
parser.add_option("-n", '--max', action='store', type=long, dest="maxEntries", default=1e9)
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False)
(options, args) = parser.parse_args()



fileList = []
if len(options.input) > 0:
    for f in options.input.split(','):
        if len(f) > 0: fileList += [REDIRECTOR + f]
if len(options.list) > 0:
    with open(options.list, "r") as f:
        for f in f.read().splitlines():
            if len(f) > 0: fileList += [REDIRECTOR + f]


#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAODv5/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7_ext2-v1/120000/BD54F09C-7429-9547-A169-467F4AFB9606.root"] #["root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAOD/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/40000/2CE738F9-C212-E811-BD0E-EC0D9A8222CE.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAODv5/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7_ext1-v1/260000/B2D94295-780F-AA4C-9A98-FB682EC1F910.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18NanoAODv5/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano1June2019_102X_upgrade2018_realistic_v19_ext3-v1/20000/8CC8B22B-6E75-F648-97DA-0ACB3AB5DB98.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAODv5/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1/100000/FAB838BB-6F66-6B4A-B772-A1B511F28ECB.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAODv5/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1/70000/9AB2845B-462E-8E42-8605-84606FB59047.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAODv5/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1/70000/A03EC876-43DF-8841-BD98-5E345C766D5C.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18NanoAODv5/WpWpJJ_EWK-QCD_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/Nano1June2019_102X_upgrade2018_realistic_v19-v1/40000/D2100A81-AF39-EE43-A57A-F2CCA607A38C.root"]
#fileList = ["root://xrootd-cms.infn.it///store/mc/RunIIAutumn18NanoAODv5/GluGluHToTauTau_M125_13TeV_powheg_pythia8/NANOAODSIM/Nano1June2019_102X_upgrade2018_realistic_v19-v1/130000/3D93C172-4F4A-3347-A6C0-720890E6CF19.root"]
#fileList = ["root://cms-xrd-global.cern.ch///store/mc/RunIISummer16NanoAODv5/HWplusJ_HToWW_M125_13TeV_powheg_pythia8/NANOAODSIM/PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1/40000/3EA28B4A-AA3E-B64B-A089-229A6840685A.root"]
#fileList = ["root://cms-xrd-global.cern.ch///store/mc/RunIIFall17NanoAODv5/HWplusJ_HToWW_M125_13TeV_powheg_pythia8_TuneCP5/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/F5A8EF96-FC78-DB4A-9BE7-56721DB6529A.root"]
#fileList = ["root://cms-xrd-global.cern.ch///store/data/Run2017E/Charmonium/NANOAOD/02Apr2020-v1/50000/6A6A7F98-6A6D-7348-8413-0BD676A8A3B3.root"]
#fileList = ["root://xrootd-cms.infn.it///store/data/Run2018A/MuonEG/NANOAOD/02Apr2020-v1/50000/E19DBB1E-52B2-6C48-B018-896CC218F744.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAODv7/ZToJPsiGamma_TuneCUETP8M1_13TeV-pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/260000/2101C93F-F83F-2949-88CE-6394F6DC5228.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAODv7/HToJPsiG_ToMuMuG_allProdMode_M125_TuneCUETP8M1_13TeV_Pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/270000/0ED3C659-F150-7C46-A562-188AF8F7400A.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/data/Run2018C/Charmonium/NANOAOD/UL2018_MiniAODv1_NanoAODv2-v1/280000/6121CB56-A6FD-4C4A-B34D-D770DAC84C4A.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAODv7/ZZTo4L_13TeV_powheg_pythia8_ext1/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/120000/D9C9041F-2D29-234F-8469-84268E2CE8E6.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/110000/2EA04007-C18B-0947-AD3F-A45BDFFCC287.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv2/ZH_HToJPsiG_JPsiToMuMu_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_mc2017_realistic_v8-v1/270000/41F44D62-F3AF-1F4B-93C1-9D81208080E4.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv2/GluGluToH_HToJPsiG_JPsiToMuMu_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_mc2017_realistic_v8-v1/280000/E95DA0F6-7458-354C-A087-773D7DDFC7AF.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv2/ZH_HToJPsiG_JPsiToMuMu_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/280000/87E2D29A-CCFA-274D-AA42-E42DE0A7132B.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv2/WH_HToJPsiG_JPsiToMuMu_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/120000/DB5BD773-328F-A645-B3E7-76EECAA77555.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAODv7/ZZTo4L_13TeV_powheg_pythia8/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/110000/BB59BAF3-9039-9F46-94C9-750B8035A3CC.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18NanoAODv7/GluGluHToMuMuG_M125_MLL-0To60_Dalitz_012j_13TeV_amcatnloFXFX_pythia8_PSWeight/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/40622801-4FD0-5C43-B23E-CF9D9A2CB0F0.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv2/ZH_HToJPsiG_JPsiToMuMu_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/280000/87E2D29A-CCFA-274D-AA42-E42DE0A7132B.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv2/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_mc2017_realistic_v8-v1/230000/7F690E06-A6C5-1F4B-8265-9CB811191A2D.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16NanoAODAPVv9/JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v1/250000/A8579E47-F328-0743-8B12-25C3C4E11982.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/data/Run2016F/SingleElectron/NANOAOD/UL2016_MiniAODv1_NanoAODv2-v2/230000/8CFB53D5-FC17-5240-A0DE-17A31CDDE7D7.root"]
fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16NanoAODv9/GluGluToH_HToJPsiG_JPsiToMuMu_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v17-v1/270000/EE1DA73C-B9E0-1846-8E5E-783784AC1C7F.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/GluGluHToMuMuG_M125_Dalitz_012j_TuneCP5_13TeV_amcatnlo_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/2820000/200E14A4-9A65-3B4A-9210-00378775AA7D.root"]
#fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18NanoAODv7/VBFHToMuMuG_M125_MLL-0To60_TuneCP5_Dalitz_13TeV_amcatnlo_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/120000/29831CA7-465E-1540-8665-FB09B06E2D2E.root"]


jsonFile, muonRoccoR = None, muonScaleRes2018
if len(fileList) == 0:
    print "- Filelist is empty"
    exit()
#elif "Run2016" in fileList[0]: jsonFile = MAINDIR + 'json/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt' # 36.773 /fb
#elif "Run2017" in fileList[0]: jsonFile = MAINDIR + 'json/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt' # 41.86 /fb
#elif "Run2018" in fileList[0]: jsonFile = MAINDIR + 'json/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt' # 58.83 /fb
elif "Run2016" in fileList[0]:
    jsonFile = MAINDIR + 'json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt' # 35.93 /fb
    muonRoccoR = muonScaleRes2016
elif "Run2017" in fileList[0]:
    jsonFile = MAINDIR + 'json/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt' # 41.48 /fb
    muonRoccoR = muonScaleRes2017
elif "Run2018" in fileList[0]:
    jsonFile = MAINDIR + 'json/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt' # 59.83 /fb
    muonRoccoR = muonScaleRes2018

#p = PostProcessor(outputDir=options.output, inputFiles=fileList, cut=None, branchsel=None, modules=[SameSign()], jsonInput=jsonFile, histFileName=None, histDirName=None, outputbranchsel=MAINDIR + "keep_and_drop.txt", maxEntries=long(options.maxEntries))

#p = PostProcessor(outputDir=options.output, inputFiles=fileList, cut="(nMuon + nElectron >= 4)", branchsel=None, modules=[HZZ4L()], jsonInput=jsonFile, histFileName=None, histDirName=None, outputbranchsel=MAINDIR + "keep_and_drop_higgs.txt", maxEntries=long(options.maxEntries))

p = PostProcessor(outputDir=options.output, inputFiles=fileList, cut=(None if 'JPsiG' in fileList[0] else "(nMuon >= 2 && nPhoton >= 1)"), branchsel=None, modules=[muonRoccoR(),Higgs()], jsonInput=jsonFile, histFileName=None, histDirName=None, outputbranchsel=MAINDIR + "keep_and_drop_higgs.txt", maxEntries=long(options.maxEntries))

#p = PostProcessor(outputDir=options.output, inputFiles=fileList, cut="(nMuon>=2 && Muon_pt[1]>5.)", branchsel=None, modules=[TriggerJPsi()], jsonInput=jsonFile, histFileName=None, histDirName=None, outputbranchsel=MAINDIR + "keep_and_drop_trigger.txt", maxEntries=long(options.maxEntries))

#p = PostProcessor(outputDir=options.output, inputFiles=fileList, cut="(nElectron>=1 && Electron_pt[0]>25. && nMuon>=1 && Muon_pt[0]>10. && nPhoton>=1 && Photon_pt[0]>15.)", branchsel=None, modules=[TriggerMuPhoton()], jsonInput=jsonFile, histFileName=None, histDirName=None, outputbranchsel=MAINDIR + "keep_and_drop_trigger.txt", maxEntries=long(options.maxEntries))

#p = PostProcessor(outputDir=options.output, inputFiles=fileList, cut="(nMuon>=2 && Muon_pt[0]>30. && Muon_pt[1]>10. && nPhoton>=1 && Photon_pt[0]>15.)", branchsel=None, modules=[TriggerPhoton()], jsonInput=jsonFile, histFileName=None, histDirName=None, outputbranchsel=MAINDIR + "keep_and_drop_trigger.txt", maxEntries=long(options.maxEntries))

#p = PostProcessor(outputDir=options.output, inputFiles=fileList, cut="(nElectron>=1 && Electron_pt[0]>30. && nMuon>=1 && Muon_pt[0]>10.)", branchsel=None, modules=[TriggerMu()], jsonInput=jsonFile, histFileName=None, histDirName=None, outputbranchsel=MAINDIR + "keep_and_drop_trigger.txt", maxEntries=long(options.maxEntries))

p.run()

print "+ Done."


#def __init__(self,outputDir,inputFiles,cut=None,branchsel=None,modules=[],compression="LZMA:9",friend=False,postfix=None,jsonInput=None,noOut=False,justcount=False,provenance=False,haddFileName=None,fwkJobReport=False,histFileName=None,histDirName=None,outputbranchsel=None,maxEntries=None,firstEntry=0,prefetch=False,longTermCache=False):

# voms-proxy-init --rfc --voms cms --valid 168:00

