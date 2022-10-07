#!/bin/env python

import os, sys
from global_paths import *
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *

from higgs_marino import *


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

fileList = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/GluGluToH_HToJPsiG_JPsiToMuMu_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/70000/FF49F3B9-055A-7942-8E5B-2B3028B33CB6.root"]
# GluGluToH_HToJPsi     ->   /store/mc/RunIISummer20UL18NanoAODv9/GluGluToH_HToJPsiG_JPsiToMuMu_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/70000/FF49F3B9-055A-7942-8E5B-2B3028B33CB6.root
# QCD                   ->   /store/mc/RunIISummer20UL16NanoAODv9/QCD_Pt-30_MuEnrichedPt4_TuneCP5_13TeV_pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v17-v1/130000/4DA5E108-3E3A-E443-A4A4-B6A16A7F265C.root
p = PostProcessor(outputDir=options.output, inputFiles=fileList, cut=(None if 'JPsiG' in fileList[0] else "(nMuon >= 2 && nPhoton >= 1)"), branchsel=None, modules=[Higgs_marino()], jsonInput=None, histFileName=None, histDirName=None, outputbranchsel="keep_and_drop_higgs_marino.txt", maxEntries=long(options.maxEntries))

p.run()

print "+ Done."