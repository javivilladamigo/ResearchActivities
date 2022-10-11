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

fileList = ["/lustre/cmsdata/zucchett/Dataset/QCD_Pt-30_MuEnrichedPt4_TuneCP5_13TeV_pythia8_RunIISummer20UL16/4DA5E108-3E3A-E443-A4A4-B6A16A7F265C.root"]
# GluGluToH_HToJPsi     ->   /lustre/cmsdata/zucchett/Dataset/GluGluToH_HToJPsiG_JPsiToMuMu_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16/147A56C1-E191-3842-84BB-81A243C2A4AD.root
# ZToJPsiG              ->   /lustre/cmsdata/zucchett/Dataset/ZToJPsiG_JPsiToMuMu_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16/3CFD568D-76B9-3D4D-8750-5A29F1E2C278.root
# QCD                   ->   /lustre/cmsdata/zucchett/Dataset/QCD_Pt-30_MuEnrichedPt4_TuneCP5_13TeV_pythia8_RunIISummer20UL16/BB9CEC64-DAB0-424C-B36E-840FB52120F3.root

dir_file = "/lustre/cmsdata/zucchett/Dataset/QCD_Pt-30_MuEnrichedPt4_TuneCP5_13TeV_pythia8_RunIISummer20UL16/"
fileList = [os.path.join(dir_file, f) for f in os.listdir(dir_file) if (os.path.isfile(os.path.join(dir_file, f)) and ".root" in f)]

p = PostProcessor(outputDir=options.output, inputFiles=fileList, cut=(None if 'JPsiG' in fileList[0] else "(nMuon >= 2 && nPhoton >= 1)"), branchsel=None, modules=[Higgs_marino()], jsonInput=None, histFileName=None, histDirName=None, outputbranchsel="keep_and_drop_higgs_marino.txt", maxEntries=long(options.maxEntries))

p.run()

print "+ Done."