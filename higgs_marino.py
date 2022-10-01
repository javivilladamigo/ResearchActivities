#!/usr/bin/env python
import os, sys, math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *

#from CorrectionTools.PileupWeightTool import *
#from CorrectionTools.BTaggingTool import BTagWeightTool, BTagWPs
#from CorrectionTools.MuonSFs import *
#from CorrectionTools.ElectronSFs import *
#from CorrectionTools.RecoilCorrectionTool import getTTptWeight, getTTPt
#from CorrectionTools.DYCorrection import *
#from PileupWeightTool import PileupWeightTool
#from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeightProducer
from CorrectionTools.PileupWeightTool import *
from CorrectionTools.MuonSFs import *
from CorrectionTools.PhotonSFs import *
#from metadataSF import *
from global_paths import *

from ROOT import TLorentzVector, TVector3, TVector2
from utils import EV, XS, getNameFromFile



class Higgs(Module):
    def __init__(self):
        self.writeHistFile=True
        


    def beginJob(self,histFile=None,histDirName=None):
        Module.beginJob(self,histFile,histDirName)
        self.event = 0
        self.hists = {}
        self.hists["Nevents"] = ROOT.TH1F("Nevents", "Nevents", 1, 0, 1)
        self.hists["Acceptance"] = ROOT.TH1F("Acceptance", "Acceptance", 5, -0.5, 4.5)
        self.hists["genH_pt"] = ROOT.TH1F("genH_pt", ";generator H p_{T} (GeV)", 100, 0., 100.)
        self.hists["genH_eta"] = ROOT.TH1F("genH_eta", ";generator H #eta", 100, -5., 5.)
        self.hists["genH_phi"] = ROOT.TH1F("genH_phi", ";generator H #phi", 100, -3.15, 3.15)
        self.hists["genH_mass"] = ROOT.TH1F("genH_mass", ";generator H mass", 300, 0., 150.)
        self.hists["genJPsi_pt"] = ROOT.TH1F("genJPsi_pt", ";generator J/#Psi p_{T} (GeV)", 100, 0., 100.)
        self.hists["genJPsi_eta"] = ROOT.TH1F("genJPsi_eta", ";generator J/#Psi #eta", 100, -5., 5.)
        self.hists["genJPsi_phi"] = ROOT.TH1F("genJPsi_phi", ";generator J/#Psi #phi", 100, -3.15, 3.15)
        self.hists["genJPsi_mass"] = ROOT.TH1F("genJPsi_mass", ";generator J/#Psi mass", 200, 0., 4.)
        self.hists["genPhoton_pt"] = ROOT.TH1F("genPhoton_pt", ";generator #gamma p_{T} (GeV)", 100, 0., 100.)
        self.hists["genPhoton_eta"] = ROOT.TH1F("genPhoton_eta", ";generator #gamma #eta", 100, -5., 5.)
        self.hists["genMuon1_pt"] = ROOT.TH1F("genMuon1_pt", ";generator #mu^{-} p_{T} (GeV)", 100, 0., 100.)
        self.hists["genMuon1_eta"] = ROOT.TH1F("genMuon1_eta", ";generator #mu^{-} #eta", 100, -5., 5.)
        self.hists["genMuon2_pt"] = ROOT.TH1F("genMuon2_pt", ";generator #mu^{+} p_{T} (GeV)", 100, 0., 100.)
        self.hists["genMuon2_eta"] = ROOT.TH1F("genMuon2_eta", ";generator #mu^{+} #eta", 100, -5., 5.)
        
        self.hists["genCosThetaStar"] = ROOT.TH1F("genCosThetaStar", ";cos #theta^{*}", 100, -1., 1.)
        self.hists["genCosTheta1"] = ROOT.TH1F("genCosTheta1", ";cos #theta_{1}", 100, -1., 1.)
        self.hists["genPhi1"] = ROOT.TH1F("genPhi1", ";#Phi_{1}", 100, -3.1415, 3.1415)
        self.hists["genCosThetaStarZtoMM"] = ROOT.TH1F("genCosThetaStarZtoMM", ";cos #theta^{*}", 100, -1., 1.)
        
        self.hists["central"] = ROOT.TH1F("central", "", 125, 76, 200)
        self.hists["phScaleUp"] = ROOT.TH1F("phScaleUp", "", 125, 76, 200)
        self.hists["phScaleDown"] = ROOT.TH1F("phScaleDown", "", 125, 76, 200)
        self.hists["phResUp"] = ROOT.TH1F("phResUp", "", 125, 76, 200)
        self.hists["phResDown"] = ROOT.TH1F("phResDown", "", 125, 76, 200)
        
        self.hists["Cutflow"] = ROOT.TH1F("Cutflow", "Cutflow", 15, -0.5, 14.5)
        self.hists["Cutflow"].GetXaxis().SetBinLabel(1, "All events")
        self.hists["Cutflow"].GetXaxis().SetBinLabel(2, "Acceptance")
        self.hists["Cutflow"].GetXaxis().SetBinLabel(3, "2 reco muons")
        self.hists["Cutflow"].GetXaxis().SetBinLabel(4, "J/#Psi cand")
        self.hists["Cutflow"].GetXaxis().SetBinLabel(5, "reco #gamma")
        self.hists["Cutflow"].GetXaxis().SetBinLabel(6, "preselections")
        self.hists["Cutflow"].GetXaxis().SetBinLabel(7, "muon id")
        self.hists["Cutflow"].GetXaxis().SetBinLabel(8, "muon iso")
        self.hists["Cutflow"].GetXaxis().SetBinLabel(9, "#gamma id")
        self.hists["Cutflow"].GetXaxis().SetBinLabel(10, "#gamma pixel veto")
        self.hists["Cutflow"].GetXaxis().SetBinLabel(11, "trigger kin")
        self.hists["Cutflow"].GetXaxis().SetBinLabel(12, "trigger")
        self.hists["Cutflow"].GetXaxis().SetBinLabel(13, "signal region")
        
        self.hists["Matching"] = ROOT.TH1F("Matching", "Matching", 10, -0.5, 9.5)
        self.hists["Matching"].GetXaxis().SetBinLabel(1, "all events")
        self.hists["Matching"].GetXaxis().SetBinLabel(2, "gen match")
        self.hists["Matching"].GetXaxis().SetBinLabel(3, "in acceptance")
        self.hists["Matching"].GetXaxis().SetBinLabel(4, "2 reco muons")
        self.hists["Matching"].GetXaxis().SetBinLabel(5, "2 matched muons")
        self.hists["Matching"].GetXaxis().SetBinLabel(6, "reco #gamma")
        self.hists["Matching"].GetXaxis().SetBinLabel(7, "matched #gamma")
        self.hists["Matching"].GetXaxis().SetBinLabel(8, "All match")
        
        
#        self.addObject(self.h_events)
#        self.h_events.SetDirectory
        self.verbose = 0
    
    def endJob(self):
        Module.endJob(self)
        if self.verbose >= 0: print "+ Module ended successfully."#, self.h_events.GetEntries(), "events analyzed"
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("isMC", "I")
        self.out.branch("is2016", "I")
        self.out.branch("is2017", "I")
        self.out.branch("is2018", "I")
        self.out.branch("isSingleMuonTrigger", "I")
        self.out.branch("isSingleMuonPhotonTrigger", "I")
        self.out.branch("isSingleMuonNoFiltersPhotonTrigger", "I")
        self.out.branch("isDoubleMuonTrigger", "I")
        self.out.branch("isDoubleMuonPhotonTrigger", "I")
        self.out.branch("isJPsiTrigger", "I")
        self.out.branch("isDisplacedTrigger", "I")
        self.out.branch("passedMETFilters", "I")
        self.out.branch("nCleanElectron", "I")
        self.out.branch("nCleanMuon", "I")
        #self.out.branch("nCleanTau", "I")
        self.out.branch("nCleanPhoton", "I")
        #self.out.branch("nCleanCentralJet", "I")
        #self.out.branch("nCleanForwardJet", "I")
        #self.out.branch("nCleanBTagJet", "I")
        self.out.branch("HT30", "F")
        self.out.branch("iPhoton", "I")
        self.out.branch("iMuon1", "I")
        self.out.branch("iMuon2", "I")
        #
#        self.out.branch("isINC", "I")
#        self.out.branch("isVBF", "I")
#        self.out.branch("isHFL", "I")
#        self.out.branch("isLEP", "I")
        #
        self.out.branch("Muon1_pt", "F")
        self.out.branch("Muon1_eta", "F")
        self.out.branch("Muon2_pt", "F")
        self.out.branch("Muon2_eta", "F")
        self.out.branch("Muon1_pfRelIso03", "F")
        self.out.branch("Muon2_pfRelIso03", "F")
        self.out.branch("Muon1_pfRelIso04", "F")
        self.out.branch("Muon2_pfRelIso04", "F")
        self.out.branch("Muon1_mediumId", "I")
        self.out.branch("Muon2_mediumId", "I")
        self.out.branch("Muon1_mediumPromptId", "I")
        self.out.branch("Muon2_mediumPromptId", "I")
        self.out.branch("Muon1_ip3d", "F")
        self.out.branch("Muon2_ip3d", "F")
        self.out.branch("minMuonPfIso", "F")
        self.out.branch("maxMuonPfIso", "F")
        self.out.branch("minMuonPfSubIso", "F")
        self.out.branch("maxMuonPfSubIso", "F")
        self.out.branch("minMuonTrkSubIso", "F")
        self.out.branch("maxMuonTrkSubIso", "F")
        self.out.branch("Muon12_diffdxy", "F")
        self.out.branch("Muon12_diffdz", "F")
        self.out.branch("Muon12_signdxy", "F")
        self.out.branch("Muon12_signdz", "F")
        self.out.branch("Photon1_pt", "F")
        self.out.branch("Photon1_eta", "F")
        self.out.branch("Photon1_mvaID_WP80", "I")
        self.out.branch("Photon1_mvaID_WP90", "I")
        self.out.branch("Photon1_r9", "F")
        self.out.branch("Photon1_pfRelIso03", "F")
        self.out.branch("Photon1_pixelSeed", "I")
        
        self.out.branch("JPsi_pt", "F")
        self.out.branch("JPsi_eta", "F")
        self.out.branch("JPsi_phi", "F")
        self.out.branch("JPsi_mass", "F")
        self.out.branch("JPsi_dEta", "F")
        self.out.branch("JPsi_dPhi", "F")
        self.out.branch("JPsi_dR", "F")
        self.out.branch("H_pt", "F")
        self.out.branch("H_eta", "F")
        self.out.branch("H_phi", "F")
        self.out.branch("H_mass", "F")
        self.out.branch("H_dEta", "F")
        self.out.branch("H_dPhi", "F")
        self.out.branch("H_dR", "F")
        
        self.out.branch("JPsi_mass_muScaleUp", "F")
        self.out.branch("JPsi_mass_muScaleDown", "F")
        self.out.branch("H_mass_muScaleUp", "F")
        self.out.branch("H_mass_muScaleDown", "F")
        #self.out.branch("H_mass_phScaleUp", "F")
        #self.out.branch("H_mass_phScaleDown", "F")
        self.out.branch("H_mass_phResUp", "F")
        self.out.branch("H_mass_phResDown", "F")
        self.out.branch("H_mass_noCorr", "F")
        
        self.out.branch("minMuonMetDPhi", "F")
        self.out.branch("maxMuonMetDPhi", "F")
        self.out.branch("photonMetDPhi", "F")
        self.out.branch("metPlusPhotonDPhi", "F")
        #self.out.branch("vbfMjj", "F")
        #self.out.branch("vbfDR", "F")
        #self.out.branch("centrality", "F")
        #self.out.branch("trCentr", "F")
        #self.out.branch("cosThetaStar", "F")
        #self.out.branch("cosTheta1", "F")
        #self.out.branch("phi1", "F")
        #self.out.branch("genCosThetaStar", "F")
        #self.out.branch("genCosTheta1", "F")
        #self.out.branch("genPhi1", "F")
        self.out.branch("lumiWeight", "F")
        self.out.branch("lheWeight", "F")
        self.out.branch("stitchWeight", "F")
        #self.out.branch("topWeight", "F")
        #self.out.branch("qcdnloWeight", "F")
        #self.out.branch("qcdnnloWeight", "F")
        #self.out.branch("ewknloWeight", "F")
        #self.out.branch("puWeight", "F")
        #self.out.branch("puWeightUp", "F")
        #self.out.branch("puWeightDown", "F")
        self.out.branch("triggerWeight", "F")
        #self.out.branch("triggerWeightUp", "F")
        #self.out.branch("triggerWeightDown", "F")
        #self.out.branch("triggerMuonEGWeight", "F")
        #self.out.branch("triggerMuonEGWeightUp", "F")
        #self.out.branch("triggerMuonEGWeightDown", "F")
        self.out.branch("triggerJPsiWeight", "F")
        self.out.branch("triggerJPsiWeightUp", "F")
        self.out.branch("triggerJPsiWeightDown", "F")
        #self.out.branch("leptonWeight", "F")
        #self.out.branch("leptonWeightUp", "F")
        #self.out.branch("leptonWeightDown", "F")
        #self.out.branch("photonWeight", "F")
        #self.out.branch("photonWeightUp", "F")
        #self.out.branch("photonWeightDown", "F")
        #self.out.branch("eventMuonEGWeightLumi", "F")
        #self.out.branch("eventJPsiWeightLumi", "F")
        self.out.branch("eventWeightLumi", "F")
        
        self.fileName = inputFile.GetName()
        self.sampleName = getNameFromFile(self.fileName)
        
        self.isMC = not "Run201" in self.fileName
        
        self.N = 0
        
        if self.verbose >= 0: print "+ Opening file", self.fileName
        
        self.thMuons = 5.
        self.thPhoton = 10.
        self.thMatching = 0.01
        self.muonCriteria = "DR"
        
        # b-tagging working points for DeepCSV
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
        if "Run2016" in self.fileName or "UL16" in self.fileName or "Summer16" in self.fileName:
            self.year = 2016
            self.lumi = 35930.
            self.btagLoose = 0.2217 #0.0614
            self.btagMedium = 0.6321 #0.3093
            self.btagTight = 0.8953 #0.7221
        elif "Run2017" in self.fileName or "UL17" in self.fileName or "Fall17" in self.fileName:
            self.year = 2017
            self.lumi = 41480.
            self.btagLoose = 0.1522 #0.0521
            self.btagMedium = 0.4941 #0.3033
            self.btagTight = 0.8001 #0.7489
        elif "Run2018" in self.fileName or "UL18" in self.fileName or "Autumn18" in self.fileName:
            self.year = 2018
            self.lumi = 59830.
            self.btagLoose = 0.1241 #0.0494
            self.btagMedium = 0.4184 #0.2770
            self.btagTight = 0.7527 #0.7264
        else:
            if self.verbose >= 0: print "- Unknown year, aborting module"
            sys.exit()
        
        if not self.sampleName in XS:
            if self.verbose >= 0: print "- No XS available for sample ", self.sampleName, ", aborting"
            sys.exit()
        if not self.sampleName in EV:
            if self.verbose >= 0: print "- No number of events available for sample ", self.sampleName, ", aborting"
            sys.exit()
        self.xs = XS[self.sampleName]
        self.nevents = EV[self.sampleName]
        self.xsWeight = self.xs / self.nevents
        self.lumiWeight = self.xsWeight * self.lumi if self.isMC else 1.
        self.isLO = self.nevents < 1e9 and abs(self.nevents % 1) < 1.e-6 # if True, the event count is integer, so the weight should be normalized (+1)
        self.isSignal = ('JPsiG' in self.sampleName or 'PsiPrime' in self.sampleName)
        self.isZSignal = ('ZToJPsiG' in self.sampleName or 'ZToPsiPrime' in self.sampleName)
        self.isHSignal = ('HToJPsiG' in self.sampleName or 'HToPsiPrime' in self.sampleName)
        
        if self.verbose >= 0: print "+ Module parameters: isMC", self.isMC, ", isLO", self.isLO, ", year", self.year, ", lumi", self.lumi, "pb"
        if self.verbose >= 0: print "+ Sample:", self.sampleName, ", XS:", self.xs, ", events:", self.nevents
        if self.verbose >= 0: print "+ LumiWeight:", self.lumiWeight
        #if self.isMC and self.isLO and self.verbose >= 0: print "+ Sample is LO, gen weight will be set to 1"
#        self.puTool = PileupWeightTool(year = year) if self.isMC else None

        #self.SingleMuonTriggers = ["HLT_IsoMu27"]
        self.SingleMuonPhotonTriggers = ["HLT_Mu17_Photon30_CaloIdL_L1ISO", "HLT_Mu17_Photon30_IsoCaloId", "HLT_Mu17_Photon30_CaloIdL"] # 27.13 in 2017
        #self.SingleMuonNoFiltersPhotonTriggers = ["HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL", "HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL", "HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL"]
        #self.DoubleMuonTriggers = ["HLT_Mu17_Mu8", "HLT_Mu17_Mu8_DZ", "HLT_Mu17_TkMu8_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass", "HLT_Mu37_TkMu27", ] #, "HLT_DoubleMu33NoFiltersNoVtxDisplaced"]
        #self.DoubleMuonPhotonTriggers = ["HLT_DoubleMu20_7_Mass0to30_Photon23"]
        #self.JPsiTriggers = ["HLT_Dimuon16_Jpsi", "HLT_Dimuon18_PsiPrime", "HLT_Dimuon18_PsiPrime_noCorrL1", "HLT_Dimuon25_Jpsi", "HLT_Dimuon25_Jpsi_noCorrL1", "HLT_Dimuon20_Jpsi", ]
        #self.DisplacedTriggers = ["HLT_DoubleMu4_JpsiTrk_Displaced", "HLT_DoubleMu4_PsiPrimeTrk_Displaced"]
        
        if self.isMC:
            self.muSFs  = None #MuonSFs(year = self.year)
            self.elSFs  = None #ElectronSFs(year = self.year)
            self.puTool = PileupWeightTool(year = self.year)
        
        # trigger SF
        self.triggerSFdict = {
        #    "SingleMuonPhoton" : {
        #        2016 : 0.982,
        #        2017 : 0.919,
        #        2018 : 0.972,
        #    },
            "SingleMuonPhotonUp" : {
                2016 : 0.982 + 0.059,
                2017 : 0.919 + 0.061,
                2018 : 0.972 + 0.035,
            },
            "SingleMuonPhotonDown" : {
                2016 : 0.982 - 0.059,
                2017 : 0.919 - 0.061,
                2018 : 0.972 - 0.035,
            },
        #    "JPsi" : {
        #        2016 : 1.004,
        #        2017 : 1.021,
        #        2018 : 1.006,
        #    },
        #    "JPsiUp" : {
        #        2016 : 1.004 + 0.032,
        #        2017 : 1.021 + 0.010,
        #        2018 : 1.006 + 0.008,
        #    },
        #    "JPsiDown" : {
        #        2016 : 1.004 - 0.032,
        #        2017 : 1.021 - 0.010,
        #        2018 : 1.006 - 0.008,
        #    },
        }
            
        # muon SF tools
        self.muSFdict = {
            "loose" : {
                2016 : MuonSFs(MUONSFROOTFILES, MUONSFHISTNAMES["loose"],  MUONSFIS_PTVSETA, year=2016, isPreVFP=False),
                2017 : MuonSFs(MUONSFROOTFILES, MUONSFHISTNAMES["loose"],  MUONSFIS_PTVSETA, year=2017, isPreVFP=False),
                2018 : MuonSFs(MUONSFROOTFILES, MUONSFHISTNAMES["loose"],  MUONSFIS_PTVSETA, year=2018, isPreVFP=False),
            },
            "tight" : {
                2016 : MuonSFs(MUONSFROOTFILES, MUONSFHISTNAMES["tight"],  MUONSFIS_PTVSETA, year=2016, isPreVFP=False),
                2017 : MuonSFs(MUONSFROOTFILES, MUONSFHISTNAMES["tight"],  MUONSFIS_PTVSETA, year=2017, isPreVFP=False),
                2018 : MuonSFs(MUONSFROOTFILES, MUONSFHISTNAMES["tight"],  MUONSFIS_PTVSETA, year=2018, isPreVFP=False),
            },
            "prompt": {
                2016 : MuonSFs(MUONSFROOTFILES, MUONSFHISTNAMES["prompt"], MUONSFIS_PTVSETA, year=2016, isPreVFP=False),
                2017 : MuonSFs(MUONSFROOTFILES, MUONSFHISTNAMES["prompt"], MUONSFIS_PTVSETA, year=2017, isPreVFP=False),
                2018 : MuonSFs(MUONSFROOTFILES, MUONSFHISTNAMES["prompt"], MUONSFIS_PTVSETA, year=2018, isPreVFP=False),
            },
        }

        # photon SF tools
        self.phoSFdict = {
            "MVA80" : {
                2016 : PhotonSFs(PHOTONSFROOTFILES["MVA80"], PHOTONSFHISTNAMES, PHOTONSFIS_PTVSETA, year=2016, isPreVFP=False),
                2017 : PhotonSFs(PHOTONSFROOTFILES["MVA80"], PHOTONSFHISTNAMES, PHOTONSFIS_PTVSETA, year=2017, isPreVFP=False),
                2018 : PhotonSFs(PHOTONSFROOTFILES["MVA80"], PHOTONSFHISTNAMES, PHOTONSFIS_PTVSETA, year=2018, isPreVFP=False),
            },
            "MVA90" : {
                2016 : PhotonSFs(PHOTONSFROOTFILES["MVA90"], PHOTONSFHISTNAMES, PHOTONSFIS_PTVSETA, year=2016, isPreVFP=False),
                2017 : PhotonSFs(PHOTONSFROOTFILES["MVA90"], PHOTONSFHISTNAMES, PHOTONSFIS_PTVSETA, year=2017, isPreVFP=False),
                2018 : PhotonSFs(PHOTONSFROOTFILES["MVA90"], PHOTONSFHISTNAMES, PHOTONSFIS_PTVSETA, year=2018, isPreVFP=False),
            },
        }

        # pileup tools
        self.puWeightSFdict = {
            "central" : {
                2016 : PileupWeightTool(year=2016, sigma="central"),
                2017 : PileupWeightTool(year=2017, sigma="central"),
                2018 : PileupWeightTool(year=2018, sigma="central"),
            },
            "up" : {
                2016 : PileupWeightTool(year=2016, sigma="up"),
                2017 : PileupWeightTool(year=2017, sigma="up"),
                2018 : PileupWeightTool(year=2018, sigma="up"),
            },
            "down" : {
                2016 : PileupWeightTool(year=2016, sigma="down"),
                2017 : PileupWeightTool(year=2017, sigma="down"),
                2018 : PileupWeightTool(year=2018, sigma="down"),
            },
        }
            
            
    
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        print "Efficiency:\t%.4f" % ( (self.hists["Cutflow"].GetBinContent(13) / self.hists["Cutflow"].GetBinContent(1)) if self.hists["Cutflow"].GetBinContent(1) > 0 else 0., )
        #
        outputFile.mkdir("Hists")
        outputFile.cd("Hists")
        for histname, hist in self.hists.iteritems():
            hist.Write()
        outputFile.cd("..")
        if self.verbose >= 0: print "+ File closed successfully, written", self.N, "events to tree"
        pass
        
    
    def analyze(self, event):
        eventWeightLumi, eventMuonEGWeightLumi, eventJPsiWeightLumi, lheWeight, stitchWeight, qcdnloWeight, qcdnnloWeight, ewknloWeight, topWeight = 1., 1., 1., 1., 1., 1., 1., 1., 1.
        puWeight, puWeightUp, puWeightDown, leptonWeight, leptonWeightUp, leptonWeightDown, photonWeight, photonWeightUp, photonWeightDown = 1., 1., 1., 1., 1., 1., 1., 1., 1.
        triggerWeight, triggerWeightUp, triggerWeightDown, triggerMuonEGWeight, triggerMuonEGWeightUp, triggerMuonEGWeightDown, triggerJPsiWeight, triggerJPsiWeightUp, triggerJPsiWeightDown = 1., 1., 1., 1., 1., 1., 1., 1., 1.
        isSingleMuonTrigger, isSingleMuonPhotonTrigger, isSingleMuonNoFiltersPhotonTrigger, isDoubleMuonTrigger, isDoubleMuonPhotonTrigger, isJPsiTrigger, isDisplacedTrigger = False, False, False, False, False, False, False
        nCleanElectron, nCleanMuon, nCleanTau, nCleanPhoton, nCleanCentralJet, nCleanForwardJet, nCleanBTagJet, HT30 = 0, 0, 0, 0, 0, 0, 0, 0
        cosThetaStar, cosTheta1, phi1 = -2., -2., -4.
        genCosThetaStar, genCosTheta1, genPhi1 = -2., -2., -4.
        centrality, trCentr = -1., -1.

        for t in self.SingleMuonTriggers:
            if hasattr(event, t) and getattr(event, t): isSingleMuonTrigger = True
        for t in self.SingleMuonPhotonTriggers:
            if hasattr(event, t) and getattr(event, t): isSingleMuonPhotonTrigger = True
        for t in self.SingleMuonNoFiltersPhotonTriggers:
            if hasattr(event, t) and getattr(event, t): isSingleMuonNoFiltersPhotonTrigger = True
        for t in self.DoubleMuonTriggers:
            if hasattr(event, t) and getattr(event, t): isDoubleMuonTrigger = True
        for t in self.DoubleMuonPhotonTriggers:
            if hasattr(event, t) and getattr(event, t): isDoubleMuonPhotonTrigger = True
        for t in self.JPsiTriggers:
            if hasattr(event, t) and getattr(event, t): isJPsiTrigger = True
        for t in self.DisplacedTriggers:
            if hasattr(event, t) and getattr(event, t): isDisplacedTrigger = True
        
        lheWeight = 1.
        
        if self.isMC: 
            # Event weight
            if not self.isLO and hasattr(event, "LHEWeight_originalXWGTUP"): lheWeight = event.LHEWeight_originalXWGTUP
            else: lheWeight = 1.
            
            # PU weight
            #puWeight = self.puTool.getWeight(event.Pileup_nTrueInt)
            
            if self.verbose >= 1: print "+ Event LHE Weight:", lheWeight
            if self.verbose >= 1: print "+ Event Weight:", lheWeight * self.lumiWeight
        
        self.hists["Nevents"].Fill(0, lheWeight)
        self.hists["Acceptance"].Fill(0, lheWeight)
        self.hists["Cutflow"].Fill(0, lheWeight)
        
        genHIdx, genJPsiIdx, genMuon1Idx, genMuon2Idx, genPhotonIdx = -1, -1, -1, -1, -1
        genH, genJPsi, genMuon1, genMuon2, genPhoton = TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector()
        isGenMatched, isGenAcceptance, nAddGenMuons = False, False, 0
        genMuonsMatching, genPhotonMatching, genEventMatching = False, False, False
        # Gen studies
        if self.isMC and self.isSignal:
#            print "-"*40
            self.hists["Matching"].Fill(0, lheWeight)
            for i in range(event.nGenPart):
#                print i, "\t", event.GenPart_pdgId[i], "\t", event.GenPart_status[i], "\t", event.GenPart_statusFlags[i], "\t", event.GenPart_pt[i]
                if event.GenPart_pdgId[i] == 25 or event.GenPart_pdgId[i] == 23: genHIdx = i
                if event.GenPart_pdgId[i]  in [443, 100443] and event.GenPart_genPartIdxMother[i] >= 0 and event.GenPart_pdgId[event.GenPart_genPartIdxMother[i]] in [23, 25]: genJPsiIdx = i
                if event.GenPart_pdgId[i] == 22 and event.GenPart_status[i] in [1, 23] and event.GenPart_genPartIdxMother[i] >= 0 and event.GenPart_pdgId[event.GenPart_genPartIdxMother[i]] in [23, 25]: genPhotonIdx = i #and (genPhotonIdx < 0 or event.GenPart_pt[i] > event.GenPart_pt[genPhotonIdx])
                if event.GenPart_pdgId[i] == -13 and event.GenPart_status[i] in [1, 23] and event.GenPart_genPartIdxMother[i] >= 0 and event.GenPart_pdgId[event.GenPart_genPartIdxMother[i]] in [443, 100443] and event.GenPart_genPartIdxMother[event.GenPart_genPartIdxMother[i]] >= 0 and event.GenPart_pdgId[event.GenPart_genPartIdxMother[event.GenPart_genPartIdxMother[i]]] in [23, 25, 443, 100443]: genMuon1Idx = i
                if event.GenPart_pdgId[i] == +13 and event.GenPart_status[i] in [1, 23] and event.GenPart_genPartIdxMother[i] >= 0 and event.GenPart_pdgId[event.GenPart_genPartIdxMother[i]] in [443, 100443] and event.GenPart_genPartIdxMother[event.GenPart_genPartIdxMother[i]] >= 0 and event.GenPart_pdgId[event.GenPart_genPartIdxMother[event.GenPart_genPartIdxMother[i]]] in [23, 25, 443, 100443]: genMuon2Idx = i 
                if abs(event.GenPart_pdgId[i]) in [13] and event.GenPart_status[i] in [1, 23] and event.GenPart_genPartIdxMother[i] >= 0 and event.GenPart_pdgId[event.GenPart_genPartIdxMother[i]] in [23, 24]: nAddGenMuons += 1
                
            #if nAddGenMuons < 1: return False # Only for associated production gen filtering
            if genHIdx >= 0: genH.SetPtEtaPhiM(event.GenPart_pt[genHIdx], event.GenPart_eta[genHIdx], event.GenPart_phi[genHIdx], event.GenPart_mass[genHIdx])
            if genJPsiIdx >= 0: genJPsi.SetPtEtaPhiM(event.GenPart_pt[genJPsiIdx], event.GenPart_eta[genJPsiIdx], event.GenPart_phi[genJPsiIdx], event.GenPart_mass[genJPsiIdx])
            if genPhotonIdx >= 0: genPhoton.SetPtEtaPhiM(event.GenPart_pt[genPhotonIdx], event.GenPart_eta[genPhotonIdx], event.GenPart_phi[genPhotonIdx], event.GenPart_mass[genPhotonIdx])
            if genMuon1Idx >= 0: genMuon1.SetPtEtaPhiM(event.GenPart_pt[genMuon1Idx], event.GenPart_eta[genMuon1Idx], event.GenPart_phi[genMuon1Idx], event.GenPart_mass[genMuon1Idx])
            if genMuon2Idx >= 0: genMuon2.SetPtEtaPhiM(event.GenPart_pt[genMuon2Idx], event.GenPart_eta[genMuon2Idx], event.GenPart_phi[genMuon2Idx], event.GenPart_mass[genMuon2Idx])
            
            if genHIdx >= 0 and genJPsiIdx >= 0 and genPhotonIdx >= 0 and genMuon1Idx >= 0 and genMuon2Idx >= 0:
                isGenMatched = True
                # Recalculate candidate 4-vectors for consistent calculation of angular variables
                genJPsi = genMuon1 + genMuon2
                genH = genJPsi + genPhoton
                
                genCosThetaStar, genCosTheta1, genPhi1 = self.returnCosThetaStar(genH, genJPsi), self.returnCosTheta1(genJPsi, genMuon1, genMuon2, genPhoton), self.returnPhi1(genH, genMuon1, genMuon2)
                
                self.hists["genH_pt"].Fill(genH.Pt(), lheWeight)
                self.hists["genH_eta"].Fill(genH.Eta(), lheWeight)
                self.hists["genH_phi"].Fill(genH.Phi(), lheWeight)
                self.hists["genH_mass"].Fill(genH.M(), lheWeight)
                self.hists["genJPsi_pt"].Fill(genJPsi.Pt(), lheWeight)
                self.hists["genJPsi_eta"].Fill(genJPsi.Eta(), lheWeight)
                self.hists["genJPsi_phi"].Fill(genJPsi.Phi(), lheWeight)
                self.hists["genJPsi_mass"].Fill(genJPsi.M(), lheWeight)
                self.hists["genPhoton_pt"].Fill(genPhoton.Pt(), lheWeight)
                self.hists["genPhoton_eta"].Fill(genPhoton.Eta(), lheWeight)
                self.hists["genMuon1_pt"].Fill(max(genMuon1.Pt(), genMuon2.Pt()), lheWeight)
                self.hists["genMuon1_eta"].Fill(genMuon1.Eta(), lheWeight)
                self.hists["genMuon2_pt"].Fill(min(genMuon1.Pt(), genMuon2.Pt()), lheWeight)
                self.hists["genMuon2_eta"].Fill(genMuon2.Eta(), lheWeight)
                self.hists["genCosThetaStar"].Fill(genCosThetaStar, lheWeight)
                self.hists["genCosTheta1"].Fill(genCosTheta1, lheWeight)
                self.hists["genPhi1"].Fill(genPhi1, lheWeight)
                self.hists["genCosThetaStarZtoMM"].Fill(self.returnCosThetaStar(genH, genMuon1), lheWeight)
                
                # Reweight
                #topWeight = 3./4. * (1. + genCosTheta1**2) # Transverse polarization (H, Z)
                #if 'ZToJPsiG' in self.sampleName:
                #    stitchWeight = 3./2. * (1. - genCosTheta1**2) # Longitudinal polarization (Z)
                
                # Acceptance
                if abs(genPhoton.Eta()) < 2.5:
                    self.hists["Acceptance"].Fill(1, lheWeight)
                    if abs(genMuon1.Eta()) < 2.4:
                        self.hists["Acceptance"].Fill(2, lheWeight)
                        if abs(genMuon2.Eta()) < 2.4:
                            self.hists["Acceptance"].Fill(3, lheWeight)
                            self.hists["Cutflow"].Fill(1, lheWeight)
                            if genPhoton.Pt() > 10. and genMuon1.Pt() > 5. and genMuon2.Pt() > 5.:
                                self.hists["Acceptance"].Fill(4, lheWeight)
                                isGenAcceptance = True
                self.hists["Matching"].Fill(1, lheWeight)
                if isGenAcceptance: self.hists["Matching"].Fill(2, lheWeight)
#            else:
#                print "-"*60
#                print genHIdx, genJPsiIdx, genPhotonIdx, genMuon1Idx, genMuon2Idx
#                for i in range(event.nGenPart):
#                    print i, "\t", event.GenPart_pdgId[i], "\t", event.GenPart_status[i], "\t", (event.GenPart_pdgId[event.GenPart_genPartIdxMother[i]] if event.GenPart_genPartIdxMother[i] >= 0 else "-")
        #for i in range(event.nMuon):
        #    print i, event.Muon_pt[i], event.Muon_corrected_pt[i], event.Muon_correctedUp_pt[i], event.Muon_correctedDown_pt[i]
        
        # Muons
        m1, m2 = -1, -1
        if self.muonCriteria == "PT":
            for i in range(event.nMuon):
                if event.Muon_pt[i] > self.thMuons and abs(event.Muon_eta[i]) < 2.4 and event.Muon_looseId[i]:
                    if m1 < 0: m1 = i
                    if m2 < 0 and m1 >= 0 and event.Muon_charge[m1] != event.Muon_charge[i]: m2 = i
        
        elif self.muonCriteria == "MASS":
            minMinv = 1.e6
            for i in range(event.nMuon):
                for j in range(event.nMuon):
                    if i == j or event.Muon_charge[i] == event.Muon_charge[j] or event.Muon_pt[i] < self.thMuons or event.Muon_pt[j] < self.thMuons or abs(event.Muon_eta[i]) > 2.4  or abs(event.Muon_eta[j]) > 2.4 or not event.Muon_looseId[i] or not event.Muon_looseId[j]: continue
                    tmuon1, tmuon2 = TLorentzVector(), TLorentzVector()
                    tmuon1.SetPtEtaPhiM(event.Muon_pt[i], event.Muon_eta[i], event.Muon_phi[i], event.Muon_mass[i])
                    tmuon2.SetPtEtaPhiM(event.Muon_pt[j], event.Muon_eta[j], event.Muon_phi[j], event.Muon_mass[j])
                    Minv = (tmuon1 + tmuon2).M()
                    if Minv < minMinv:
                        minMinv, m1, m2 = Minv, i, j
        
        elif self.muonCriteria == "DR":
            minDR = 1.e6
            for i in range(event.nMuon):
                for j in range(event.nMuon):
                    if i == j or event.Muon_charge[i] == event.Muon_charge[j] or event.Muon_pt[i] < self.thMuons or event.Muon_pt[j] < self.thMuons or abs(event.Muon_eta[i]) > 2.4  or abs(event.Muon_eta[j]) > 2.4 or not event.Muon_looseId[i] or not event.Muon_looseId[j]: continue
                    tmuon1, tmuon2 = TLorentzVector(), TLorentzVector()
                    tmuon1.SetPtEtaPhiM(event.Muon_pt[i], event.Muon_eta[i], event.Muon_phi[i], event.Muon_mass[i])
                    tmuon2.SetPtEtaPhiM(event.Muon_pt[j], event.Muon_eta[j], event.Muon_phi[j], event.Muon_mass[j])
                    DR = tmuon1.DeltaR(tmuon2)
                    if DR < minDR:
                        minDR, m1, m2 = DR, i, j
        
        else:
            if self.verbose >= 0: print "- Select a valid muon selection criterium"
            return False
        
        if m1 < 0 or m2 < 0:
            if self.verbose >= 2: print "- No OS loose muons in acceptance"
            return False
            
        self.hists["Cutflow"].Fill(2, lheWeight)
        
        _mu1pt, _mu2pt, _mu1ptUp, _mu2ptUp, _mu1ptDown, _mu2ptDown = event.Muon_pt[m1], event.Muon_pt[m2], event.Muon_pt[m1], event.Muon_pt[m2], event.Muon_pt[m1], event.Muon_pt[m2]
        if hasattr(event, "Muon_corrected_pt"): _mu1pt, _mu2pt, _mu1ptUp, _mu2ptUp, _mu1ptDown, _mu2ptDown = event.Muon_corrected_pt[m1], event.Muon_corrected_pt[m2], event.Muon_correctedUp_pt[m1], event.Muon_correctedUp_pt[m2], event.Muon_correctedDown_pt[m1], event.Muon_correctedDown_pt[m2]
        
        muon1, muon2, muon1NoCorr, muon2NoCorr, muon1ScaleUp, muon1ScaleDown, muon2ScaleUp, muon2ScaleDown = TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector()
        muon1.SetPtEtaPhiM(_mu1pt, event.Muon_eta[m1], event.Muon_phi[m1], event.Muon_mass[m1])
        muon2.SetPtEtaPhiM(_mu2pt, event.Muon_eta[m2], event.Muon_phi[m2], event.Muon_mass[m2])
        muon1NoCorr.SetPtEtaPhiM(event.Muon_pt[m1], event.Muon_eta[m1], event.Muon_phi[m1], event.Muon_mass[m1])
        muon2NoCorr.SetPtEtaPhiM(event.Muon_pt[m2], event.Muon_eta[m2], event.Muon_phi[m2], event.Muon_mass[m2])
        muon1ScaleUp.SetPtEtaPhiM(_mu1ptUp, event.Muon_eta[m1], event.Muon_phi[m1], event.Muon_mass[m1])
        muon1ScaleDown.SetPtEtaPhiM(_mu1ptDown, event.Muon_eta[m1], event.Muon_phi[m1], event.Muon_mass[m1])
        muon2ScaleUp.SetPtEtaPhiM(_mu2ptUp, event.Muon_eta[m2], event.Muon_phi[m2], event.Muon_mass[m2])
        muon2ScaleDown.SetPtEtaPhiM(_mu2ptDown, event.Muon_eta[m2], event.Muon_phi[m2], event.Muon_mass[m2])
        muonP, muonM = (muon1.Clone(), muon2.Clone()) if event.Muon_charge[m1] > event.Muon_charge[m2] else (muon2.Clone(), muon1.Clone())
        muon1v2, muon2v2 = TVector2(muon1.Px(), muon1.Py()), TVector2(muon2.Px(), muon2.Py())
        
        # Gen matching
        if self.isMC and self.isSignal and isGenMatched:
            self.hists["Matching"].Fill(3, lheWeight)
            if (genMuon1.DeltaR(muon1) < self.thMatching and genMuon2.DeltaR(muon2) < self.thMatching) or (genMuon1.DeltaR(muon2) < self.thMatching and genMuon2.DeltaR(muon1) < self.thMatching):
                genMuonsMatching = True
                self.hists["Matching"].Fill(4, lheWeight)
        
        jpsi = muon1 + muon2
        jpsi_muScaleUp = muon1ScaleUp + muon2ScaleUp
        jpsi_muScaleDown = muon1ScaleDown + muon2ScaleDown
        
        if jpsi.M() < 2. or jpsi.M() > 12.:
            if self.verbose >= 2: print "- Dimuon invariant mass < 2 or > 12 GeV"
            return False
        
        self.hists["Cutflow"].Fill(3, lheWeight)
        
        # Photons
        p0 = -1
        for i in range(event.nPhoton):
            if event.Photon_pt[i] > self.thPhoton and abs(event.Photon_eta[i]) < 2.5 and event.Photon_pfRelIso03_all[i] < 0.5: # and event.Photon_mvaID_WP80[i]:
                if p0 < 0: p0 = i
        
        if p0 < 0:
            if self.verbose >= 2: print "- No isolated photons in acceptance"
            return False
        
        self.hists["Cutflow"].Fill(4, lheWeight)
        
        _phScaleUp = + event.Photon_energyErr[p0]
        _phScaleDown = - event.Photon_energyErr[p0]
        _phResUp = + event.Photon_dEsigmaUp[p0] if hasattr(event, "Photon_dEsigmaUp") else 0.
        _phResDown = + event.Photon_dEsigmaDown[p0] if hasattr(event, "Photon_dEsigmaDown") else 0.
        
        photon, photonScaleUp, photonScaleDown, photonResUp, photonResDown = TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector()
        photon.SetPtEtaPhiM(event.Photon_pt[p0], event.Photon_eta[p0], event.Photon_phi[p0], event.Photon_mass[p0])
        photonScaleUp.SetPtEtaPhiM(event.Photon_pt[p0] + _phScaleUp, event.Photon_eta[p0], event.Photon_phi[p0], event.Photon_mass[p0]) #event.Photon_dEscaleUp[p0]
        photonScaleDown.SetPtEtaPhiM(event.Photon_pt[p0] + _phScaleDown, event.Photon_eta[p0], event.Photon_phi[p0], event.Photon_mass[p0]) #event.Photon_dEscaleDown[p0]
        photonResUp.SetPtEtaPhiM(event.Photon_pt[p0] + _phResUp, event.Photon_eta[p0], event.Photon_phi[p0], event.Photon_mass[p0])
        photonResDown.SetPtEtaPhiM(event.Photon_pt[p0] + _phResDown, event.Photon_eta[p0], event.Photon_phi[p0], event.Photon_mass[p0])
        photonv2 = TVector2(photon.Px(), photon.Py())
        
        # Gen matching
        if self.isMC and self.isSignal and isGenMatched:
            self.hists["Matching"].Fill(5, lheWeight)
            if genPhoton.DeltaR(photon) < (self.thMatching*5.):
                genPhotonMatching = True
                self.hists["Matching"].Fill(6, lheWeight)
            genEventMatching = genMuonsMatching and genPhotonMatching
            if genEventMatching: self.hists["Matching"].Fill(7, lheWeight)
        
        # MET
        met, metPlusPhoton = TVector2(), TVector2()
        met.SetMagPhi(event.MET_pt, event.MET_phi)
        metPlusPhoton.Set(met.Px() + photon.Px(), met.Py() + photon.Py())
        
        h = jpsi + photon
        
        h_noCorr = muon1NoCorr + muon2NoCorr + photon
        h_muScaleUp = jpsi_muScaleUp + photon
        h_muScaleDown = jpsi_muScaleDown + photon
        h_phScaleUp = jpsi + photonScaleUp
        h_phScaleDown = jpsi + photonScaleDown
        h_phResUp = jpsi + photonResUp
        h_phResDown = jpsi + photonResDown

        jpsi_pt = jpsi.Pt()
        jpsi_eta = jpsi.Eta()
        jpsi_phi = jpsi.Phi()
        jpsi_mass = jpsi.M()
        jpsi_dEta = abs(muon1.Eta() - muon2.Eta())
        jpsi_dPhi = abs(muon1.DeltaPhi(muon2))
        jpsi_dR = muon1.DeltaR(muon2)
        
        h_pt = h.Pt()
        h_eta = h.Eta()
        h_phi = h.Phi()
        h_mass = h.M()
        h_dEta = abs(jpsi.Eta() - photon.Eta())
        h_dPhi = abs(jpsi.DeltaPhi(photon))
        h_dR = jpsi.DeltaR(photon)
        
        minMuonPfIso = min(event.Muon_pfRelIso04_all[m1], event.Muon_pfRelIso04_all[m2])
        maxMuonPfIso = max(event.Muon_pfRelIso04_all[m1], event.Muon_pfRelIso04_all[m2])
        
        Muon1TrkIso, Muon2TrkIso = event.Muon_tkRelIso[m1], event.Muon_tkRelIso[m2]
        if jpsi_dR < 0.3:
            Muon1TrkIso = max(0., (Muon1TrkIso * event.Muon_pt[m1] * event.Muon_tunepRelPt[m1] - event.Muon_pt[m2] * event.Muon_tunepRelPt[m2])) / (event.Muon_pt[m1] * event.Muon_tunepRelPt[m1])
            Muon2TrkIso = max(0., (Muon2TrkIso * event.Muon_pt[m2] * event.Muon_tunepRelPt[m2] - event.Muon_pt[m1] * event.Muon_tunepRelPt[m1])) / (event.Muon_pt[m2] * event.Muon_tunepRelPt[m2])
        minMuonTrkSubIso = min(Muon1TrkIso, Muon2TrkIso)
        maxMuonTrkSubIso = max(Muon1TrkIso, Muon2TrkIso)
        
        Muon1PfSubIso, Muon2PfSubIso = event.Muon_pfRelIso04_all[m1], event.Muon_pfRelIso04_all[m2]
        if jpsi_dR < 0.4:
            Muon1PfSubIso = max(0., (Muon1PfSubIso * event.Muon_pt[m1] - event.Muon_pt[m2])) / event.Muon_pt[m1]
            Muon2PfSubIso = max(0., (Muon2PfSubIso * event.Muon_pt[m2] - event.Muon_pt[m1])) / event.Muon_pt[m2]
        minMuonPfSubIso = min(Muon1PfSubIso, Muon2PfSubIso)
        maxMuonPfSubIso = max(Muon1PfSubIso, Muon2PfSubIso)
        
        minMuonMetDPhi = min(abs(muon1v2.DeltaPhi(met)), abs(muon2v2.DeltaPhi(met)))
        maxMuonMetDPhi = max(abs(muon1v2.DeltaPhi(met)), abs(muon2v2.DeltaPhi(met)))
        photonMetDPhi = abs(photonv2.DeltaPhi(met))
        metPlusPhotonDPhi = abs(met.DeltaPhi(metPlusPhoton))
        
        centrality = self.returnCentrality(h, muonP, muonM, photon)
        trCentr = min(math.log(1./(1.-centrality))/10., 1.)
        
        cosThetaStar = self.returnCosThetaStar(h, jpsi)
        cosTheta1 = self.returnCosTheta1(jpsi, muonP, muonM, photon)
        phi1 = self.returnPhi1(h, muonP, muonM)

        # Weights
#        if self.isMC:
#            triggerWeight = self.muSFs.getTriggerSF(event.Muon_pt[m1], event.Muon_eta[m1])
#            IdSF1 = self.muSFs.getIdSF(event.Muon_pt[0], event.Muon_eta[0], 2)
#            IsoSF1 = self.muSFs.getIsoSF(event.Muon_pt[0], event.Muon_eta[0])
#            IdSF2 = self.muSFs.getIdSF(event.Muon_pt[1], event.Muon_eta[1], 2)
#            IsoSF2 = self.muSFs.getIsoSF(event.Muon_pt[1], event.Muon_eta[1])
#            IdIsoSF3 = self.elSFs.getIdIsoSF(event.Electron_pt[0], event.Electron_eta[0])
#            leptonWeight = IdSF1 * IsoSF1 * IdSF2 * IsoSF2 * IdIsoSF3
        
        passedMETFilters = True
        filters = ["Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_BadPFMuonFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_ecalBadCalibFilter", "Flag_ecalBadCalibFilterV2"]
        if not self.isMC: filters += ["Flag_eeBadScFilter"]
        for f in filters:
            if hasattr(event, f) and getattr(event, f) == False: passedMETFilters = False
#        try:
##            if event.Flag_goodVertices: print "Flag_goodVertices"
##            if event.Flag_globalSuperTightHalo2016Filter: print "Flag_globalSuperTightHalo2016Filter"
##            if event.Flag_BadPFMuonFilter: print "Flag_BadPFMuonFilter"
##            if event.Flag_EcalDeadCellTriggerPrimitiveFilter: print "Flag_EcalDeadCellTriggerPrimitiveFilter"
##            if event.Flag_HBHENoiseFilter: print "Flag_HBHENoiseFilter"
##            if event.Flag_HBHENoiseIsoFilter: print "Flag_HBHENoiseIsoFilter"
###            if (self.isMC or event.Flag_eeBadScFilter): print "Flag_eeBadScFilter"
##            if event.Flag_ecalBadCalibFilter: print "Flag_ecalBadCalibFilterV2"
#            if event.Flag_goodVertices and event.Flag_globalSuperTightHalo2016Filter and event.Flag_BadPFMuonFilter and event.Flag_EcalDeadCellTriggerPrimitiveFilter and event.Flag_HBHENoiseFilter and event.Flag_HBHENoiseIsoFilter: # and event.Flag_ecalBadCalibFilter: #and (self.isMC or event.Flag_eeBadScFilter): FIXME
#                passedMETFilters = True
##            if not self.isMC:
##                if not event.Flag_eeBadScFilter:
##                    passedMETFilters = False
#        except:
#            passedMETFilters = False
        
        
        
        ### Event variables ###
        tLV, tJ1, tJ2 = TLorentzVector(), TLorentzVector(), TLorentzVector()
        
        # Muons
        for i in range(event.nMuon):
            if i != m1 and i != m2 and event.Muon_pt[i] > 10. and abs(event.Muon_eta[i]) < 2.4 and event.Muon_mediumId[i] and event.Muon_pfRelIso04_all[i] < 0.15:
                tLV.SetPtEtaPhiM(event.Muon_pt[i], event.Muon_eta[i], event.Muon_phi[i], event.Muon_mass[i])
                if muon1.DeltaR(tLV) > 0.4 and muon2.DeltaR(tLV) > 0.4 and photon.DeltaR(tLV) > 0.3:
                    nCleanMuon += 1
        
        # Electrons
        for i in range(event.nElectron):
            if event.Electron_pt[i] > 10. and abs(event.Electron_eta[i]) < 2.5 and event.Electron_cutBased[i] >= 2:
                tLV.SetPtEtaPhiM(event.Electron_pt[i], event.Electron_eta[i], event.Electron_phi[i], event.Electron_mass[i])
                if muon1.DeltaR(tLV) > 0.4 and muon2.DeltaR(tLV) > 0.4 and photon.DeltaR(tLV) > 0.3:
                    nCleanElectron += 1
        
        # Taus
        for i in range(event.nTau):
            if event.Tau_pt[i] > 20. and abs(event.Tau_eta[i]) < 2.5 and event.Tau_idDeepTau2017v2p1VSe[i] >= 16 and event.Tau_idDeepTau2017v2p1VSmu[i] >= 8 and event.Tau_idDeepTau2017v2p1VSjet[i] >= 16 and event.Tau_rawIsodR03[i] < 0.15:
                tLV.SetPtEtaPhiM(event.Tau_pt[i], event.Tau_eta[i], event.Tau_phi[i], event.Tau_mass[i])
                if muon1.DeltaR(tLV) > 0.4 and muon2.DeltaR(tLV) > 0.4 and photon.DeltaR(tLV) > 0.3:
                    nCleanTau += 1
        
        # Photons
        for i in range(event.nPhoton):
            if i != p0 and event.Photon_pt[i] > 15. and abs(event.Photon_eta[i]) < 2.5 and event.Photon_pfRelIso03_all[i] < 0.15 and event.Photon_mvaID_WP90[i]:
                tLV.SetPtEtaPhiM(event.Photon_pt[i], event.Photon_eta[i], event.Photon_phi[i], event.Photon_mass[i])
                if muon1.DeltaR(tLV) > 0.4 and muon2.DeltaR(tLV) > 0.4 and photon.DeltaR(tLV) > 0.3:
                    nCleanPhoton += 1
        
        # Jets and Event variables
        idxJets = []
        nJet = event.nJet if hasattr(event, "nJet") else 0 # should not happen, but it does
        for i in range(nJet):
            if event.Jet_jetId[i] >= 6 and event.Jet_pt[i] > 30. and abs(event.Jet_eta[i]) < 5.:
                tLV.SetPtEtaPhiM(event.Jet_pt[i], event.Jet_eta[i], event.Jet_phi[i], event.Jet_mass[i])
                if muon1.DeltaR(tLV) > 0.4 and muon2.DeltaR(tLV) > 0.4 and photon.DeltaR(tLV) > 0.4:
                    idxJets.append(i)
                    if abs(event.Jet_eta[i]) < 2.5:
                        HT30 += event.Jet_pt[i]
                        nCleanCentralJet += 1
                        if event.Jet_btagDeepB[i] >= self.btagMedium: nCleanBTagJet += 1
                    else:
                        nCleanForwardJet += 1
        
        # max invariant mass of all jets
        vbfmjj, vbfdR = -1., -1.
        for i in range(len(idxJets)):
            for j in range(len(idxJets)):
                if i != j:
                    tJ1.SetPtEtaPhiM(event.Jet_pt[i], event.Jet_eta[i], event.Jet_phi[i], event.Jet_mass[i])
                    tJ2.SetPtEtaPhiM(event.Jet_pt[j], event.Jet_eta[j], event.Jet_phi[j], event.Jet_mass[j])
                    tmjj = (tJ1 + tJ2).M()
                    if tmjj > vbfmjj:
                        vbfmjj = tmjj
                        vbfdR = tJ1.DeltaR(tJ2)
        
        
#        # Categorization
#        isINC, isVBF, isHFL, isLEP = False, False, False, False
#        if nCleanMuon + nCleanElectron >= 1: isLEP = True
#        elif nCleanBTagJet >= 1: isHFL = True
#        elif vbfmjj > 350.: isVBF = True
#        else: isINC = True
        
        
        
        if self.isMC:
            # muon Id SF and errors
            # Muon1_mediumId_SF          = self.muSFdict["loose"][self.year].getIdSF     (event.Muon_pt[m1], event.Muon_eta[m1])
            # Muon1_mediumId_SFError     = self.muSFdict["loose"][self.year].getIdSFerror(event.Muon_pt[m1], event.Muon_eta[m1])
            # Muon2_mediumId_SF          = self.muSFdict["loose"][self.year].getIdSF     (event.Muon_pt[m2], event.Muon_eta[m2])
            # Muon2_mediumId_SFError     = self.muSFdict["loose"][self.year].getIdSFerror(event.Muon_pt[m2], event.Muon_eta[m2])
            # Muon1_tightId_SF           = self.muSFdict["tight"][self.year].getIdSF     (event.Muon_pt[m1], event.Muon_eta[m1])
            # Muon1_tightId_SFError      = self.muSFdict["tight"][self.year].getIdSFerror(event.Muon_pt[m1], event.Muon_eta[m1])
            # Muon2_tightId_SF           = self.muSFdict["tight"][self.year].getIdSF     (event.Muon_pt[m2], event.Muon_eta[m2])
            # Muon2_tightId_SFError      = self.muSFdict["tight"][self.year].getIdSFerror(event.Muon_pt[m2], event.Muon_eta[m2])
            Muon1_mediumPromptId_SF      = self.muSFdict["prompt"][self.year].getIdSF     (event.Muon_pt[m1], event.Muon_eta[m1])
            Muon1_mediumPromptId_SFError = self.muSFdict["prompt"][self.year].getIdSFerror(event.Muon_pt[m1], event.Muon_eta[m1])
            Muon2_mediumPromptId_SF      = self.muSFdict["prompt"][self.year].getIdSF     (event.Muon_pt[m2], event.Muon_eta[m2])
            Muon2_mediumPromptId_SFError = self.muSFdict["prompt"][self.year].getIdSFerror(event.Muon_pt[m2], event.Muon_eta[m2])

            # muon Iso SF and errors
            # Muon1_looseIso_SF          = self.muSFdict["loose"][self.year].getIsoSF     (event.Muon_pt[m1], event.Muon_eta[m1])
            # Muon1_looseIso_SFError     = self.muSFdict["loose"][self.year].getIsoSFerror(event.Muon_pt[m1], event.Muon_eta[m1])
            # Muon2_looseIso_SF          = self.muSFdict["loose"][self.year].getIsoSF     (event.Muon_pt[m2], event.Muon_eta[m2])
            # Muon2_looseIso_SFError     = self.muSFdict["loose"][self.year].getIsoSFerror(event.Muon_pt[m2], event.Muon_eta[m2])
            Muon1_tightIso_SF          = self.muSFdict["prompt"][self.year].getIsoSF     (event.Muon_pt[m1], event.Muon_eta[m1])
            Muon1_tightIso_SFError     = self.muSFdict["prompt"][self.year].getIsoSFerror(event.Muon_pt[m1], event.Muon_eta[m1])
            Muon2_tightIso_SF          = self.muSFdict["prompt"][self.year].getIsoSF     (event.Muon_pt[m2], event.Muon_eta[m2])
            Muon2_tightIso_SFError     = self.muSFdict["prompt"][self.year].getIsoSFerror(event.Muon_pt[m2], event.Muon_eta[m2])

            # photon Id+Iso SF and errors
            Photon1_mvaID_WP80_SF      = self.phoSFdict["MVA80"][self.year].getIdIsoSF     (event.Photon_pt[p0], event.Photon_eta[p0])
            Photon1_mvaID_WP80_SFError = self.phoSFdict["MVA80"][self.year].getIdIsoSFerror(event.Photon_pt[p0], event.Photon_eta[p0])
            Photon1_mvaID_WP90_SF      = self.phoSFdict["MVA90"][self.year].getIdIsoSF     (event.Photon_pt[p0], event.Photon_eta[p0])
            Photon1_mvaID_WP90_SFError = self.phoSFdict["MVA90"][self.year].getIdIsoSFerror(event.Photon_pt[p0], event.Photon_eta[p0])

            # pileup SF
            puWeight            = self.puWeightSFdict["central"][self.year].getWeight(event.Pileup_nTrueInt)
            puWeightUp            = self.puWeightSFdict["up"     ][self.year].getWeight(event.Pileup_nTrueInt)
            puWeightDown            = self.puWeightSFdict["down"   ][self.year].getWeight(event.Pileup_nTrueInt)

            # trigger SF
            triggerWeight   = 1.0
            triggerWeightUp   = 1.0
            triggerWeightDown   = 1.0
            
            triggerMuonEGWeight   = self.triggerSFdict["SingleMuonPhoton"][self.year]
            triggerMuonEGWeightUp   = self.triggerSFdict["SingleMuonPhotonUp"][self.year]
            triggerMuonEGWeightDown   = self.triggerSFdict["SingleMuonPhotonDown"][self.year]

            triggerJPsiWeight   = self.triggerSFdict["JPsi"][self.year]
            triggerJPsiWeightUp   = self.triggerSFdict["JPsiUp"][self.year]
            triggerJPsiWeightDown   = self.triggerSFdict["JPsiDown"][self.year]
            
            if self.year == 2017:
                triggerMuonEGWeight *= 1. - 14400./self.lumi
                triggerMuonEGWeightUp *= 1. - 14400./self.lumi
                triggerMuonEGWeightDown *= 1. - 14400./self.lumi

            # group muon SF
            leptonWeight        = Muon1_mediumPromptId_SF * Muon2_mediumPromptId_SF * Muon1_tightIso_SF * Muon2_tightIso_SF
            leptonWeightUp        = (Muon1_mediumPromptId_SF+Muon1_mediumPromptId_SFError) * \
                                             (Muon2_mediumPromptId_SF+Muon2_mediumPromptId_SFError) * \
                                             (Muon1_tightIso_SF      +Muon1_tightIso_SFError      ) * \
                                             (Muon2_tightIso_SF      +Muon2_tightIso_SFError      )
            leptonWeightDown        = (Muon1_mediumPromptId_SF-Muon1_mediumPromptId_SFError) * \
                                             (Muon2_mediumPromptId_SF-Muon2_mediumPromptId_SFError) * \
                                             (Muon1_tightIso_SF      -Muon1_tightIso_SFError      ) * \
                                             (Muon2_tightIso_SF      -Muon2_tightIso_SFError      )
            # group photon SF
            photonWeight       = Photon1_mvaID_WP80_SF
            photonWeightUp       = Photon1_mvaID_WP80_SF + Photon1_mvaID_WP80_SFError
            photonWeightDown       = Photon1_mvaID_WP80_SF - Photon1_mvaID_WP80_SFError
            
            # global event weight
            eventWeightLumi = self.lumiWeight * lheWeight * puWeight * topWeight * qcdnloWeight * qcdnnloWeight * ewknloWeight * triggerWeight * leptonWeight
            eventMuonEGWeightLumi = self.lumiWeight * lheWeight * puWeight * topWeight * qcdnloWeight * qcdnnloWeight * ewknloWeight * triggerMuonEGWeight * leptonWeight
            eventJPsiWeightLumi = self.lumiWeight * lheWeight * puWeight * topWeight * qcdnloWeight * qcdnnloWeight * ewknloWeight * triggerJPsiWeight * leptonWeight
            
            
        
        self.out.fillBranch("isMC", self.isMC)
        self.out.fillBranch("is2016", (self.year == 2016))
        self.out.fillBranch("is2017", (self.year == 2017))
        self.out.fillBranch("is2018", (self.year == 2018))
        self.out.fillBranch("isSingleMuonTrigger", isSingleMuonTrigger)
        self.out.fillBranch("isSingleMuonPhotonTrigger", isSingleMuonPhotonTrigger)
        self.out.fillBranch("isSingleMuonNoFiltersPhotonTrigger", isSingleMuonNoFiltersPhotonTrigger)
        self.out.fillBranch("isDoubleMuonTrigger", isDoubleMuonTrigger)
        self.out.fillBranch("isDoubleMuonPhotonTrigger", isDoubleMuonPhotonTrigger)
        self.out.fillBranch("isJPsiTrigger", isJPsiTrigger)
        self.out.fillBranch("passedMETFilters", passedMETFilters)
        self.out.fillBranch("nCleanElectron", nCleanElectron)
        self.out.fillBranch("nCleanMuon", nCleanMuon)
        self.out.fillBranch("nCleanTau", nCleanTau)
        self.out.fillBranch("nCleanPhoton", nCleanPhoton)
        self.out.fillBranch("nCleanCentralJet", nCleanCentralJet)
        self.out.fillBranch("nCleanForwardJet", nCleanForwardJet)
        self.out.fillBranch("nCleanBTagJet", nCleanBTagJet)
        self.out.fillBranch("HT30", HT30)
        self.out.fillBranch("iPhoton", p0)
        self.out.fillBranch("iMuon1", m1)
        self.out.fillBranch("iMuon2", m2)
        #
#        self.out.fillBranch("isINC", isINC)
#        self.out.fillBranch("isVBF", isVBF)
#        self.out.fillBranch("isHFL", isHFL)
#        self.out.fillBranch("isLEP", isLEP)
        #
        self.out.fillBranch("Muon1_pt", muon1.Pt())
        self.out.fillBranch("Muon1_eta", event.Muon_eta[m1])
        self.out.fillBranch("Muon2_pt", muon2.Pt())
        self.out.fillBranch("Muon2_eta", event.Muon_eta[m2])
        self.out.fillBranch("Muon1_pfRelIso03", event.Muon_pfRelIso03_all[m1])
        self.out.fillBranch("Muon2_pfRelIso03", event.Muon_pfRelIso03_all[m2])
        self.out.fillBranch("Muon1_pfRelIso04", event.Muon_pfRelIso04_all[m1])
        self.out.fillBranch("Muon2_pfRelIso04", event.Muon_pfRelIso04_all[m2])
        self.out.fillBranch("Muon1_mediumId", event.Muon_mediumId[m1])
        self.out.fillBranch("Muon2_mediumId", event.Muon_mediumId[m2])
        self.out.fillBranch("Muon1_mediumPromptId", event.Muon_mediumPromptId[m1])
        self.out.fillBranch("Muon2_mediumPromptId", event.Muon_mediumPromptId[m2])
        self.out.fillBranch("Muon1_ip3d", event.Muon_ip3d[m1])
        self.out.fillBranch("Muon2_ip3d", event.Muon_ip3d[m2])
        self.out.fillBranch("minMuonPfIso", minMuonPfIso)
        self.out.fillBranch("maxMuonPfIso", maxMuonPfIso)
        self.out.fillBranch("minMuonPfSubIso", minMuonPfSubIso)
        self.out.fillBranch("maxMuonPfSubIso", maxMuonPfSubIso)
        self.out.fillBranch("minMuonTrkSubIso", minMuonTrkSubIso)
        self.out.fillBranch("maxMuonTrkSubIso", maxMuonTrkSubIso)
        self.out.fillBranch("Muon12_diffdxy", abs(event.Muon_dxy[m1]-event.Muon_dxy[m2]))
        self.out.fillBranch("Muon12_diffdz", abs(event.Muon_dz[m1]-event.Muon_dz[m2]))
        self.out.fillBranch("Muon12_signdxy", abs(event.Muon_dxy[m1]-event.Muon_dxy[m2]) / math.sqrt(event.Muon_dxyErr[m1]**2 + event.Muon_dxyErr[m2]**2))
        self.out.fillBranch("Muon12_signdz", abs(event.Muon_dz[m1]-event.Muon_dz[m2]) / math.sqrt(event.Muon_dzErr[m1]**2 + event.Muon_dzErr[m2]**2))
        self.out.fillBranch("Photon1_pt", event.Photon_pt[p0])
        self.out.fillBranch("Photon1_eta", event.Photon_eta[p0])
        self.out.fillBranch("Photon1_mvaID_WP80", event.Photon_mvaID_WP80[p0])
        self.out.fillBranch("Photon1_mvaID_WP90", event.Photon_mvaID_WP90[p0])
        self.out.fillBranch("Photon1_r9", event.Photon_r9[p0])
        self.out.fillBranch("Photon1_pfRelIso03", event.Photon_pfRelIso03_all[p0])
        self.out.fillBranch("Photon1_pixelSeed", event.Photon_pixelSeed[p0])
        #
        self.out.fillBranch("JPsi_pt", jpsi_pt)
        self.out.fillBranch("JPsi_eta", jpsi_eta)
        self.out.fillBranch("JPsi_phi", jpsi_phi)
        self.out.fillBranch("JPsi_mass", jpsi_mass)
        self.out.fillBranch("JPsi_dEta", jpsi_dEta)
        self.out.fillBranch("JPsi_dPhi", jpsi_dPhi)
        self.out.fillBranch("JPsi_dR", jpsi_dR)
        self.out.fillBranch("H_pt", h_pt)
        self.out.fillBranch("H_eta", h_eta)
        self.out.fillBranch("H_phi", h_phi)
        self.out.fillBranch("H_mass", h_mass)
        self.out.fillBranch("H_dEta", h_dEta)
        self.out.fillBranch("H_dPhi", h_dPhi)
        self.out.fillBranch("H_dR", h_dR)
        #
        self.out.fillBranch("JPsi_mass_muScaleUp", jpsi_muScaleUp.M())
        self.out.fillBranch("JPsi_mass_muScaleDown", jpsi_muScaleDown.M())
        self.out.fillBranch("H_mass_muScaleUp", h_muScaleUp.M())
        self.out.fillBranch("H_mass_muScaleDown", h_muScaleDown.M())
        #self.out.fillBranch("H_mass_phScaleUp", h_phScaleUp.M())
        #self.out.fillBranch("H_mass_phScaleDown", h_phScaleDown.M())
        self.out.fillBranch("H_mass_phResUp", h_phResUp.M())
        self.out.fillBranch("H_mass_phResDown", h_phResDown.M())
        self.out.fillBranch("H_mass_noCorr", h_noCorr.M())
        #
        self.out.fillBranch("minMuonMetDPhi", minMuonMetDPhi)
        self.out.fillBranch("maxMuonMetDPhi", maxMuonMetDPhi)
        self.out.fillBranch("photonMetDPhi", photonMetDPhi)
        self.out.fillBranch("metPlusPhotonDPhi", metPlusPhotonDPhi)
        self.out.fillBranch("vbfMjj", vbfmjj)
        self.out.fillBranch("vbfDR", vbfdR)
        self.out.fillBranch("centrality", centrality)
        self.out.fillBranch("trCentr", trCentr)
        self.out.fillBranch("cosThetaStar", cosThetaStar)
        self.out.fillBranch("cosTheta1", cosTheta1)
        self.out.fillBranch("phi1", phi1)
        self.out.fillBranch("genCosThetaStar", genCosThetaStar)
        self.out.fillBranch("genCosTheta1", genCosTheta1)
        self.out.fillBranch("genPhi1", genPhi1)
        self.out.fillBranch("lumiWeight", self.lumiWeight)
        self.out.fillBranch("lheWeight", lheWeight)
        self.out.fillBranch("stitchWeight", stitchWeight)
        self.out.fillBranch("topWeight", topWeight)
        self.out.fillBranch("qcdnloWeight", qcdnloWeight)
        self.out.fillBranch("qcdnnloWeight", qcdnnloWeight)
        self.out.fillBranch("ewknloWeight", ewknloWeight)
        self.out.fillBranch("puWeight", puWeight)
        self.out.fillBranch("puWeightUp", puWeightUp)
        self.out.fillBranch("puWeightDown", puWeightDown)
        self.out.fillBranch("triggerWeight", triggerWeight)
        self.out.fillBranch("triggerWeightUp", triggerWeightUp)
        self.out.fillBranch("triggerWeightDown", triggerWeightDown)
        self.out.fillBranch("triggerMuonEGWeight", triggerMuonEGWeight)
        self.out.fillBranch("triggerMuonEGWeightUp", triggerMuonEGWeightUp)
        self.out.fillBranch("triggerMuonEGWeightDown", triggerMuonEGWeightDown)
        self.out.fillBranch("triggerJPsiWeight", triggerJPsiWeight)
        self.out.fillBranch("triggerJPsiWeightUp", triggerJPsiWeightUp)
        self.out.fillBranch("triggerJPsiWeightDown", triggerJPsiWeightDown)
        self.out.fillBranch("leptonWeight", leptonWeight)
        self.out.fillBranch("leptonWeightUp", leptonWeightUp)
        self.out.fillBranch("leptonWeightDown", leptonWeightDown)
        self.out.fillBranch("photonWeight", photonWeight)
        self.out.fillBranch("photonWeightUp", photonWeightUp)
        self.out.fillBranch("photonWeightDown", photonWeightDown)
        self.out.fillBranch("eventMuonEGWeightLumi", eventMuonEGWeightLumi)
        self.out.fillBranch("eventJPsiWeightLumi", eventJPsiWeightLumi)
        self.out.fillBranch("eventWeightLumi", eventWeightLumi)
        
        if self.verbose >= 2: print "+ Tree filled"
        self.N += 1
        
        # Uncertainties
        self.hists["central"].Fill(h.M(), 1.)
        self.hists["phScaleUp"].Fill(h_phScaleUp.M(), 1.)
        self.hists["phScaleDown"].Fill(h_phScaleDown.M(), 1.)
        self.hists["phResUp"].Fill(h_phResUp.M(), 1.)
        self.hists["phResDown"].Fill(h_phResDown.M(), 1.)
        
        # Cutflows
        
        self.hists["Cutflow"].Fill(5, lheWeight)
        
        if event.Muon_mediumPromptId[m1] and event.Muon_mediumPromptId[m2]:
            self.hists["Cutflow"].Fill(6, lheWeight)
            if minMuonPfIso < 0.15 and maxMuonPfIso < 0.15:
                self.hists["Cutflow"].Fill(7, lheWeight)
                if event.Photon_mvaID_WP80[p0] > 0:
                    self.hists["Cutflow"].Fill(8, lheWeight)
                    if event.Photon_pixelSeed[p0] == 0:
                        self.hists["Cutflow"].Fill(9, lheWeight)
                        # signal-dependent cuts and triggers
                        if self.isZSignal:
                            if event.Photon_pt[p0] > 15. and jpsi_pt > 26.:
                                self.hists["Cutflow"].Fill(10, lheWeight)
                                if isJPsiTrigger:
                                     self.hists["Cutflow"].Fill(11, lheWeight)
                                     if (jpsi_mass> 3.0 and jpsi_mass < 3.2) or (jpsi_mass > 3.6 and jpsi_mass < 3.8):
                                        self.hists["Cutflow"].Fill(12, lheWeight)
                        elif self.isHSignal:
                            if event.Photon_pt[p0] > 32. and muon1.Pt() > 18.:
                                self.hists["Cutflow"].Fill(10, lheWeight)
                                if isSingleMuonPhotonTrigger:
                                     self.hists["Cutflow"].Fill(11, lheWeight)
                                     if (jpsi_mass> 3.0 and jpsi_mass < 3.2) or (jpsi_mass > 3.6 and jpsi_mass < 3.8):
                                        self.hists["Cutflow"].Fill(12, lheWeight)
        
        return True





    def returnCentrality(self, theH, theL1, theL2, thePhoton):
        h, l1, l2, g = TLorentzVector(theH), TLorentzVector(theL1), TLorentzVector(theL2), TLorentzVector(thePhoton)
        # Boost objects to the H rest frame
        l1.Boost( -h.BoostVector() )
        l2.Boost( -h.BoostVector() )
        g.Boost( -h.BoostVector() )
        value = (l1.Pt() + l2.Pt() + g.Pt())/(l1.Energy() + l2.Energy() + g.Energy())
        if value != value or math.isinf(value): return -1.
        return value
    
    def returnCosThetaStar(self, theH, theJPsi):
        h, j = TLorentzVector(theH), TLorentzVector(theJPsi)
        # Boost the Z to the A rest frame
        j.Boost( -h.BoostVector() )
        value = j.CosTheta()
        if value != value or math.isinf(value): return -2.
        return value

    def returnCosTheta1(self, theJPsi, theL1, theL2, thePhoton):
        j, l1, l2, g = TLorentzVector(theJPsi), TLorentzVector(theL1), TLorentzVector(theL2), TLorentzVector(thePhoton)
        # Boost objects to the JPsi rest frame
        l1.Boost( -j.BoostVector() )
        l2.Boost( -j.BoostVector() )
        g.Boost( -j.BoostVector() )
        # cos theta = Gamma dot L1 / (|Gamma|*|L1|)
        value = g.Vect().Dot( l1.Vect() ) / ( (g.Vect().Mag())*(l1.Vect().Mag()) ) if g.Vect().Mag() > 0. and l1.Vect().Mag() > 0. else -2
        if value != value or math.isinf(value): return -2.
        return value

    def returnPhi1(self, theH, theL1, theL2):
        h, l1, l2 = TLorentzVector(theH), TLorentzVector(theL1), TLorentzVector(theL2)
        beamAxis = TVector3(0., 0., 1.)
        # Boost objects to the A rest frame
        l1.Boost( -h.BoostVector() )
        l2.Boost( -h.BoostVector() )
        # Reconstruct JPsi in H rest frame
        j = l1 + l2
        # Build unit vectors orthogonal to the decay planes
        Zplane = TVector3(l1.Vect().Cross( l2.Vect() )) # L1 x L2
        Bplane = TVector3(beamAxis.Cross( j.Vect() )) # Beam x JPsi, beam/JPsi plane
        if Zplane.Mag() == 0. or Bplane.Mag() == 0.: return -4.
        Zplane.SetMag(1.)
        Bplane.SetMag(1.)
        # Sign of Phi1
        sgn = j.Vect().Dot( Zplane.Cross(Bplane) )
        sgn /= abs(sgn)
        if abs(Zplane.Dot(Bplane)) > 1.: return -5.
        value = sgn * math.acos( Zplane.Dot(Bplane) )
        if value != value or math.isinf(value): return -5.
        return value


