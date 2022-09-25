import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

def sig_eff(fileList, muon_pt_cut = 17, photon_pt_cut = 30):
    eff = []

    for f in fileList:
        events = NanoEventsFactory.from_root(f, schemaclass = NanoAODSchema.v7, treepath='Events').events()

        den = ak.sum( ak.any(events.Muon.pt > muon_pt_cut, axis = -1) & ak.any(events.Photon.pt > photon_pt_cut, axis = -1))
        num = ak.sum( ak.any(events.Muon[events.HLT.Mu17_Photon30_CaloIdL_L1ISO].pt > muon_pt_cut, axis = -1) & ak.any(events.Photon[events.HLT.Mu17_Photon30_CaloIdL_L1ISO].pt > photon_pt_cut, axis = -1))

        eff.append(num / den)
    
    mean_eff = sum(eff) / len(eff)

    return mean_eff

fileList = ["GluGluToH_HToJPsiG_JPsiToMuMu_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16_Skim.root"]

print("Signal efficiency is: %.2f%%" % (100 * sig_eff(fileList)))