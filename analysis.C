#include <stdexcept>


void analysis(const char *fileName = "GluGluToH_HToJPsiG_JPsiToMuMu_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16_Skim.root") {
    
    gROOT->SetStyle("Default");
    //gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(0);
	gStyle->SetStatBorderSize(1);
    ROOT::EnableImplicitMT();

    // variables
    Int_t bins = 256;



    TFile *f = new TFile(fileName);
    TTree *tree = (TTree*)f->Get("Events");

    TBranch *muon_pt = tree->GetBranch("Muon_pt");

    


    TCanvas *c1 = new TCanvas("c1", "Muon_pt");
    TH1D *h_muon_pt = new TH1D("h_muon_pt", "", bins, 0, 40);

    for (Int_t i = 0; i<muon_pt->GetEntries(); i++) // parsing the file
        {
            h_muon_pt->Fill(muon_pt->GetEntry(i));
        }

    c1->cd();
    h_muon_pt->GetYaxis()->SetTitle("counts"); h_muon_pt->GetXaxis()->SetTitle("p_T [GeV]");
    h_muon_pt->GetYaxis()->SetRangeUser(1., 1e5);
    h_muon_pt->SetLineColor(kBlue);
    h_muon_pt->SetFillColorAlpha(kBlue, 0.4);
    h_muon_pt->Draw();
    h_muon_pt->Draw("SAME");
}