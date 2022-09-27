void mergeTrees (std::string data_dir = "data.nosync/") {
    
    const int nFiles = 3;
    std::string treesToMerge[nFiles] = {"6E527BF2-E7F7-DF4E-B1C2-5582052AD7E6_Skim.root", "75F407F6-DDE1-6644-A678-F6F95A52C31D_Skim.root", "CBD8B1C1-40C8-B744-84F0-92F471F0F946_Skim.root"};
        
    
    TList *list = new TList;
    TTree *tree;
    TFile *f;

    cout << endl << "Root files to be merged:" << endl;
    for (Int_t i = 0; i < nFiles; i++)
    {

        cout << data_dir + treesToMerge[i] << endl;
        f = new TFile((data_dir + treesToMerge[i]).c_str());
        tree = (TTree*)f->Get("Events");
        list->Add(tree);
    }

    
    TFile *mergedFile = new TFile("data.nosync/mergedFile.root", "RECREATE");
    TTree *mergedTree = TTree::MergeTrees(list);
    mergedTree->SetName("mergedEvents");
    mergedTree->Write();
    mergedFile->Close();
    cout << endl << "Merging successfully completed. Trees have been saved as mergedTree in mergedFile.root" << endl;
}