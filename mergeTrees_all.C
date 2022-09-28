

void mergeTrees_all(const char *dirname="data.nosync/", const char *ext=".root", TString inputTree = "Events", TString outputTree = "mergedEvents", TString outputFile = "mergedFile.root")
{
    
    

    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();

    if (files)
    {   
        TFile *f;
        TList *list = new TList;
        TSystemFile *file;
        TString fname;

        Int_t kFile = 0;
        Int_t nEntries = 0;
        Int_t totEntries;
        bool correctFiles_found = FALSE;

        cout << "Merging files in " << dirname << endl << endl;
        TIter next(files);
        while ((file=(TSystemFile*)next()))
        {
            fname = file->GetName();
            
            if (!file->IsDirectory() && fname.EndsWith(ext) && fname != ".DS_Store" && fname != outputFile)
            {
                cout << " * Adding " << fname.Data() << " for merging..." << endl;
                f = new TFile(dirname + fname);
                TTree *tree = (TTree*)f->Get(inputTree);
                cout << "nEvents: " << tree->GetEntries() << endl;
                nEntries += tree->GetEntries();
                list->Add(tree);
                kFile += 1;
                correctFiles_found = TRUE;
            }
        }
        
        if (correctFiles_found)
        {
            TFile *mergedFile = new TFile(dirname + outputFile, "RECREATE");
            cout << endl << "Merging...";
            TTree *mergedTree = TTree::MergeTrees(list);
            totEntries = mergedTree->GetEntries();

            mergedTree->SetName(outputTree);
            mergedTree->Write();
            mergedFile->Close();

            if (nEntries == totEntries) { cout << endl << endl << "Merging of " << kFile << " files successfully completed. Total of " << totEntries << " events have been saved as a `"<< outputTree << "` TTree in " << outputFile << endl << endl; }
        }
        else { cout << "No adequate files found, make sure dirname is properly written and there are files with the extension `" << ext << "` inside." << endl; }
    }
}