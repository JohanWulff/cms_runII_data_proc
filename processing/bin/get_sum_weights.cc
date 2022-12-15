#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1.h"


std::cout << "Reading from file: " << fname << "\n";

double get_sum_weight(std::string fname){
    TFile* in_file = TFile::Open(fname.c_str());
    TH1F *h_eff = (TH1F*)in_file.Get("h_eff");
    Double_t binContent = heff->GetBinContent(1);
    return binContent;
}