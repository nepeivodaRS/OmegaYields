#include "ptDistr.h"

// .L analyze_eff_check.C
// eff_check("./outputEff/mc_Eff_2march_injected.root");

TFile* FindFileFresh(const Char_t* fileName){
  // Find file
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName);
  if(file) {
    file->Close();
    delete file;
  }

  file = TFile::Open(fileName, "READ");

  if(!file)
    cout << "File : " << fileName << " was not found" << endl;

  return file;
}

void InitHists(){

}

void WriteToFile(TFile* outFile){
}

void eff_check(const Char_t* fileNameEff){
  fileEff = FindFileFresh(fileNameEff);
  fileEff->cd();
  hEffOmegaMB = (TH1D*)fileEff->Get("hEffOmegaMB");
  Int_t EffBin = hEffOmegaMB->GetXaxis()->FindBin(1.0);
  std::cout << hEffOmegaMB->GetBinContent(EffBin) << std::endl;
}
