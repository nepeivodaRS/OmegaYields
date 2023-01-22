#include "inconsistencyXi.h"

/*
.L inconsistencyXi.C

inconsistency("mc_anal_inconsistency_after_cuts_wider_range.root")
*/

TFile* FindFileFresh(const Char_t* fileName)
{
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

void inconsistency(const Char_t* fileNameData)
{
  TFile* fileData = FindFileFresh(fileNameData);
  if(!fileData)
    std::cout << " no data file was found" << std::endl;
  fileData->cd();
  TList *ListOfHists = (TList*)fileData->Get("ListOfEInvMassHists");
  TH2D* hOmegaInvMassVsPt = (TH2D*)ListOfHists->FindObject("hOmegaInvMassVsPt");
  hOmegaInvMassVsPt->SetStats("e");
  gStyle->SetOptStat("e");
  hOmegaInvMassVsPt->SetTitle("");
  if(hOmegaInvMassVsPt){
    hOmegaInvMassVsPt->GetXaxis()->SetRangeUser(-0.05, 0.1);
    hOmegaInvMassVsPt->GetYaxis()->SetRangeUser(-0.06, 0.06);
  }
  else
  {
    std::cout << " no histo was found" << std::endl;
  }
  hOmegaInvMassVsPt->Draw("colz");
}