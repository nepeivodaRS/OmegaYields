#include "TopoSelectionStatistics.h"

/*
.L TopoSelectionStatistics.C

TopoSelectionStatistics("mc_anal_22.root")
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

void TopoSelectionStatistics(const Char_t* fileNameData)
{
  TFile* fileData = FindFileFresh(fileNameData);
  if(!fileData)
    std::cout << " no data file was found" << std::endl;
  fileData->cd();
  TH1I* hStat = (TH1I*)fileData->Get("hCascStat");
  gStyle->SetOptStat(0);
  if(hStat){
    hStat->GetXaxis()->SetBinLabel(7, "Standard");
    hStat->SetStats(0);
    hStat->SetTitle("");
  }
  else
  {
    std::cout << " no histo was found" << std::endl;
  }
  hStat->Draw("colz");
}