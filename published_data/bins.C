#include "bins.h"

/*
.L bins.C

printBinningFromFile("fileOmegaStatOnly.root")
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

void PrintBinning(TH1D* hist)
{
  const Int_t nBinsX = hist->GetXaxis()->GetNbins();

  cout << "const Int_t nBins = " << nBinsX << ";" << endl
       << "Double_t bins[nBins+1] = " << endl
       << "{";
  for(Int_t bin = 1; bin <= nBinsX; bin++) {
   
    printf("%.2f, ", hist->GetXaxis()->GetBinLowEdge(bin));
    if (bin%10 == 0)
      cout << endl;
  }
 
  printf("%.2f", hist->GetXaxis()->GetBinUpEdge(nBinsX));
  cout << "}" << endl; 
}

void printBinningFromFile(const Char_t* fileNameData)
{
  TFile* fileData = FindFileFresh(fileNameData);
  if(!fileData)
    std::cout << " no data file was found" << std::endl;
  fileData->cd();

  //TH1D* hPtCorrData = (TH1D*)fileData->Get("hptcorr_norm_5");
  TH1D* hPtCorrData = (TH1D*)fileData->Get("hPtOmegaStatOnly_V0M_00000to00500");
  
  if(hPtCorrData){
    PrintBinning(hPtCorrData);
  }
  else
  {
    std::cout << " no histo was found" << std::endl;
  }
}