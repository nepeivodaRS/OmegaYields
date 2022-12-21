#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TROOT.h>

#include "THStack.h"
#include "TCanvas.h"
#include <TLatex.h>

#include <iostream>
#include <fstream>

// #include "helpFunc.C"

/*
  .L ptDistr.C+

  make_results("./outputAnal/mc_anal.root")
  
 */

const Double_t massOmega = 1.67245;
const Int_t nMinvBins = 180;
// Binning used in MB analysis
const Int_t nPtBinsMB = 7;
Double_t xBinsMB[nPtBinsMB+1] = {1.00, 1.40, 1.80, 2.30, 2.80, 3.30, 3.80, 4.80};

TCanvas * c1 = 0;

TList *PtHistList;
TFile* outFile;

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