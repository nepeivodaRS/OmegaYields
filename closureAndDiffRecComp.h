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
#include <TKey.h>

TFile* outFile;
TCanvas* cRatioPlot;
// Binning used in MB analysis
TH1D* hOmegaMBPublished;

TH1D* hOmegaMBdef;
TH1D* hOmegaMBSideBand;
TH1D* hOmegaMBBGfix;
TH1D* hOmegaMBCombined;

const Int_t nPtBinsMB = 7;
Double_t xBinsMB[nPtBinsMB+1] = {1.00, 1.40, 1.80, 2.30, 2.80, 3.30, 3.80, 4.80};