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

// // Bins for HM analysis
// const Int_t nPtBinsHM = 13;
// Double_t xBinsHM[nPtBinsHM+1] = {0.60, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.50, 2.90, 3.40, 4.00, 5.00, 6.50 };

// Binning used in HM analysis for test
const Int_t nPtBinsHM = 7;
Double_t xBinsHM[nPtBinsHM+1] = {1.00, 1.40, 1.80, 2.30, 2.80, 3.30, 3.80, 4.80};

TCanvas * c1 = 0;

TList *PtHistList;
TList *ListOfPtHists;
TFile* outFile;
TFile* fileEff;

TH1D* hOmegaMBdef;
TH1D* hOmegaMBSideBand;
TH1D* hOmegaMBBGfix;
TH1D* hOmegaMBCombined;

TH1D* hOmegaMBMC;
TH1D* hOmegaHM;
TH1D* hOmegaVHM;

TH1D* hGenOmegaMB;
TH1D* hGenOmegaHM;
TH1D* hGenOmegaVHM;

TH1D* hEffOmegaMB;
TH1D* hEffOmegaHM;
TH1D* hEffOmegaVHM;

TH3D* hInvMassSumMC;
TH3D* hInvMassSumBGMC;