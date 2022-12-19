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

#include <AliAnalysisPIDCascadeEvent.h>
#include <AliAnalysisPIDCascadeTrack.h>
#include <AliAnalysisPIDCascadeV0.h>
#include <AliAnalysisPIDCascadeParticle.h>

#include "AliAnalysisPIDCascade.h"
#include "THStack.h"
#include "TCanvas.h"

#include "helpFunc.C"

#include <iostream>
#include <fstream>

// Binning used
const Int_t nPtBins = 19;
Double_t xBins[nPtBins+1] = 
  {0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4,
   2.6, 2.8, 3.0, 3.3, 3.6, 4.0, 4.5, 5.0, 5.5, 6.5};


const Int_t nMinvBins = 180;
Double_t minvBins[nMinvBins+1] = { 0 };

// Binning used in MB analysis
const Int_t nPtBinsMB = 7;
Double_t xBinsMB[nPtBinsMB+1] = 
  {1.00, 1.40, 1.80, 2.30, 2.80, 3.30, 3.80, 4.80};

// Bins for HM analysis

const Int_t nPtBinsHM = 13;
Double_t xBinsHM[nPtBinsHM+1] = 
  {0.60, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.50, 2.90, 
   3.40, 4.00, 5.00, 6.50 };

Bool_t tightPileUpEventCut = kFALSE; // not needed
Int_t minFastTracks = 1;

TFile* outFile;
TTree* tree;

// Event statistics
TH1I* hEventStat;

// Cascade statistics
TH1I* hCascStat;

Int_t nMB = 0, nHM = 0, nVHM = 0;

TH2D* hOmegaInvMassVsPtMB;
TH2D* hOmegaInvMassVsPtHM;
TH2D* hOmegaInvMassVsPtVHM;

TH1D* hGenOmegaMB;
TH1D* hGenOmegaHM;
TH1D* hGenOmegaVHM;

TH1D* hRecOmegaHMSd;
TH1D* hRecOmegaMBSd;
TH1D* hRecOmegaVHMSd;

TH1D* hEffOmegaMB;
TH1D* hEffOmegaHM;
TH1D* hEffOmegaVHM;
