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

Bool_t tightPileUpEventCut = kFALSE; // not needed
Int_t minFastTracks = 1;

// Pt hists
TH1D *hOmegaMB, *hOmegaHM, *hOmegaVHM;
TH1D* hGenOmegaMB;
TH1D* hGenOmegaHM;
TH1D* hGenOmegaVHM;

Int_t nMB = 0, nHM = 0, nVHM = 0;

const Int_t nSpecies = 2;
const Char_t * pdgNameOmega[nSpecies] = {"Omega", "OmegaBar"};

// Cut hists
const Int_t nSignalTypes = 12;
const Char_t *SignalTypeName[nSignalTypes] = {"Generated", "BackgroundAllCasc", "AllCandidates", "RecLoose", "RecLooseFake", "RecLooseReal", "RecStandard", "RecStandardFake", "RecStandardReal", "RecTight", "RecTightFake", "RecTightReal"};

TH3D *hOmegaInvMassVsPtCuts[nSignalTypes] = { 0 };
TH3D *hLambdaInvMassVsPtCuts[nSignalTypes] = { 0 };

TH2D *hV0BachVsDCApn;

// Bach cuts
TH1D *hBachDCA[nSignalTypes] = { 0 };

// Cascade cuts
TH1D *hCascPA[nSignalTypes] = { 0 };
TH1D *hcascR[nSignalTypes] = { 0 };
TH1D *hCascPVDCA[nSignalTypes] = { 0 };
TH1D *hCascROverPt[nSignalTypes] = { 0 };
TH1D *hCascV0DCA[nSignalTypes] = { 0 };

// V0 cuts
TH1D *hPosDCA[nSignalTypes] = { 0 };
TH1D *hNegDCA[nSignalTypes] = { 0 };
TH1D *hV0DaughtersDCA[nSignalTypes] = { 0 };
TH1D *hV0PA[nSignalTypes] = { 0 };
TH1D *hV0R[nSignalTypes] = { 0 };
TH1D *hV0PVDCA[nSignalTypes] = { 0 };
TH1D *hV0ROverPt[nSignalTypes] = { 0 };
TH1D *hV0BachDCA[nSignalTypes] = { 0 };

TList *CutList;

// Inv mass 3d hist
TH3D *hOmegaInvMassVsPt[nSpecies] = { 0 };
TH3D *hOmegaInvMassVsPtTrue[nSpecies] = { 0 };
TH2D *hOmegaInconsistencyXi;

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

// // Bins for HM analysis

// const Int_t nPtBinsHM = 13;
// Double_t xBinsHM[nPtBinsHM+1] = 
//   {0.60, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.50, 2.90, 
//    3.40, 4.00, 5.00, 6.50 };

// Binning used in HM analysis for test
const Int_t nPtBinsHM = 7;
Double_t xBinsHM[nPtBinsHM+1] = 
  {1.00, 1.40, 1.80, 2.30, 2.80, 3.30, 3.80, 4.80};

// // Binning used in HM analysis for test
// const Int_t nPtBinsHM = 6;
// Double_t xBinsHM[nPtBinsHM+1] = 
//   {0.90, 1.60, 2.20, 2.60, 3.00, 3.80, 5.50};

// Bins for centrality analysis
const Int_t nCentrBins = 11;
Double_t xCentrBins[nCentrBins+1] = 
{0, 1, 5, 10, 20, 40, 50, 60, 70, 80, 90, 100};

// Event statistics
TH1I* hEventStat;

// Cascade statistics
TH1I* hCascStat;

TH1D* hVtxStatus;

TH1D* hNorm;

TList *PtHistList;
TList *InvMassList;

TFile* outFile;
TTree* tree;
