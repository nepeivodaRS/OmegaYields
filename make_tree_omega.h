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

const Int_t nSpecies = 2;
const Char_t * pdgNameOmega[nSpecies] = {"Omega", "OmegaBar"};

const Int_t nMinvBins = 180;
Double_t minvBins[nMinvBins+1] = { 0 };

// Binning INV mass
const Int_t nPtBins = 19;
Double_t xBins[nPtBins+1] = 
{0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.3, 3.6, 4.0, 4.5, 5.0, 5.5, 6.5};

// Binning in MB pt analysis
const Int_t nPtBinsMB = 7;
Double_t xBinsMB[nPtBinsMB+1] = {1.00, 1.40, 1.80, 2.30, 2.80, 3.30, 3.80, 4.80};

// Binning in HM pt analysis
const Int_t nPtBinsHM = 13;
Double_t xBinsHM[nPtBinsHM+1] = {0.60, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.50, 2.90, 3.40, 4.00, 5.00, 6.50 };

TH2D* hInvMassOmegaVsPt[nSpecies] = { 0 };
TH1D* hVtxStatus;
TH1D* hVtxZ;

// Event statistics
TH1I* hEventStat;

TTree* treeOut;

TClonesArray *mcArrayOut, *cascadeArrayOut;

TTree* tree = 0;

AliAnalysisPIDCascadeEvent* event = 0;
TClonesArray* trackArray = 0;
TClonesArray* cascadeArray = 0;
TClonesArray* mcTrackArray = 0;

TH1D* hGenOmegaMB[nSpecies] = { 0 };
TH1D* hGenOmegaHM[nSpecies] = { 0 };
TH1D* hGenOmegaVHM[nSpecies] = { 0 };
