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
#include "TRatioPlot.h"
#include "TCanvas.h"
#include "TString.h"

#include <iostream>
#include <fstream>

TFile* outFile;
TFile* inFile;

const Int_t nSignalTypes = 12;
const Char_t *SignalTypeName[nSignalTypes] = {"Generated", "BackgroundAllCasc", "AllCandidates", "RecLoose", "RecLooseFake", "RecLooseReal", "RecStandard", "RecStandardFake", "RecStandardReal", "RecTight", "RecTightFake", "RecTightReal"};

const Int_t nXAxisNames = 14;
const Char_t *XAxisNames[nXAxisNames] = {" cos(PA)", " R_{Casc}", " DCA_{+}", " DCA_{-}", " DCA_{CascPV}", " DCA_{V0Daughters}", " cos(V0PA)", " R_{V0}", " DCA_{V0PV}", " DCA_{Bach}", " R_{casc}/P_{t}^{casc}", " R_{V0}/P_{t}^{casc}", " DCA_{CascV0}", " DCA_{V0Bach}"};
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
TH1D *hV0BachDCA = { 0 };

TList *ListOfCutHists;
TList *ListFakeToRecRatiosTight;
TList *ListRecToGenRatiosStandart;
TList *ListFakeToAll;
TList *ListRealToFakeStandard;
TList *ListGenToRecTight;
TList *ListRealToFakeTight;
TLatex *uLatex;



