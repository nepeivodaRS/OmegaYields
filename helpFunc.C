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

#include <iostream>
#include <fstream>

const Double_t massOmega = 1.67245;
const Double_t massXi = 1.32171;
const Double_t massLambda = 1.115683;

Bool_t CheckCascOmegaToXiMass(AliAnalysisPIDCascade* cascade) {
  // reject Omegas that are consistent with Xis
  const Double_t deltaM = cascade->GetIMXi() - massXi;
  if(TMath::Abs(deltaM)<0.008)
    return kFALSE;
  return kTRUE;
}

Bool_t IsRealOmegaCascade(AliAnalysisPIDCascade* cascade){
  AliAnalysisPIDCascadeV0* v0 = cascade->GetV0();
  AliAnalysisPIDCascadeTrack* bachelor = cascade->GetBachAnalysisTrack();
  const Int_t q = cascade->GetCharge();

  if(bachelor->GetMCPdgCode()    != q*321   || // K+
  bachelor->GetMCMotherPdgCode() != q*-3334 || // Omega+
  bachelor->GetMCMotherPrimary() != kTRUE)
    return kFALSE;

  const Int_t omegaLabel = bachelor->GetMCMotherLabel();

  // checks that daughters have same Lambda mother (pdg and label)
  if(v0->GetMCPdgCode() != q*-3122) // lambda
    return kFALSE;

  // Check that V0 daughters have the correct Omega grandmother
  // just need to check one daughter as they have common mother
  if(v0->GetPosAnalysisTrack()->GetMCPrimaryLabel() != omegaLabel)
    return kFALSE;
  return kTRUE;
}

Int_t CascadeHasFastSignal(AliAnalysisPIDCascade* casc) 
{
  AliAnalysisPIDCascadeTrack* tr[3] = {casc->GetBachAnalysisTrack(), 
                                       casc->GetV0()->GetPosAnalysisTrack(),
                                       casc->GetV0()->GetNegAnalysisTrack()};
  Int_t hasFast = 0;
  for (Int_t i = 0; i < 3; ++i) {
    if(tr[i]->GetStatus() & AliESDtrack::kITSrefit) {
      hasFast++;
      continue;
    }
    if(tr[i]->HasTOFPID())
      hasFast++;
  }
  return hasFast;
}

Double_t GetAbsImpactParameterXY(AliAnalysisPIDCascadeTrack* track)
{
  Double_t d1 = track->GetImpactParameter(0); 
  return TMath::Abs(d1);
}

Int_t AcceptEvent(AliAnalysisPIDCascadeEvent* eventIn, TClonesArray* trackArrayIn){
  // Apply all Vytautas cuts
  // Reject pileup
  //  + kNotPileupInSPD = 1,
  //  ! kNotPileupInMV = 2,
  //  ! kNotPileupInMB = 4,
  //  ! kINELgtZERO = 8,
  //  ! kNoInconsistentVtx = 16,
  //  ! kNoV0Asym = 32,
  //  ! kVertexSelected2015pp=64,
  //  ! kSPDandTrkVtxExists=128,
  //  ! kPassProximityCut=256,
  AliAnalysisPIDCascadeEvent::SetCheckFlag(1);
  if(!eventIn->AcceptEvent(kFALSE))
    return -2;

  //Check if has a vertex
  if(!eventIn->HasVertex())
    return -1;

  //  ! kNotPileupInSPD = 1,
  //  ! kNotPileupInMV = 2,
  //  ! kNotPileupInMB = 4,
  //  ! kINELgtZERO = 8,
  //  ! kNoInconsistentVtx = 16,
  //  ! kNoV0Asym = 32,
  //  + kVertexSelected2015pp=64,
  //  + kSPDandTrkVtxExists=128,
  //  + kPassProximityCut=256,
  AliAnalysisPIDCascadeEvent::SetCheckFlag(448);
  if(!eventIn->AcceptEvent(kFALSE))
    return -1;

  //Check if has a vertex within 10 cm
  if(TMath::Abs(eventIn->GetVertexZ())>10)
    return 0;

  //If no tracks in the event, also throw it away
  if(trackArrayIn) {
    if(trackArrayIn->GetEntries()<1) {
    return 1;
    }
  }
  return 1;
}

Bool_t CheckCascLooseCuts(AliAnalysisPIDCascade* cascade) {
  // Pt check
  if(cascade->GetPtCasc() < 0.15)
    return kFALSE;

  AliAnalysisPIDCascadeV0* v0 = cascade->GetV0();
  AliAnalysisPIDCascadeTrack* bachelor = cascade->GetBachAnalysisTrack();
  const Double_t bachDCA = GetAbsImpactParameterXY(bachelor);
  const Double_t negDCA  = GetAbsImpactParameterXY(v0->GetNegAnalysisTrack());
  const Double_t posDCA  = GetAbsImpactParameterXY(v0->GetPosAnalysisTrack());
  const Double_t cascPA     = cascade->GetCascCosinePA();
  const Double_t cascR      = cascade->GetCascRadius();
  const Double_t dMassOmega     = cascade->GetIMO() - massOmega;
  const Double_t dMassLambda = v0->GetIML() - massLambda;
  AliAnalysisPIDCascadeTrack* tr[3] = {cascade->GetBachAnalysisTrack(), 
                                       cascade->GetV0()->GetPosAnalysisTrack(),
                                       cascade->GetV0()->GetNegAnalysisTrack()};
  if(cascade->GetCharge() < 0) {
    AliAnalysisPIDCascadeTrack* dummy = tr[1];
    tr[1] = tr[2];
    tr[2] = dummy;
  }

  // Track selection
  for (Int_t i = 0; i < 3; ++i) {
    if(tr[i]->GetTPCNcls() < 70 ||
      !(tr[i]->GetStatus() & AliESDtrack::kTPCrefit) ||
      TMath::Abs(tr[i]->GetEta()) > 0.8 ||
      !tr[i]->HasTPCPID())
      return kFALSE;
  }
  // hCascStat->Fill("Track Loose", 1);

  // TPC PID cuts
  if(TMath::Abs(tr[0]->GetNSigmaKaonTPC()) > 5 ||
    TMath::Abs(tr[1]->GetNSigmaPionTPC()) > 5 ||
    TMath::Abs(tr[2]->GetNSigmaProtonTPC()) > 5)
    return kFALSE;
  // hCascStat->Fill("TPC Loose", 1);

  // TOF PID bachelor cut
  if(bachelor->HasTOFPID())
  if(TMath::Abs(bachelor->GetNSigmaKaonTOF()) > 5)
    return kFALSE;
  // hCascStat->Fill("K TOF Loose", 1);

  // Topo cuts
  if(cascPA < 0.97 || 
    cascR < 0.6 || cascR > 100 ||
    posDCA < 0.02 ||
    negDCA < 0.02 ||
    //cascade->GetCascDCA() > 1.5 ||
    TMath::Abs(dMassLambda) > 0.008 ||
    v0->GetDCAV0Daughters() > 1.7 ||
    v0->GetV0CosinePA() < 0.97 ||
    v0->GetRadius() < 1.3 || v0->GetRadius() > 100 ||
    v0->GetDCAPV() < 0.06 ||
    bachDCA < 0.05 ||
    cascade->GetCascRadius()/cascade->GetPtCasc() > 16 ||
    v0->GetRadius()/cascade->GetPtCasc() > 45 ||
    //cascade->GetCascDCAPV() > 1.0 ||
    cascade->GetV0DCA() > 11)
    return kFALSE;
  return kTRUE;
}

Bool_t CheckCascStandardCuts(AliAnalysisPIDCascade* cascade) {
  // Pt check
  if(cascade->GetPtCasc() < 0.15)
    return kFALSE;

  AliAnalysisPIDCascadeV0* v0 = cascade->GetV0();
  AliAnalysisPIDCascadeTrack* bachelor = cascade->GetBachAnalysisTrack();
  const Double_t bachDCA = GetAbsImpactParameterXY(bachelor);
  const Double_t negDCA  = GetAbsImpactParameterXY(v0->GetNegAnalysisTrack());
  const Double_t posDCA  = GetAbsImpactParameterXY(v0->GetPosAnalysisTrack());
  const Double_t cascPA     = cascade->GetCascCosinePA();
  const Double_t cascR      = cascade->GetCascRadius();
  const Double_t dMassOmega     = cascade->GetIMO() - massOmega;
  const Double_t dMassLambda = v0->GetIML() - massLambda;
  AliAnalysisPIDCascadeTrack* tr[3] = {cascade->GetBachAnalysisTrack(), 
                                       cascade->GetV0()->GetPosAnalysisTrack(),
                                       cascade->GetV0()->GetNegAnalysisTrack()};
  if(cascade->GetCharge() < 0) {
    AliAnalysisPIDCascadeTrack* dummy = tr[1];
    tr[1] = tr[2];
    tr[2] = dummy;
  }

  // Track selection
  for (Int_t i = 0; i < 3; ++i) {
    if(tr[i]->GetTPCNcls() < 80 ||
      !(tr[i]->GetStatus() & AliESDtrack::kTPCrefit) ||
      TMath::Abs(tr[i]->GetEta()) > 0.8 ||
      !tr[i]->HasTPCPID())
      return kFALSE;
  }
  // hCascStat->Fill("Track Standard", 1);

  // TPC PID cuts
  if(TMath::Abs(tr[0]->GetNSigmaKaonTPC()) > 4 ||
    TMath::Abs(tr[1]->GetNSigmaPionTPC()) > 4 ||
    TMath::Abs(tr[2]->GetNSigmaProtonTPC()) > 4)
    return kFALSE;
  // hCascStat->Fill("TPC Standard", 1);

  // TOF PID bachelor cut
  if(bachelor->HasTOFPID())
  if(TMath::Abs(bachelor->GetNSigmaKaonTOF()) > 4)
    return kFALSE;
  // hCascStat->Fill("K TOF Standard", 1);

  // Topo cuts
  if(cascPA < 0.98 || 
    cascR < 0.6 || cascR > 100 ||
    posDCA < 0.03 ||
    negDCA < 0.03 ||
    //cascade->GetCascDCA() > 1.5 ||
    TMath::Abs(dMassLambda) > 0.006 ||
    v0->GetDCAV0Daughters() > 1.6 ||
    v0->GetV0CosinePA() < 0.98 ||
    v0->GetRadius() < 1.4 || v0->GetRadius() > 100 ||
    v0->GetDCAPV() < 0.07 ||
    bachDCA < 0.05 ||
    cascade->GetCascRadius()/cascade->GetPtCasc() > 15 ||
    v0->GetRadius()/cascade->GetPtCasc() > 40 ||
    //cascade->GetCascDCAPV() > 1.0 ||
    cascade->GetV0DCA() > 10)
    return kFALSE;
  return kTRUE;
}

Bool_t CheckCascTightCuts(AliAnalysisPIDCascade* cascade) {
  // Pt check
  if(cascade->GetPtCasc() < 0.15)
    return kFALSE;

  AliAnalysisPIDCascadeV0* v0 = cascade->GetV0();
  AliAnalysisPIDCascadeTrack* bachelor = cascade->GetBachAnalysisTrack();
  const Double_t bachDCA = GetAbsImpactParameterXY(bachelor);
  const Double_t negDCA  = GetAbsImpactParameterXY(v0->GetNegAnalysisTrack());
  const Double_t posDCA  = GetAbsImpactParameterXY(v0->GetPosAnalysisTrack());
  const Double_t cascPA     = cascade->GetCascCosinePA();
  const Double_t cascR      = cascade->GetCascRadius();
  const Double_t dMassOmega     = cascade->GetIMO() - massOmega;
  const Double_t dMassLambda = v0->GetIML() - massLambda;
  AliAnalysisPIDCascadeTrack* tr[3] = {cascade->GetBachAnalysisTrack(), 
                                       cascade->GetV0()->GetPosAnalysisTrack(),
                                       cascade->GetV0()->GetNegAnalysisTrack()};
  if(cascade->GetCharge() < 0) {
    AliAnalysisPIDCascadeTrack* dummy = tr[1];
    tr[1] = tr[2];
    tr[2] = dummy;
  }

  // Track selection
  for (Int_t i = 0; i < 3; ++i) {
    if(tr[i]->GetTPCNcls() < 80 ||
      !(tr[i]->GetStatus() & AliESDtrack::kTPCrefit) ||
      TMath::Abs(tr[i]->GetEta()) > 0.8 ||
      !tr[i]->HasTPCPID())
      return kFALSE;
  }
  // hCascStat->Fill("Track Tight", 1);

  // TPC PID cuts
  if(TMath::Abs(tr[0]->GetNSigmaKaonTPC()) > 4 ||
    TMath::Abs(tr[1]->GetNSigmaPionTPC()) > 4 ||
    TMath::Abs(tr[2]->GetNSigmaProtonTPC()) > 4)
    return kFALSE;
  // hCascStat->Fill("TPC Tight", 1);

  // TOF PID bachelor cut
  if(bachelor->HasTOFPID())
  if(TMath::Abs(bachelor->GetNSigmaKaonTOF()) > 4)
    return kFALSE;
  // hCascStat->Fill("K TOF Tight", 1);

  // Topo cuts
  if(cascPA < 0.99 || 
    cascR < 0.6 || cascR > 100 ||
    posDCA < 0.04 ||
    negDCA < 0.04 ||
    //cascade->GetCascDCA() > 1.5 ||
    TMath::Abs(dMassLambda) > 0.005 ||
    v0->GetDCAV0Daughters() > 1.5 ||
    v0->GetV0CosinePA() < 0.99 ||
    v0->GetRadius() < 1.5 || v0->GetRadius() > 100 ||
    v0->GetDCAPV() < 0.08 ||
    bachDCA < 0.05 ||
    cascade->GetCascRadius()/cascade->GetPtCasc() > 14 ||
    v0->GetRadius()/cascade->GetPtCasc() > 35 ||
    //ascade->GetCascDCAPV() > 1.0 ||
    cascade->GetV0DCA() > 9)
    return kFALSE;
  return kTRUE;
}

TFile* FindFile(const Char_t* fileName){
  // Find file
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName);
  if(file) {
    return file;
  }
  file = TFile::Open(fileName, "READ");
  if(!file)
    cout << "File : " << fileName << " was not found" << endl;
  return file;
}

TChain* ReadChainFromFile(const char *fileIn, const char *treeName, const char *fName, Int_t maxFiles = -1, Int_t startFile = 0){
  // Create the chain
  TChain* chain = new TChain(treeName);
  // Open the input stream
  ifstream in;
  in.open(fileIn);
  // Read the input list of files and add them to the chain
  TString currentFile;
  Int_t counter=0;
  while(in.good()) {
        in >> currentFile;
    if (fName) {
      currentFile+="#";
      currentFile+=fName;
    }
    if (!currentFile.Contains("root")) continue; // protection
    counter++;
    if (counter<startFile) continue;
    if (maxFiles>0 && counter>maxFiles+startFile) break;
    TFile * f = FindFile(currentFile.Data());
    if (f){
      chain->Add(currentFile.Data());
    }
    delete f;
  }

  in.close();

  return chain;
}

void ReadHistoFromFile(const char *fileIn, TH1D *inHist, Int_t maxFiles = -1, Int_t startFile = 0){
  // Open the input stream
  ifstream in;
  in.open(fileIn);
  // Read the input list of files and add them to the chain
  TString currentFile;
  Int_t counter=0;
  while(in.good()) {
        in >> currentFile;
    if (!currentFile.Contains("root")) continue; // protection
    counter++;
    if (counter<startFile) continue;
    if (maxFiles>0 && counter>maxFiles+startFile) break;
    TFile * f = FindFile(currentFile.Data());
    if (f){
      TH1D* tempHist = (TH1D*)f->Get("hVtxStatus");
      TH1D* hVtxStatus = (TH1D*)tempHist->Clone();
      inHist->Add(hVtxStatus);
    }
    delete f;
  }
  in.close();
}