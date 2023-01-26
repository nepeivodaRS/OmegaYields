#include "mc_efficiency.h"

/*
  To run code:
  ============

  root
  .L mc_efficiency.C+

  mc_efficiency("./outputTreesMC/mc_tree_all.dat", "./outputEff/mc_Eff_sameBinning.root", 0)

*/

void InitHists(){ 
  for(Int_t i = 0; i <= nMinvBins; i++) {
    minvBins[i] = -0.03 + i * 0.06/nMinvBins;
  }

  hGenOmegaMB = new TH1D("hGenOmegaMB", ";  p_{T} [GeV/c]",
        nPtBinsMB, xBinsMB);
  hGenOmegaMB->Sumw2();

  hGenOmegaHM = new TH1D("hGenOmegaHM", ";  p_{T} [GeV/c]",
          nPtBinsHM, xBinsHM);
  hGenOmegaHM->Sumw2();

  hGenOmegaVHM = new TH1D("hGenOmegaVHM", ";  p_{T} [GeV/c]",
          nPtBinsHM, xBinsHM);
  hGenOmegaVHM->Sumw2();

  hRecOmegaMBSd = new TH1D("hRecOmegaMBSd", ";  p_{T} [GeV/c]",
        nPtBinsMB, xBinsMB);
  hRecOmegaMBSd->Sumw2();
  
  hRecOmegaHMSd = new TH1D("hRecOmegaHMSd", ";  p_{T} [GeV/c]",
          nPtBinsHM, xBinsHM);
  hRecOmegaHMSd->Sumw2();
  
  hRecOmegaVHMSd = new TH1D("hRecOmegaVHMSd", ";  p_{T} [GeV/c]",
          nPtBinsHM, xBinsHM);
  hRecOmegaVHMSd->Sumw2();
  
  hEffOmegaMB = new TH1D("hEffOmegaMB", ";  p_{T} [GeV/c]",
        nPtBinsMB, xBinsMB);
  hEffOmegaMB->Sumw2();
  
  hEffOmegaHM = new TH1D("hEffOmegaHM", ";  p_{T} [GeV/c]",
          nPtBinsHM, xBinsHM);
  hEffOmegaHM->Sumw2();
  
  hEffOmegaVHM = new TH1D("hEffOmegaVHM", ";  p_{T} [GeV/c]",
          nPtBinsHM, xBinsHM);
  hEffOmegaVHM->Sumw2();

  hOmegaInvMassVsPtMB = new TH2D("hOmegaInvMassVsPtMB", "M_{inv} vs p_{T}; p_{T}^{casc} [GeV/c]; M_{inv} - M_{#Omega} [GeV/c^{2}]",
          nPtBinsMB, xBinsMB, nMinvBins, minvBins);
  hOmegaInvMassVsPtMB->Sumw2();

  hOmegaInvMassVsPtHM = new TH2D("hOmegaInvMassVsPtHM", "M_{inv} vs p_{T}; p_{T}^{V0} [GeV/c]; M_{inv} - M_{#Omega}[GeV/c^{2}]",
          nPtBinsHM, xBinsHM, nMinvBins, minvBins);
  hOmegaInvMassVsPtHM->Sumw2();

  hOmegaInvMassVsPtVHM = new TH2D("hOmegaInvMassVsPtVHM", "M_{inv} vs p_{T}; p_{T}^{V0} [GeV/c]; M_{inv} - M_{#Omega} [GeV/c^{2}]",
          nPtBinsHM, xBinsHM, nMinvBins, minvBins);
  hOmegaInvMassVsPtVHM->Sumw2();

  hEventStat = new TH1I("hEventStat","",3,0,3);
  hCascStat = new TH1I("hCascStat","",4,0,4);
}

void WriteToFile(TFile* outFile){
  // Write to file
  outFile->cd();
  outFile->Write();
  outFile->Close();
}

void mc_efficiency(const Char_t* inFileName,
		   const Char_t* outFileName,
		   const Int_t maxEvents = 0)
{ 
  // Open the input file and set up the classes.
  if(strstr(inFileName, ".dat")) {
    tree = ReadChainFromFile(inFileName, "tree", 0);
  } else {
    TFile* inFile = TFile::Open(inFileName);
    if(!inFile){
      return;
    }
    tree = (TTree*)inFile->Get("tree");
  }

  AliAnalysisPIDCascadeEvent* event = 0;
  TClonesArray* allCascades = 0;
  TClonesArray* generatedOmega = 0;
  tree->SetBranchAddress("event",   &event);
  tree->SetBranchAddress("allCascades", &allCascades);
  tree->SetBranchAddress("generatedOmega", &generatedOmega);

  outFile = new TFile(outFileName, "RECREATE");

  InitHists();

  // Loop over events

  Int_t nEvents = tree->GetEntries();
  cout << "Number of events: " << nEvents << endl;
  
  if(maxEvents>0 && maxEvents < nEvents) {
    nEvents = maxEvents;
    cout << "N events was reduced to: " << maxEvents << endl;
  }
  
  for(Int_t n = 0; n < nEvents; n++) {
    tree->GetEntry(n);
    if((n+1)%100000==0)
      cout << "Event: " << n+1 << "/" << nEvents << endl;

    Bool_t HMevent = kFALSE;
    Bool_t VHMevent = kFALSE;

    // increase event copunters
    hEventStat->Fill("MB", 1);
    nMB++;

    if(event->GetV0Mmultiplicity() < 10){
      nHM++;
      HMevent = kTRUE;
      hEventStat->Fill("HM", 1);
    }
    if(event->GetV0Mmultiplicity() < 1){
      nVHM++;
      VHMevent = kTRUE;
      hEventStat->Fill("VHM", 1);
    }

    const Int_t nGenTracks = generatedOmega->GetEntries();
    for(Int_t i = 0; i < nGenTracks; i++) {
      AliAnalysisPIDCascadeParticle* trackMC = (AliAnalysisPIDCascadeParticle*)generatedOmega->At(i);
      if(TMath::Abs(trackMC->GetPt()) < 1.0 || TMath::Abs(trackMC->GetPt()) > 4.8)
        continue;
      if(TMath::Abs(trackMC->GetEta()) > 0.8) // overcheck
        continue;
      hGenOmegaMB->Fill(trackMC->GetPt());
      if(HMevent)
        hGenOmegaHM->Fill(trackMC->GetPt());
      if(VHMevent)
        hGenOmegaVHM->Fill(trackMC->GetPt());
    }
    
    //Analyze omegas

    const Int_t nCascades = allCascades->GetEntries();
    for(Int_t i = 0; i < nCascades; i++) {
      
      Bool_t bPassedStandard = kFALSE;

      AliAnalysisPIDCascade* cascade = (AliAnalysisPIDCascade*)allCascades->At(i);
      if(cascade->GetPtCasc() < 1.0 || cascade->GetPtCasc() > 4.80)
        continue;
      if(cascade->GetEtaCasc() > 0.8)
        continue;
      if(!CheckCascOmegaToXiMass(cascade))
        continue;

      hCascStat->Fill("#Xi check", 1);

      const Double_t deltaM     = cascade->GetIMO() - massOmega;
      if(TMath::Abs(deltaM) > 0.029999999)
        continue;

      hCascStat->Fill("M_{#Omega} check", 1);

      Int_t hasFast = CascadeHasFastSignal(cascade);
    
      if(!hasFast)
        continue;

      hCascStat->Fill("Fast signal", 1);

      bPassedStandard = CheckCascStandardCuts(cascade);
      if(!bPassedStandard)
        continue;

      hCascStat->Fill("Standard cuts", 1);

      // Was reconstrution correct?
      if(!IsRealOmegaCascade(cascade))
        continue;

      hOmegaInvMassVsPtMB->Fill(cascade->GetPtCasc(), deltaM);
      if(HMevent)
        hOmegaInvMassVsPtHM->Fill(cascade->GetPtCasc(), deltaM);
      if(VHMevent)
        hOmegaInvMassVsPtVHM->Fill(cascade->GetPtCasc(), deltaM);

      hRecOmegaMBSd->Fill(cascade->GetPtCasc());
      if(HMevent)
        hRecOmegaHMSd->Fill(cascade->GetPtCasc());
      if(VHMevent)
        hRecOmegaVHMSd->Fill(cascade->GetPtCasc());
    }
  }

  hEffOmegaMB->Divide(hRecOmegaMBSd, hGenOmegaMB, 1, 1, "B");
  hEffOmegaHM->Divide(hRecOmegaHMSd, hGenOmegaHM, 1, 1, "B");
  hEffOmegaVHM->Divide(hRecOmegaVHMSd, hGenOmegaVHM, 1, 1, "B");

  WriteToFile(outFile);
}