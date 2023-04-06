#include "analyze_tree_withEff.h"

/*
  To run code:
  ============

  aliroot -l

  .L analyze_tree_withEff.C+

  analyze_tree_withEff("./outputTreesDatLists/mc_tree_all.dat", "./outputEff/mc_Eff_6april_injected.root", "./outputAnal/mc_6april_effCorr.root", 0, 1)

  analyze_tree_withEff("./outputTreesDatLists/data_tree_all.dat", "./outputEff/mc_Eff_6april_injected.root", "./outputAnal/data_6april_effCorr.root", 0, 0)
*/

void InitHists(){
  // Initialization of histograms
  for(Int_t i = 0; i <= nMinvBins; i++) {
    minvBins[i] = -0.03 + i * 0.06/nMinvBins;
  }

  PtHistList = new TList();

  // hOmegaMB = new TH1D("hOmegaMB", "; p_{T} [GeV/c]",
  //      nPtBinsMB, xBinsMB);
  // hOmegaMB->Sumw2();

  // hOmegaHM = new TH1D("hOmegaHM", "; p_{T} [GeV/c]",
  //      nPtBinsHM, xBinsHM);
  // hOmegaHM->Sumw2();
  
  // hOmegaVHM = new TH1D("hOmegaVHM", "; p_{T} [GeV/c]",
  //       nPtBinsHM, xBinsHM);
  // hOmegaVHM->Sumw2();

  hGenOmegaMB = new TH1D("hGenOmegaMB", ";  p_{T} [GeV/c]",
        nPtBinsMB, xBinsMB);
  hGenOmegaMB->Sumw2();

  hGenOmegaHM = new TH1D("hGenOmegaHM", ";  p_{T} [GeV/c]",
          nPtBinsHM, xBinsHM);
  hGenOmegaHM->Sumw2();

  hGenOmegaVHM = new TH1D("hGenOmegaVHM", ";  p_{T} [GeV/c]",
          nPtBinsHM, xBinsHM);
  hGenOmegaVHM->Sumw2();

  // PtHistList->Add(hOmegaMB);
  // PtHistList->Add(hOmegaHM);
  // PtHistList->Add(hOmegaVHM);

  PtHistList->Add(hGenOmegaMB);
  PtHistList->Add(hGenOmegaHM);
  PtHistList->Add(hGenOmegaVHM);

  InvMassList = new TList();

  hOmegaInconsistencyXi = new TH2D("hOmegaInvMassVsPt", "M_{inv} inconsistency; M_{inv} - M_{#Xi} [GeV/c^{2}]; M_{inv} - M_{#Omega} [GeV/c^{2}]",
          nMinvBins, minvBins, nMinvBins, minvBins);
  hOmegaInconsistencyXi->Sumw2();
  InvMassList->Add(hOmegaInconsistencyXi);

  for(Int_t i = 0; i < nSpecies; i++) {
    hOmegaInvMassVsPt[i] = new TH3D(Form("hOmegaInvMassVsPt_%s", pdgNameOmega[i]), "M_{inv} vs p_{T}; p_{T}^{casc} [GeV/c]; M_{inv} - M_{#Omega} [GeV/c^{2}]; Centrality V0M",
            nPtBinsMB, xBinsMB, nMinvBins, minvBins, nCentrBins, xCentrBins);
    hOmegaInvMassVsPt[i]->Sumw2();

    hOmegaInvMassVsPtTrue[i] = new TH3D(Form("hOmegaInvMassVsPtTrue_%s", pdgNameOmega[i]), "M_{inv} vs p_{T}; p_{T}^{casc} [GeV/c]; M_{inv} - M_{#Omega} [GeV/c^{2}]; Centrality V0M",
            nPtBinsMB, xBinsMB, nMinvBins, minvBins, nCentrBins, xCentrBins);
    hOmegaInvMassVsPtTrue[i]->Sumw2();

    hOmegaInvMassVsPtTrueEffCorr[i] = new TH3D(Form("hOmegaInvMassVsPtTrueEffCorr_%s", pdgNameOmega[i]), "M_{inv} vs p_{T}; p_{T}^{casc} [GeV/c]; M_{inv} - M_{#Omega} [GeV/c^{2}]; Centrality V0M",
            nPtBinsMB, xBinsMB, nMinvBins, minvBins, nCentrBins, xCentrBins);
    hOmegaInvMassVsPtTrueEffCorr[i]->Sumw2();

    InvMassList->Add(hOmegaInvMassVsPt[i]);
    InvMassList->Add(hOmegaInvMassVsPtTrue[i]);
    InvMassList->Add(hOmegaInvMassVsPtTrueEffCorr[i]);
  }

  hEventStat = new TH1I("hEventStat","",3,0,3);

  hCascStat = new TH1I("hCascStat","",9,0,9);
}

void WriteToFile(TFile* outFile){
  // Write to file
  outFile->cd();
  outFile->WriteObject(PtHistList, "ListOfPtHists");
  outFile->WriteObject(InvMassList, "ListOfEInvMassHists");
  outFile->WriteObject(hEventStat, "EventStatistics");
  outFile->WriteObject(hCascStat, "hCascStat");
  outFile->WriteObject(hNorm, "hNorm");
  outFile->WriteObject(hEffOmegaMB, "hEffOmegaMB");
  outFile->Close();
}

void analyze_tree_withEff(const Char_t* inFileName,
    const Char_t* EffFileName,
		const Char_t* outFileName, 
		const Int_t maxEvents = 0,
    Bool_t isMC = kFALSE)
{
  // Open the input file and set up the classes.
  if(strstr(inFileName, ".dat")) {
    tree = ReadChainFromFile(inFileName, "tree", 0);
    hVtxStatus = ReadHistoVtxStatusFromFile(inFileName, 0);
  } else {
    TFile* inFile = TFile::Open(inFileName);
    if(!inFile){
      return;
    }
    tree = (TTree*)inFile->Get("tree");
    hVtxStatus = (TH1D*)inFile->Get("hVtxStatus");
  }

  // Get Data
  fileEff = FindFileFresh(EffFileName);
  fileEff->cd();
  hEffOmegaMB = (TH1D*)fileEff->Get("hEffOmegaMB");

  AliAnalysisPIDCascadeEvent* event = 0;
  TClonesArray* allCascades = 0;
  TClonesArray* generatedOmega = 0;
  tree->SetBranchAddress("event",   &event);
  tree->SetBranchAddress("allCascades", &allCascades);
  if(isMC)
    tree->SetBranchAddress("generatedOmega", &generatedOmega);

  outFile = new TFile(outFileName, "RECREATE");

  InitHists();

  // Loop over events

  Int_t nEvents = tree->GetEntries();
  cout << "Number of events: " << nEvents << std::endl;
  
  if(maxEvents>0 && maxEvents < nEvents) {
    nEvents = maxEvents;
    cout << "N events was reduced to: " << maxEvents << std::endl;
  }
  
  for(Int_t n = 0; n < nEvents; n++) {
    tree->GetEntry(n);

    if((n+1)%100000==0)
      cout << "Event: " << n+1 << "/" << nEvents << std::endl;

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

    if(isMC){
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
    }

    const Int_t nCascades = allCascades->GetEntries();

    for(Int_t i = 0; i < nCascades; i++) {

      Bool_t bIsRealOmegaCascade = kFALSE;
      Bool_t bPassedLoose = kFALSE;
      Bool_t bPassedStandard = kFALSE;
      Bool_t bPassedTight = kFALSE;

      hCascStat->Fill("Total", 1);

      AliAnalysisPIDCascade* cascade = (AliAnalysisPIDCascade*)allCascades->At(i);
      if(cascade->GetPtCasc() < 1.0 || cascade->GetPtCasc() > 4.80)
        continue;
      if(cascade->GetEtaCasc() > 0.8)
        continue;

      const Double_t dMassOmega     = cascade->GetIMO() - massOmega;
      const Double_t dMassXi = cascade->GetIMXi() - massXi;

      if(!CheckCascOmegaToXiMass(cascade))
        continue;

      hCascStat->Fill("#Xi check", 1);

      if(TMath::Abs(dMassOmega) > 0.029999999)
        continue;

      hCascStat->Fill("M_{#Omega} check", 1);

      if(isMC){
        bIsRealOmegaCascade = IsRealOmegaCascade(cascade);
        if(bIsRealOmegaCascade){
          hCascStat->Fill("Gen #Omega", 1);
        }
        else{
          hCascStat->Fill("Not #Omega", 1);
        }
      }

      Int_t hasFast = CascadeHasFastSignal(cascade);
      if(!hasFast)
        continue;

      hCascStat->Fill("Fast signal", 1);

      bPassedLoose = CheckCascLooseCuts(cascade);
      if(!bPassedLoose)
        continue;

      hCascStat->Fill("Loose", 1);

      bPassedStandard = CheckCascStandardCuts(cascade);
      if(bPassedStandard){
        hCascStat->Fill("Standard", 1);
        Int_t bin = 1;
        Int_t EffBin = hEffOmegaMB->GetXaxis()->FindBin(cascade->GetPtCasc());
        if(cascade->GetCharge() < 0){bin = 0;}{
          hOmegaInvMassVsPt[bin]->Fill(cascade->GetPtCasc(), dMassOmega, event->GetV0Mmultiplicity(), 1./hEffOmegaMB->GetBinContent(EffBin));
        }
        hOmegaInconsistencyXi->Fill(dMassXi, dMassOmega);
        // MC closure for signal inside reconstructed
        if(isMC){
          if(IsRealOmegaCascade(cascade)){
            if(cascade->GetCharge() < 0){bin = 0;}
              hOmegaInvMassVsPtTrue[bin]->Fill(cascade->GetPtCasc(), dMassOmega, event->GetV0Mmultiplicity());
              hOmegaInvMassVsPtTrueEffCorr[bin]->Fill(cascade->GetPtCasc(), dMassOmega, event->GetV0Mmultiplicity(), 1./hEffOmegaMB->GetBinContent(EffBin));            
          } 
        }
      }

      bPassedTight = CheckCascTightCuts(cascade);
      if(bPassedTight){
        hCascStat->Fill("Tight", 1);
      }
    }
  }

  if(maxEvents == 0)
    R__ASSERT(TMath::Nint(hVtxStatus->GetBinContent(3)) == nMB); // xcheck for broken logic
  
  outFile->cd();
  hNorm = new TH1D("hNorm", "MB: No vtx=-1, vtx rej=0, N MB=1, N HM=2, N VHM=3",
       5, -1.5, 3.5);
  hNorm->SetBinContent(1, hVtxStatus->GetBinContent(1));
  hNorm->SetBinContent(2, hVtxStatus->GetBinContent(2));
  hNorm->SetBinContent(3, nMB);
  hNorm->SetBinContent(4, nHM);
  hNorm->SetBinContent(5, nVHM);
  WriteToFile(outFile);
}
