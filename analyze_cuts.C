#include "analyze_cuts.h"

/*
  To run code:
  ============

  aliroot -l

  .L analyze_cuts.C+

  analyze_cuts("./outputTreesMC/mc_tree_all.dat", "./outputCuts/mc_cuts_27.root", 0, 1)
*/
void FillCutHists(Int_t SigType, AliAnalysisPIDCascade* cascade, AliAnalysisPIDCascadeEvent* event){
  // Fill CutList hists with signals of specific type
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

  hOmegaInvMassVsPtCuts[SigType]->Fill(cascade->GetPtCasc(), dMassOmega, event->GetV0Mmultiplicity());
  hLambdaInvMassVsPtCuts[SigType]->Fill(v0->GetPt(), dMassLambda, event->GetV0Mmultiplicity());

  hCascPA[SigType]->Fill(cascPA);
  hcascR[SigType]->Fill(cascR);
  hPosDCA[SigType]->Fill(posDCA);
  hNegDCA[SigType]->Fill(negDCA);
  hBachDCA[SigType]->Fill(bachDCA);
  hCascPVDCA[SigType]->Fill(cascade->GetCascDCAPV());
  hCascROverPt[SigType]->Fill(cascade->GetCascRadius()/cascade->GetPtCasc());
  hCascV0DCA[SigType]->Fill(cascade->GetV0DCA());
  hV0DaughtersDCA[SigType]->Fill(v0->GetDCAV0Daughters());
  hV0PA[SigType]->Fill(v0->GetV0CosinePA());
  hV0R[SigType]->Fill(v0->GetRadius());
  hV0PVDCA[SigType]->Fill(v0->GetDCAPV());
  hV0ROverPt[SigType]->Fill(v0->GetRadius()/cascade->GetPtCasc());
  hV0BachDCA[SigType]->Fill(cascade->GetCascDCA()); // the distance between the V0 and the bachelor track at the Secondary Vertex. Which should be small for true cascades.
}

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

    InvMassList->Add(hOmegaInvMassVsPt[i]);
    InvMassList->Add(hOmegaInvMassVsPtTrue[i]);
  }

  CutList = new TList();

  for(Int_t i = 0; i < nSignalTypes; i++) {
    hCascPA[i] = new TH1D(Form("hCascPA_%s", SignalTypeName[i]), "; cos(PA)",
          60, 0.945, 1.005);
    hCascPA[i]->Sumw2();

    hcascR[i] = new TH1D(Form("hcascR_%s", SignalTypeName[i]), "; R_{Casc}",
          100, 0, 20);
    hcascR[i]->Sumw2();

    hPosDCA[i] = new TH1D(Form("hPosDCA_%s", SignalTypeName[i]), "; DCA_{+}",
          50, 0, 10);
    hPosDCA[i]->Sumw2();

    hNegDCA[i] = new TH1D(Form("hNegDCA_%s", SignalTypeName[i]), "; DCA_{-}",
          50, 0, 10);
    hNegDCA[i]->Sumw2();

    hCascPVDCA[i] = new TH1D(Form("hCascPVDCA_%s", SignalTypeName[i]), "; DCA_{CascPV}",
          100, 0, 10);
    hCascPVDCA[i]->Sumw2();

    hV0DaughtersDCA[i] = new TH1D(Form("hV0DaughtersDCA_%s", SignalTypeName[i]), "; DCA_{V0Daughters}",
          10, 0, 2);
    hV0DaughtersDCA[i]->Sumw2();

    hV0PA[i] = new TH1D(Form("hV0PA_%s", SignalTypeName[i]), "; cos(V0PA)",
          60, 0.945, 1.005);
    hV0PA[i]->Sumw2();

    hV0R[i] = new TH1D(Form("hV0R_%s", SignalTypeName[i]), "; R_{V0}",
          100, 0, 30);
    hV0R[i]->Sumw2();

    hV0PVDCA[i] = new TH1D(Form("hV0PVDCA_%s", SignalTypeName[i]), "; DCA_{V0PV}",
          50, 0, 5);
    hV0PVDCA[i]->Sumw2();

    hBachDCA[i] = new TH1D(Form("hBachDCA_%s", SignalTypeName[i]), "; DCA_{Bach}",
          50, 0, 5);
    hBachDCA[i]->Sumw2();

    hCascROverPt[i] = new TH1D(Form("hCascROverPt_%s", SignalTypeName[i]), "; R_{casc}/P_{t}^{casc}",
          60, 0, 30);
    hCascROverPt[i]->Sumw2();

    hV0ROverPt[i] = new TH1D(Form("hV0ROverPt_%s", SignalTypeName[i]), "; R_{V0}/P_{t}^{casc}",
          100, 0, 50);
    hV0ROverPt[i]->Sumw2();

    hCascV0DCA[i] = new TH1D(Form("hCascV0DCA_%s", SignalTypeName[i]), "; DCA_{CascV0}",
          30, 0, 15);
    hCascV0DCA[i]->Sumw2();

    hV0BachDCA[i] = new TH1D(Form("hV0BachDCA_%s", SignalTypeName[i]), "; DCA_{V0Bach}",
          30, 0, 15);
    hV0BachDCA[i]->Sumw2();

    CutList->Add(hCascPA[i]);
    CutList->Add(hcascR[i]);
    CutList->Add(hPosDCA[i]);
    CutList->Add(hNegDCA[i]);
    CutList->Add(hCascPVDCA[i]);
    CutList->Add(hV0DaughtersDCA[i]);
    CutList->Add(hV0PA[i]);
    CutList->Add(hV0R[i]);
    CutList->Add(hV0PVDCA[i]);
    CutList->Add(hBachDCA[i]);
    CutList->Add(hCascROverPt[i]);
    CutList->Add(hV0ROverPt[i]);
    CutList->Add(hCascV0DCA[i]);
    CutList->Add(hV0BachDCA[i]);
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
  outFile->Close();
}

void analyze_cuts(const Char_t* inFileName,
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
      if(cascade->GetEtaCasc() > 0.8) // needed to check
        continue;

      const Double_t dMassOmega     = cascade->GetIMO() - massOmega;
      const Double_t dMassXi = cascade->GetIMXi() - massXi;

      if(!CheckCascOmegaToXiMass(cascade))
        continue;

      hCascStat->Fill("#Xi check", 1);

      const Double_t deltaM     = cascade->GetIMO() - massOmega;
      if(TMath::Abs(deltaM) > 0.029999999)
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

      FillCutHists(2, cascade, event);
      if(isMC){
        bIsRealOmegaCascade = IsRealOmegaCascade(cascade);
        if(bIsRealOmegaCascade){
          FillCutHists(0, cascade, event);
        }
        else{
          FillCutHists(1, cascade, event);
        }
      }

      bPassedLoose = CheckCascLooseCuts(cascade);
      if(!bPassedLoose)
        continue;

      FillCutHists(3, cascade, event);
      if(isMC){
        bIsRealOmegaCascade = IsRealOmegaCascade(cascade);
        if(bIsRealOmegaCascade){
          FillCutHists(5, cascade, event);
        }
        else{
          FillCutHists(4, cascade, event);
        }
      }

      hCascStat->Fill("Loose", 1);

      bPassedStandard = CheckCascStandardCuts(cascade);
      if(bPassedStandard){
        hCascStat->Fill("Standard", 1);
        FillCutHists(6, cascade, event);
        if(isMC){
          bIsRealOmegaCascade = IsRealOmegaCascade(cascade);
          if(bIsRealOmegaCascade){
            FillCutHists(8, cascade, event);
          }
          else{
            FillCutHists(7, cascade, event);
          }
        }
        Int_t bin = 1;
        if(cascade->GetCharge() < 0){bin = 0;}
          hOmegaInvMassVsPt[bin]->Fill(cascade->GetPtCasc(), dMassOmega, event->GetV0Mmultiplicity());
        hOmegaInconsistencyXi->Fill(dMassXi, dMassOmega);
        // MC closure for signal inside reconstructed
        if(isMC){
          if(IsRealOmegaCascade(cascade)){
            if(cascade->GetCharge() < 0){bin = 0;}
              hOmegaInvMassVsPtTrue[bin]->Fill(cascade->GetPtCasc(), dMassOmega, event->GetV0Mmultiplicity());                  
          }
        }
      }

      bPassedTight = CheckCascTightCuts(cascade);
      if(bPassedTight){
        hCascStat->Fill("Tight", 1);
        FillCutHists(9, cascade, event);
        if(isMC){
          bIsRealOmegaCascade = IsRealOmegaCascade(cascade);
          if(bIsRealOmegaCascade){
            FillCutHists(11, cascade, event);
          }
          else{
            FillCutHists(10, cascade, event);
          }
        }
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
