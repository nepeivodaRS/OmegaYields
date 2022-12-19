#include "analyze_tree_NR.h"

/*
  To run code:
  ============

  aliroot -l

  .L analyze_tree_NR.C+

  analyze_tree_MC("./outputMadeTrees/mc_tree_pp17j.root", "./outputAnalTrees/mc_analTree_pp17j.root", 0, 1)

  analyze_tree_MC("./outputMadeTrees/mc_tree_all.dat", "./outputAnalTrees/mc_analTree_all.root", 0, 1)
  analyze_tree_MC("./outputMadeTreesData/tree_all.dat", "./outputAnalTrees/analTreeData_all.root", 0, 0)
*/

void InitHists(){
  // Initialization of histograms
  for(Int_t i = 0; i <= nMinvBins; i++) {
    minvBins[i] = -0.03 + i * 0.06/nMinvBins;
  }

  PtHistList = new TList();

  hOmegaMB = new TH1D("hOmegaMB", "; p_{T} [GeV/c]",
       nPtBinsMB, xBinsMB);
  hOmegaMB->Sumw2();

  hOmegaHM = new TH1D("hOmegaHM", "; p_{T} [GeV/c]",
       nPtBinsHM, xBinsHM);
  hOmegaHM->Sumw2();
  
  hOmegaVHM = new TH1D("hOmegaVHM", "; p_{T} [GeV/c]",
        nPtBinsHM, xBinsHM);
  hOmegaVHM->Sumw2();

  hGenOmegaMB = new TH1D("hGenOmegaMB", ";  p_{T} [GeV/c]",
        nPtBinsMB, xBinsMB);
  hGenOmegaMB->Sumw2();

  hGenOmegaHM = new TH1D("hGenOmegaHM", ";  p_{T} [GeV/c]",
          nPtBinsHM, xBinsHM);
  hGenOmegaHM->Sumw2();

  hGenOmegaVHM = new TH1D("hGenOmegaVHM", ";  p_{T} [GeV/c]",
          nPtBinsHM, xBinsHM);
  hGenOmegaVHM->Sumw2();

  PtHistList->Add(hOmegaMB);
  PtHistList->Add(hOmegaHM);
  PtHistList->Add(hOmegaVHM);

  PtHistList->Add(hGenOmegaMB);
  PtHistList->Add(hGenOmegaHM);
  PtHistList->Add(hGenOmegaVHM);

  InvMassList = new TList();

  for(Int_t i = 0; i < nSignalTypes; i++) {
    hOmegaInvMassVsPt[i] = new TH3D(Form("hOmegaInvMassVsPt_%s", SignalTypeName[i]), "M_{inv} vs p_{T}; p_{T}^{casc} [GeV/c]; M_{inv} - M_{#Omega} [GeV/c^{2}]; Centrality V0M",
            nPtBins, xBins, nMinvBins, minvBins, nCentrBins, xCentrBins);
    hOmegaInvMassVsPt[i]->Sumw2();

    hLambdaInvMassVsPt[i] = new TH3D(Form("hLambdaInvMassVsPt_%s", SignalTypeName[i]), "M_{inv} vs p_{T}; p_{T}^{V0} [GeV/c]; M_{inv} - M_{#Lambda} [GeV/c^{2}]; Centrality V0M",
            nPtBins, xBins, nMinvBins, minvBins, nCentrBins, xCentrBins);
    hLambdaInvMassVsPt[i]->Sumw2();

    InvMassList->Add(hOmegaInvMassVsPt[i]);
    InvMassList->Add(hLambdaInvMassVsPt[i]);
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

  outFile->Close();
}

void analyze_tree_MC(const Char_t* inFileName,
		const Char_t* outFileName,
		const Int_t maxEvents = 0,
    Bool_t isMC = kFALSE)
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
        hGenOmegaMB->Fill(trackMC->GetPt());
        if(HMevent)
          hGenOmegaHM->Fill(trackMC->GetPt());
        if(VHMevent)
          hGenOmegaVHM->Fill(trackMC->GetPt());
      }
    }
/////////////////
    const Int_t nCascades = allCascades->GetEntries();

    for(Int_t i = 0; i < nCascades; i++) {

      Bool_t bIsRealOmegaCascade = kFALSE;
      Bool_t bPassedLoose = kFALSE;
      Bool_t bPassedStandard = kFALSE;
      Bool_t bPassedTight = kFALSE;

      hCascStat->Fill("Total", 1);

      AliAnalysisPIDCascade* cascade = (AliAnalysisPIDCascade*)allCascades->At(i);
      if(!CheckCascOmegaToXiMass(cascade))
        continue;

      hCascStat->Fill("#Xi check", 1);

      const Double_t deltaM     = cascade->GetIMO() - massOmega;
      if(TMath::Abs(deltaM) > 0.1)
        continue;

      hCascStat->Fill("M_{#Omega} check", 1);

      Double_t weight = 0;
      if(TMath::Abs(deltaM) < 0.01)
        weight = 1;
      else if(TMath::Abs(deltaM) < 0.02)
        weight = -1;
      else
        weight = 0; // weight = 0;

      FillCutHists(2, cascade, event); // Fill for all cascades

      if(isMC){
        bIsRealOmegaCascade = IsRealOmegaCascade(cascade);
        if(bIsRealOmegaCascade){
          FillCutHists(0, cascade, event); // real omega cascade
          hCascStat->Fill("Gen #Omega", 1);
        }
        else{
          FillCutHists(1, cascade, event); // background cascades
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
      FillCutHists(3, cascade, event); // Passed Loose cuts
      if(!bIsRealOmegaCascade){
        FillCutHists(4, cascade, event); // Passed Loose cuts but fake
      }
      else{
        FillCutHists(5, cascade, event); // Passed Loose cuts and real
      }

      bPassedStandard = CheckCascStandardCuts(cascade);

      if(bPassedStandard){
        hCascStat->Fill("Standart", 1);
        FillCutHists(6, cascade, event); // Passed Standard cuts
        if(!bIsRealOmegaCascade){
          FillCutHists(7, cascade, event); // Passed Standard cuts but fake
        }
        else{
          FillCutHists(8, cascade, event); // Passed Standard cuts and real
        }
        // Pt histograms
        hOmegaMB->Fill(cascade->GetPtCasc(), weight);
        if(HMevent == kTRUE)
          hOmegaHM->Fill(cascade->GetPtCasc(), weight);
        if(VHMevent == kTRUE)
          hOmegaVHM->Fill(cascade->GetPtCasc(), weight);
      }

      bPassedTight = CheckCascTightCuts(cascade);

      if(bPassedTight){
        hCascStat->Fill("Tight", 1);
        FillCutHists(9, cascade, event); // Passed Tight cuts
        if(!bIsRealOmegaCascade){
          FillCutHists(10, cascade, event); // Passed Tight cuts but fake
        }
        else{
          FillCutHists(11, cascade, event); // Passed Tight cuts and real
        }
      }
    }
  }

  // Rapidity correction
  TF1* fRap = new TF1("fRap", rap_correction, 0.0, 50.0, 2);
  fRap->SetParameters(0.8, massOmega);
  
  hOmegaMBmc->Scale(1./nMB);
  NormalizeHistogram(hOmegaMBmc);
  hOmegaMBmc->Divide(fRap);
  WriteToFile(outFile);
}
