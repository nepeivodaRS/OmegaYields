#include "make_tree_omega_MC.h"
/*
  To run code:
  ============

  aliroot -l

  .L make_tree_omega.C+
  
  // Data

  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp16k.dat", "./outputMadeTreesData/tree_pp16k.root", 0)
  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp17h.dat", "./outputMadeTreesData/tree_pp17h.root", 0)
  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp17i.dat", "./outputMadeTreesData/tree_pp17i.root", 0)
  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp17j.dat", "./outputMadeTreesData/tree_pp17j.root", 0)
  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp17k.dat", "./outputMadeTreesData/tree_pp17k.root", 0)
  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp17l.dat", "./outputMadeTreesData/tree_pp17l.root", 0)
  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp17m.dat", "./outputMadeTreesData/tree_pp17m.root", 0)
  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp17o.dat", "./outputMadeTreesData/tree_pp17o.root", 0)
  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp17r.dat", "./outputMadeTreesData/tree_pp17r.root", 0)

  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp18b.dat", "./outputMadeTreesData/tree_pp18b.root", 0)
  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp18d.dat", "./outputMadeTreesData/tree_pp18d.root", 0)
  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp18e.dat", "./outputMadeTreesData/tree_pp18e.root", 0)
  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp18f.dat", "./outputMadeTreesData/tree_pp18f.root", 0)
  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp18g.dat", "./outputMadeTreesData/tree_pp18g.root", 0)
  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp18l.dat", "./outputMadeTreesData/tree_pp18l.root", 0)
  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp18m.dat", "./outputMadeTreesData/tree_pp18m.root", 0)
  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp18n.dat", "./outputMadeTreesData/tree_pp18n.root", 0)
  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp18o.dat", "./outputMadeTreesData/tree_pp18o.root", 0)
  make_tree_omega("/home/pchristi/work/analysis/pp_vs_mult/xi_analysis/tree/pp18p.dat", "./outputMadeTreesData/tree_pp18p.root", 0)
  
  // MC
  make_tree_omega("./inputDatLists/mc_pp17h.dat", "./outputTreesMC/mc_tree_pp17h.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp17i.dat", "./outputTreesMC/mc_tree_pp17i.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp16k.dat", "./outputTreesMC/mc_tree_pp16k.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp17j.dat", "./outputTreesMC/mc_tree_pp17j.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp17k.dat", "./outputTreesMC/mc_tree_pp17k.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp17l.dat", "./outputTreesMC/mc_tree_pp17l.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp17m.dat", "./outputTreesMC/mc_tree_pp17m.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp17o.dat", "./outputTreesMC/mc_tree_pp17o.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp17r.dat", "./outputTreesMC/mc_tree_pp17r.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp18b.dat", "./outputTreesMC/mc_tree_pp18b.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp18d.dat", "./outputTreesMC/mc_tree_pp18d.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp18e.dat", "./outputTreesMC/mc_tree_pp18e.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp18f.dat", "./outputTreesMC/mc_tree_pp18f.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp18g.dat", "./outputTreesMC/mc_tree_pp18g.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp18l.dat", "./outputTreesMC/mc_tree_pp18l.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp18m.dat", "./outputTreesMC/mc_tree_pp18m.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp18n.dat", "./outputTreesMC/mc_tree_pp18n.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp18o.dat", "./outputTreesMC/mc_tree_pp18o.root", 0, 1)
  make_tree_omega("./inputDatLists/mc_pp18p.dat", "./outputTreesMC/mc_tree_pp18p.root", 0, 1)
  */

void InitHists(){
  for(Int_t i = 0; i <= nMinvBins; i++) {
    minvBins[i] = -0.03 + i * 0.06/nMinvBins;
  }
  
  // Xchecks histograms

  for(Int_t i = 0; i < nSpecies; i++) {
    hInvMassOmegaVsPt[i] = new TH2D(Form("hInvMassVsPt%s", pdgNameOmega[i]),
             Form("M_{inv} vs p_{T} for %s generated; p_{T} [GeV/c]; M_{inv} - M_{#Omega} [GeV/c^{2}]", pdgNameOmega[i]),
             nPtBins, xBins, nMinvBins, minvBins);
    hInvMassOmegaVsPt[i]->Sumw2();

    // Pt histograms

    hGenOmegaMB[i] = new TH1D(Form("hGenOmegaMB_%s", pdgNameOmega[i]), ";  p_{T} [GeV/c]",
          nPtBinsMB, xBinsMB);
    hGenOmegaMB[i]->Sumw2();

    hGenOmegaHM[i] = new TH1D(Form("hGenOmegaHM_%s", pdgNameOmega[i]), ";  p_{T} [GeV/c]",
            nPtBinsHM, xBinsHM);
    hGenOmegaHM[i]->Sumw2();

    hGenOmegaVHM[i] = new TH1D(Form("hGenOmegaVHM_%s", pdgNameOmega[i]), ";  p_{T} [GeV/c]",
            nPtBinsHM, xBinsHM);
    hGenOmegaVHM[i]->Sumw2();
  }
  
  hVtxStatus = new TH1D("hVtxStatus", "Vtx status - No Vtx = -1, Vtx outside cut = 0, Vtx inside = 1",
            3, -1.5, 1.5);
  hVtxStatus->Sumw2();
  
  hVtxZ = new TH1D("hVtxZ", "Accepted vtx z; z [cm]; Counts",
       90, -15, 15);
  hVtxZ->Sumw2();

  hEventStat = new TH1I("hEventStat","",3,0,3);
  
  // Output trees

  treeOut = new TTree("tree", "debug tree");

  treeOut->Branch("event", &event);

  mcArrayOut = new TClonesArray("AliAnalysisPIDCascadeParticle", 1000);
  treeOut->Bronch("generatedOmega", "TClonesArray", &mcArrayOut); // Tracks of generated omegas

  cascadeArrayOut = new TClonesArray("AliAnalysisPIDCascade", 1000);
  treeOut->Bronch("allCascades", "TClonesArray", &cascadeArrayOut); // All cascades
}

void make_tree_omega(const Char_t* inFileName,
		  const Char_t* outFileName,
		  const Int_t maxEvents = 0,
      Bool_t isMC = kFALSE)
{
  // Open the input file and set up the classes.
  if(strstr(inFileName, ".dat")) {
    tree = ReadChainFromFile(inFileName, "PIDTree", 0);
  } else {
    TFile* inFile = TFile::Open(inFileName);
    if(!inFile){
      return;
    }
    tree = (TTree*)inFile->Get("PIDTree");
  }

  TFile* outFile = new TFile(outFileName, "RECREATE");
  outFile->cd();
  
  tree->SetBranchStatus("AnalysisV0Track*",0);
  tree->SetBranchAddress("AnalysisEvent", &event);
  tree->SetBranchAddress("AnalysisTrack", &trackArray);
  tree->SetBranchAddress("AnalysisCascadeTrack", &cascadeArray);
  if(isMC)
    tree->SetBranchAddress("AnalysisParticle", &mcTrackArray);

  InitHists();

  // Loop over events

  Int_t nEvents = tree->GetEntries();
  cout << "Number of events: " << nEvents << endl;
  
  // Reduce number of events
  if(maxEvents>0 && maxEvents < nEvents) {
    nEvents = maxEvents;
    cout << "N events was reduced to: " << maxEvents << endl;
  }

  // Main loop over events
  for(Int_t n = 0; n < nEvents; n++) {
    // Cleaning
    mcArrayOut->Clear();
    cascadeArrayOut->Clear();

    tree->GetEntry(n);
    
    if((n+1)%100000==0){cout << "Event: " << n+1 << "/" << nEvents << endl;}

    // Get vertex status
    Int_t vtxStatus = AcceptEvent(event, trackArray);
    
    if(vtxStatus == -2) // reject
      continue;
    
    hVtxStatus->Fill(vtxStatus);
    
    if(vtxStatus < 1) // event was rejected
      continue;

    // Fill good vertxs
    hVtxZ->Fill(event->GetVertexZ());


    Bool_t HMevent = kFALSE;
    Bool_t VHMevent = kFALSE;

    // increase event copunters
    hEventStat->Fill("MB", 1);

    if(event->GetV0Mmultiplicity() < 10){
      HMevent = kTRUE;
      hEventStat->Fill("HM", 1);
    }
    if(event->GetV0Mmultiplicity() < 1){
      VHMevent = kTRUE;
      hEventStat->Fill("VHM", 1);
    }

    // Analyze MC generated data
    Int_t nOmegaMc = 0;
    if(isMC){
      const Int_t nMcParticles = mcTrackArray->GetEntries();
      for(Int_t i = 0; i < nMcParticles; i++) {
        AliAnalysisPIDCascadeParticle* trackMC = (AliAnalysisPIDCascadeParticle*)mcTrackArray->At(i);

        // if not omega then reject
        if(TMath::Abs(trackMC->GetPdgCode()) != 3334)
        continue;

        // Reject particles that are not primary
        if(!trackMC->GetPrimaryStatus())
        continue;

        if(TMath::Abs(trackMC->GetEta()) > 0.8)
        continue;

        if(TMath::Abs(trackMC->GetPt()) < 0.15)
        continue;

        Int_t bin = 1;
        if(trackMC->GetSign() < 0){bin = 0;}

        hGenOmegaMB[bin]->Fill(trackMC->GetPt());
        if(HMevent)
          hGenOmegaHM[bin]->Fill(trackMC->GetPt());
        if(VHMevent)
          hGenOmegaVHM[bin]->Fill(trackMC->GetPt());

        new((*mcArrayOut)[nOmegaMc]) AliAnalysisPIDCascadeParticle(*trackMC);
        nOmegaMc++;
      }
    }

    // Analyze cascades
    Int_t nCascOut = 0;
    Int_t nAllCascOut = 0;
    
    const Int_t nCascades = cascadeArray->GetEntries();

    for(Int_t i = 0; i < nCascades; i++) {
      AliAnalysisPIDCascade* cascade = (AliAnalysisPIDCascade*)cascadeArray->At(i);
      AliAnalysisPIDCascadeV0* v0 = cascade->GetV0();
      AliAnalysisPIDCascadeTrack* bachelor = cascade->GetBachAnalysisTrack();

      new((*cascadeArrayOut)[nAllCascOut]) AliAnalysisPIDCascade(*cascade);
      nAllCascOut++;

      // Inv mass of cascades
      const Double_t deltaM     = cascade->GetIMO() - massOmega;

      Int_t bin = 1;
      if(cascade->GetCharge() < 0){bin = 0;}

      hInvMassOmegaVsPt[bin]->Fill(cascade->GetPtCasc(), deltaM);
    }

    Bool_t fillTree = kTRUE;
    if(fillTree){
      treeOut->Fill();
    }
  }
  outFile->Write();
  outFile->Close();
}