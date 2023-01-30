#include "analyze_cuts_topo.h"

/*
  To run code:
  ============

  aliroot -l

  .L analyze_cuts_topo.C

  analyze_cuts_topo("./outputCuts/mc_cuts_27.root", "./outputCuts/mc_cuts_results_27.root")
*/

void DrawCutLines(Int_t nAxis, Double_t Height){
  const Int_t n = 2;
  Double_t xL[n], yL[n];
  Double_t xS[n], yS[n];
  Double_t xT[n], yT[n];
  yL[0] = yS[0] = yT[0] = 0;
  yL[1] = yS[1] = yT[1] = Height + 100;
  xL[0] = xL[1] = 0;
  xS[0] = xS[1] = 0;
  xT[0] = xT[1] = 0;

  switch (nAxis) {
    case 0:
    {
      xL[0] = xL[1] = 0.97;
      xS[0] = xS[1] = 0.98;
      xT[0] = xT[1] = 0.99;
      break;
    }
    case 1:
    {
      xL[0] = xL[1] = 0.6;
      xS[0] = xS[1] = 0.6;
      xT[0] = xT[1] = 0.6;
      break;
    }
    case 2:
    {
      xL[0] = xL[1] = 0.02;
      xS[0] = xS[1] = 0.03;
      xT[0] = xT[1] = 0.04;
      break;
    }
    case 3:
    {
      xL[0] = xL[1] = 0.02;
      xS[0] = xS[1] = 0.03;
      xT[0] = xT[1] = 0.04;
      break;
    }
    case 4:
    {
      // xL[0] = xL[1] = 1.0;
      // xS[0] = xS[1] = 1.0;
      // xT[0] = xT[1] = 1.0;
      break;
    }
    case 5:
    {
      xL[0] = xL[1] = 1.7;
      xS[0] = xS[1] = 1.6;
      xT[0] = xT[1] = 1.5;
      break;
    }
  case 6:
    {
      xL[0] = xL[1] = 0.97;
      xS[0] = xS[1] = 0.98;
      xT[0] = xT[1] = 0.99;
      break;
    }
  case 7:
    {
      xL[0] = xL[1] = 1.3;
      xS[0] = xS[1] = 1.4;
      xT[0] = xT[1] = 1.5;
      break;
    }
  case 8:
    {
      xL[0] = xL[1] = 0.06;
      xS[0] = xS[1] = 0.07;
      xT[0] = xT[1] = 0.08;
      break;
    }
  case 9:
    {
      xL[0] = xL[1] = 0.05;
      xS[0] = xS[1] = 0.05;
      xT[0] = xT[1] = 0.05;
      break;
    }
  case 10:
    {
      xL[0] = xL[1] = 16;
      xS[0] = xS[1] = 15;
      xT[0] = xT[1] = 14;
      break;
    }
  case 11:
    {
      xL[0] = xL[1] = 45;
      xS[0] = xS[1] = 40;
      xT[0] = xT[1] = 35;
      break;
    }
  case 12:
    {
      xL[0] = xL[1] = 11;
      xS[0] = xS[1] = 10;
      xT[0] = xT[1] = 9;
      break;
    }
  case 13: // hV0BachDCA
    {
      xL[0] = xL[1] = 1;
      xS[0] = xS[1] = 2;
      xT[0] = xT[1] = 3;
      break;
    }
  default:
    {
      break;
    }
  }
  auto grL = new TGraph(n,xL,yL);
  grL->SetLineColor(8);
  grL->SetLineWidth(2);
  grL->SetLineStyle(2);
  grL->Draw("SAME");

  auto grS = new TGraph(n,xS,yS);
  grS->SetLineColor(9);
  grS->SetLineWidth(2);
  grS->SetLineStyle(2);
  grS->Draw("SAME");

  auto grT = new TGraph(n,xT,yT);
  grT->SetLineColor(2);
  grT->SetLineWidth(2);
  grT->SetLineStyle(2);
  grT->Draw("SAME");
}

void AddRatioOfHistsPlotToList(TH1D *h1, TH1D *h2, TList* InList, const Char_t* CanvName, const Char_t* axisName){
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  h1->GetYaxis()->SetTitle("Entries");
  h1->SetLineColor(kMagenta+1);
  h1->SetLineWidth(1);
  h2->SetLineWidth(1);
  Double_t MaxHeightHist = (h1->GetMaximum() > h2->GetMaximum()) ? h1->GetMaximum() : h2->GetMaximum();

  auto canvas= new TCanvas(CanvName, CanvName, 1280, 720);
  auto rp = new TRatioPlot(h1, h2);
  rp->SetH1DrawOpt("P");
  rp->SetH2DrawOpt("P");
  rp->SetSeparationMargin(0.0);
  rp->Draw();
  rp->GetLowerRefYaxis()->SetTitle(axisName);
  rp->GetUpperPad()->SetLogy();
  rp->GetLowerPad()->SetFillColorAlpha(0., 0.);
  rp->GetLowerPad()->SetFillColorAlpha(0., 0.);
  rp->GetUpperRefYaxis()->SetLimits(0, MaxHeightHist);
  rp->GetUpperPad()->cd();
  TLegend*  lLegend     = new TLegend(0.22,0.66,0.42,0.86);
  lLegend->SetLineColorAlpha(0.,0.);
  lLegend->SetFillColorAlpha(0.,0.);
  lLegend->AddEntry  (h1,h1->GetName(),"lpf");
  lLegend->AddEntry  (h2,h2->GetName(),"lpf");
  lLegend->Draw();

  // Draw cut lines
  //std::cout << h1->GetXaxis()->GetTitle() << std::endl;
  for(Int_t i = 0; i < nXAxisNames; i++) {
    if(strcmp(XAxisNames[i], h1->GetXaxis()->GetTitle()) == 0){
      //std::cout << h1->GetXaxis()->GetTitle() << std::endl;
      DrawCutLines(i, MaxHeightHist);
      break;
    }
  }

  //uLatex->SetTextSize(0.07);
  uLatex->DrawLatexNDC(0.23, 0.86,"pp: 16k, 17hjklmor, 18bdefglmnop");
  canvas->Update();
  InList->Add(canvas, CanvName);
  gROOT->SetBatch(kFALSE);
}

void CalcRatios(Int_t Sig1, Int_t Sig2, TList* InList, const Char_t* axisName, const Char_t* ListId){
  AddRatioOfHistsPlotToList(hCascPA[Sig1], hCascPA[Sig2], InList, Form("CascPA_%s", ListId), axisName);
  AddRatioOfHistsPlotToList(hcascR[Sig1], hcascR[Sig2], InList, Form("CascR_%s", ListId), axisName);
  AddRatioOfHistsPlotToList(hCascPVDCA[Sig1], hCascPVDCA[Sig2], InList, Form("CascPVDCA_%s", ListId), axisName);
  AddRatioOfHistsPlotToList(hCascROverPt[Sig1], hCascROverPt[Sig2], InList, Form("CascROverPt_%s", ListId), axisName);
  AddRatioOfHistsPlotToList(hCascV0DCA[Sig1], hCascV0DCA[Sig2], InList, Form("CascV0DCA_%s", ListId), axisName);
  AddRatioOfHistsPlotToList(hPosDCA[Sig1], hPosDCA[Sig2], InList, Form("PosDCA_%s", ListId), axisName);
  AddRatioOfHistsPlotToList(hNegDCA[Sig1], hNegDCA[Sig2], InList, Form("NegDCA_%s", ListId), axisName);
  AddRatioOfHistsPlotToList(hV0DaughtersDCA[Sig1], hV0DaughtersDCA[Sig2], InList, Form("V0DaughtersDCA_%s", ListId), axisName);
  AddRatioOfHistsPlotToList(hV0PA[Sig1], hV0PA[Sig2], InList, Form("V0PA_%s", ListId), axisName);
  AddRatioOfHistsPlotToList(hV0R[Sig1], hV0R[Sig2], InList, Form("V0R_%s", ListId), axisName);
  AddRatioOfHistsPlotToList(hV0PVDCA[Sig1], hV0PVDCA[Sig2], InList, Form("V0PVDCA_%s", ListId), axisName);
  AddRatioOfHistsPlotToList(hV0ROverPt[Sig1], hV0ROverPt[Sig2], InList, Form("V0ROverPt_%s", ListId), axisName);
  AddRatioOfHistsPlotToList(hV0BachDCA[Sig1], hV0BachDCA[Sig2], InList, Form("hV0BachDCA_%s", ListId), axisName);
}

void WriteToFile(TFile* outF){
  // Write to file
  outF->cd();
  outF->WriteObject(ListFakeToRecRatiosTight, "ListFakeToRecRatiosTight", "SingleKey");
  outF->WriteObject(ListRecToGenRatiosStandart, "ListRecToGenRatiosStandart", "SingleKey");
  outF->WriteObject(ListFakeToAll, "ListFakeToAll", "SingleKey");
  outF->WriteObject(ListRealToFakeStandard, "ListRealToFakeStandard", "SingleKey");
  outF->WriteObject(ListGenToRecTight, "ListGenToRecTight", "SingleKey");
  outF->WriteObject(ListRealToFakeTight, "ListRealToFakeTight", "SingleKey");
  //outF->Close();
}

void analyze_cuts_topo(const Char_t* inFileName, const Char_t* outFileName){
  //Open the input file and set up the classes.
  inFile = TFile::Open(inFileName);
  if(!inFile)
    return;

  outFile = new TFile(outFileName, "RECREATE");

  ListOfCutHists = (TList*)inFile->Get("ListOfCuts");
  for(Int_t i = 0; i < nSignalTypes; i++) {
    hCascPA[i] = (TH1D*)ListOfCutHists->FindObject(Form("hCascPA_%s", SignalTypeName[i]));
    hcascR[i] = (TH1D*)ListOfCutHists->FindObject(Form("hcascR_%s", SignalTypeName[i]));
    hCascPVDCA[i] = (TH1D*)ListOfCutHists->FindObject(Form("hCascPVDCA_%s", SignalTypeName[i]));
    hCascROverPt[i] = (TH1D*)ListOfCutHists->FindObject(Form("hCascROverPt_%s", SignalTypeName[i]));
    hCascV0DCA[i] = (TH1D*)ListOfCutHists->FindObject(Form("hCascV0DCA_%s", SignalTypeName[i]));
    hPosDCA[i] = (TH1D*)ListOfCutHists->FindObject(Form("hPosDCA_%s", SignalTypeName[i]));
    hNegDCA[i] = (TH1D*)ListOfCutHists->FindObject(Form("hNegDCA_%s", SignalTypeName[i]));
    hV0DaughtersDCA[i] = (TH1D*)ListOfCutHists->FindObject(Form("hV0DaughtersDCA_%s", SignalTypeName[i]));
    hV0PA[i] = (TH1D*)ListOfCutHists->FindObject(Form("hV0PA_%s", SignalTypeName[i]));
    hV0R[i] = (TH1D*)ListOfCutHists->FindObject(Form("hV0R_%s", SignalTypeName[i]));
    hV0PVDCA[i] = (TH1D*)ListOfCutHists->FindObject(Form("hV0PVDCA_%s", SignalTypeName[i]));
    hV0ROverPt[i] = (TH1D*)ListOfCutHists->FindObject(Form("hV0ROverPt_%s", SignalTypeName[i]));
    hBachDCA[i] = (TH1D*)ListOfCutHists->FindObject(Form("hBachDCA_%s", SignalTypeName[i]));
    hV0BachDCA[i] = (TH1D*)ListOfCutHists->FindObject(Form("hV0BachDCA_%s", SignalTypeName[i]));
  }

  // Setup Latex var
  uLatex = new TLatex();
  uLatex->SetTextFont(60);

  // Ratio Hists
  ListFakeToRecRatiosTight = new TList();
  CalcRatios(10, 9, ListFakeToRecRatiosTight, "Fake_{tt}/Rec_{tt}", "FoR(tight)");

  // ListRecToGenRatiosStandart = new TList();
  // CalcRatios(6, 0, ListRecToGenRatiosStandart, "Rec_{sd}/Gen", "RoG(st)");

  // ListFakeToAll = new TList();
  // CalcRatios(1, 2, ListFakeToAll, "Background/All", "NoCuts");

  // ListRealToFakeStandard = new TList();
  // CalcRatios(8, 7, ListRealToFakeStandard, "Gen_{sd}/Fake_{sd}", "GoF(st)");

  // ListGenToRecTight = new TList();
  // CalcRatios(0, 9, ListGenToRecTight, "Gen/Rec_{tt}", "GoR(tt)");

  // ListRealToFakeTight = new TList();
  // CalcRatios(11, 10, ListRealToFakeTight, "Gen_{tt}/Fake_{tt}", "GoF(tt)");

  WriteToFile(outFile);
}