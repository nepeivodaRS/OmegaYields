// #include "helpFunc.C"
#include "ptDistr.h"

/*
  .L ptDistr.C+

  make_results("./outputAnal/mc_anal_testBinning.root", "./outputEff/mc_Eff_testBinning.root", "./outputPtHists/PtHist.root", 1)
  
 */

TFile* FindFileFresh(const Char_t* fileName)
{
  // Find file
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName);
  if(file) {
    file->Close();
    delete file;
  }

  file = TFile::Open(fileName, "READ");

  if(!file)
    cout << "File : " << fileName << " was not found" << endl;

  return file;
}

Double_t rap_correction(Double_t* x, Double_t* par)
{
  Double_t pt = x[0];  

  Double_t eta  = par[0];
  Double_t mass = par[1];

  const Double_t mt = TMath::Sqrt(pt*pt + mass*mass);
  
  const Double_t rap = TMath::ASinH(pt/mt*TMath::SinH(eta));
  
  return rap/eta;
}

void NormalizeHistogram(TH1D* hist){
  const Int_t nBins = hist->GetXaxis()->GetNbins();
  for(Int_t bin = 1; bin <= nBins; bin++) {
    const Double_t binwidth = hist->GetXaxis()->GetBinWidth(bin);
    Double_t scale = 1.0/binwidth/1.6;
    hist->SetBinContent(bin, scale*hist->GetBinContent(bin));
    hist->SetBinError(bin, scale*hist->GetBinError(bin));
  }
}

Double_t DoubleSidedCB2(Double_t x, Double_t mu, Double_t width, Double_t a1, Double_t p1, Double_t a2, Double_t p2)
{
  Double_t u   = (x-mu)/width;
  Double_t A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
  Double_t A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
  Double_t B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
  Double_t B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

  Double_t result(1);
  if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2);
  else            result *= A2*TMath::Power(B2+u,-p2);
  return result;
}

Double_t DoubleSidedCB(Double_t* x, Double_t *par)
{
  return(par[0] * DoubleSidedCB2(x[0], par[1],par[2],par[3],par[4],par[5],par[6]));
}

Double_t background(Double_t *x, Double_t *par)  {
return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

Double_t GAUSS(Double_t *x, Double_t *par){
return par[0]*exp(-0.5*(((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2])));
}

Double_t fitFunctionG(Double_t *x, Double_t *par) {
  return background(x,par) + GAUSS(x,&par[3]);
}

Double_t fitFunctionCB(Double_t *x, Double_t *par) {
  return background(x,par) + DoubleSidedCB(x,&par[3]);
}

void SignalExtractionPt(const Double_t *xPtBins, const Int_t nPtBins, TH3D *inHist3D, Int_t leftCentr, Int_t rightCentr, TH1D* inHist, TFile* outFile){
  gStyle->SetOptStat("me");
  // Create dir to store all the fitted hists
  TDirectory* dir = outFile->mkdir(Form("PtFitHists_%s_from_%d_to_%d_centr", inHist3D->GetName(), leftCentr, rightCentr));
  dir->cd();
  // Clone hist not to change the original one
  TH3D* invMassHist = (TH3D*)inHist3D->Clone();
  invMassHist->Write();
  invMassHist->GetZaxis()->SetRange(leftCentr, rightCentr);
  TH2D* hProfileInvMassZ = static_cast<TH2D*>(invMassHist->Project3D("xy"));
  hProfileInvMassZ->Write();
  for(Int_t i = 0; i < nPtBins; ++i) {
    gROOT->SetBatch(kFALSE);
    TH1D* hProfileInvMassX = hProfileInvMassZ->ProjectionX("_px", hProfileInvMassZ->GetYaxis()->FindBin(xPtBins[i] + 0.00001), hProfileInvMassZ->GetYaxis()->FindBin(xPtBins[i+1]-0.00001));
    hProfileInvMassX->SetTitle(Form("M_{inv} in pt bins fit, %d bin",  i+1));
    TF1 *fitFcn = new TF1("fitFcn",fitFunctionG,-0.03, 0.03,6);
    Double_t sigPickApprox = 20;
    fitFcn->SetParameters(1,1,1,sigPickApprox,0,0.001);
    fitFcn->SetNpx(1e4);
    gROOT->SetBatch(kTRUE);
    // Create canvas for fitted inv mass
    TCanvas *c1 = new TCanvas(Form("PtFit_%d", i+1),"PtHist",10,10,1200,900);
    c1->cd();
    TFitResultPtr FitResult = hProfileInvMassX->Fit("fitFcn","eprS");
    Double_t par[6];
    fitFcn->GetParameters(par);
    // Fill histogram with fitted signal: integral divided by the bin width. And eval the error of each bin
    Double_t FillValue = par[3]*TMath::Power(TMath::TwoPi(), 0.5)*TMath::Abs(par[5])/(0.06/nMinvBins);
    std::cout << par[3] << " " << FillValue << std::endl;
    inHist->Fill((xPtBins[i] + xPtBins[i+1])/2., FillValue);
    inHist->SetBinError(i+1, FillValue*std::pow((std::pow(FitResult->ParError(5)/par[5], 2) + std::pow(FitResult->ParError(3)/par[3], 2)), 0.5));

    // Addition of fitted signal and background curves to our canvas
    TF1 *backFcn = new TF1("backFcn",background,-0.03,0.03,3);
    TF1 *signalFcn = new TF1("signalFcn",GAUSS,-0.03,0.03,3);
    
    signalFcn->SetNpx(1e3);
    fitFcn->GetParameters(par);

    backFcn->SetParameters(par);
    backFcn->SetLineStyle(2);
    backFcn->SetLineColor(8);
    backFcn->Draw("same");

    signalFcn->SetLineColor(kMagenta+1);
    signalFcn->SetParameters(&par[3]);
    signalFcn->Draw("same");

    // Settings of Legend
    TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->SetLineColorAlpha(0.,0.);
    legend->SetFillColorAlpha(0.,0.);
    legend->AddEntry(hProfileInvMassX,"pp data","lpe");
    legend->AddEntry(backFcn,"pol2 background fit","l");
    legend->AddEntry(signalFcn,"gauss signal fit","l");
    legend->AddEntry(fitFcn,"global Fit","l");
    legend->Draw("same");

    // Settings of latex label of pt range
    TLatex ptRange;
    ptRange.SetTextSize(0.03);
    ptRange.SetTextFont(42);
    ptRange.SetNDC();
    ptRange.DrawLatex(0.194491, 0.74, Form("%.2f < p_{T} < %.2f", xPtBins[i], xPtBins[i+1]));

    // Settings of stat box
    gPad->Update();
    TPaveStats *st = (TPaveStats*)hProfileInvMassX->FindObject("stats");
    st->SetX1NDC(0.191152);
    st->SetX2NDC(0.333055);
    st->SetY1NDC(0.784897);
    st->SetY2NDC(0.847826);
    st->SetBorderSize(0);

    // Write canvas of fitted inv mass distr
    c1->Update();
    c1->Write();
  }
  gStyle->SetOptStat(0);
  // Go back to default dir in output file
  dir->cd("/");
}

void WriteToFile(TFile* outFile, Bool_t isMC = kFALSE){
  // Write to file
  outFile->cd();
  TDirectory* dirOut = outFile->mkdir("PtHists");
  dirOut->cd();
  hOmegaMB->Write();
  hOmegaHM->Write();
  hOmegaVHM->Write();

  if(isMC){
    TDirectory* dirOutGen = outFile->mkdir("PtHistsGen");
    dirOutGen->cd();
    hGenOmegaMB->Write();
    hGenOmegaHM->Write();
    hGenOmegaVHM->Write();
  }
  outFile->Close();
}

void make_results(const Char_t* fileNameData, const Char_t* fileNameEff, const Char_t* outputFileName, Bool_t isMC = kFALSE){
  gStyle->SetOptStat(0);
  // Get Efficiency
  fileEff = FindFileFresh(fileNameEff);
  fileEff->cd();
  hEffOmegaMB = (TH1D*)fileEff->Get("hEffOmegaMB");
  hEffOmegaHM = (TH1D*)fileEff->Get("hEffOmegaHM");
  hEffOmegaVHM = (TH1D*)fileEff->Get("hEffOmegaVHM");

  // Get Data
  TFile* fileData = FindFileFresh(fileNameData);
  fileData->cd();
  TList *ListOfHists = (TList*)fileData->Get("ListOfEInvMassHists");
  TH3D* hInvMassOmega = (TH3D*)ListOfHists->FindObject("hOmegaInvMassVsPt_Omega");
  TH3D* hInvMassOmegaBar = (TH3D*)ListOfHists->FindObject("hOmegaInvMassVsPt_OmegaBar");
  // Make summed hsitogram of Omega and OmegaBar
  TH3D* hInvMassSum = (TH3D*)hInvMassOmega->Clone();
  TH3D* histToAdd  = (TH3D*)hInvMassOmegaBar->Clone();
  hInvMassSum->Add(histToAdd);

  if(isMC){
    ListOfPtHists = (TList*)fileData->Get("ListOfPtHists");
    hGenOmegaMB = (TH1D*)ListOfPtHists->FindObject("hGenOmegaMB");
    hGenOmegaHM = (TH1D*)ListOfPtHists->FindObject("hGenOmegaHM");
    hGenOmegaVHM = (TH1D*)ListOfPtHists->FindObject("hGenOmegaVHM");
  }

  // Setup hists
  hOmegaMB = new TH1D("hOmegaMB", "; p_{T} [GeV/c]", nPtBinsMB, xBinsMB);
  hOmegaMB->Sumw2();

  hOmegaHM = new TH1D("hOmegaHM", ";  p_{T} [GeV/c]", nPtBinsHM, xBinsHM);
  hOmegaHM->Sumw2();

  hOmegaVHM = new TH1D("hOmegaVHM", ";  p_{T} [GeV/c]", nPtBinsHM, xBinsHM);
  hOmegaVHM->Sumw2();

  // Create output file
  outFile = new TFile(outputFileName, "RECREATE");
  // Extract the signal from 3d inv mass hist
  SignalExtractionPt(xBinsMB, nPtBinsMB, hInvMassSum, 1, 11, hOmegaMB, outFile); // 1 - 11 means 0 - 100 %
  SignalExtractionPt(xBinsHM, nPtBinsHM, hInvMassSum, 1, 4, hOmegaHM, outFile); // 1 - 4 means 0 - 10 %
  SignalExtractionPt(xBinsHM, nPtBinsHM, hInvMassSum, 1, 2, hOmegaVHM, outFile); // 1 - 2 means 0 - 1 %

  TF1* fRap = new TF1("fRap", rap_correction, 0.0, 50.0, 2);
  fRap->SetParameters(0.8, massOmega);

  const Int_t nSpectra = 3;
  TH1D* hEff[nSpectra] = { hEffOmegaMB, hEffOmegaHM, hEffOmegaVHM };
  TH1D* hRaw[nSpectra] = { hOmegaMB,    hOmegaHM,    hOmegaVHM };
  TH1D* hGen[nSpectra] = { hGenOmegaMB,    hGenOmegaHM,    hGenOmegaVHM};

  // Normalize results
  TH1D* hNorm = (TH1D*)fileData->Get("hNorm");
  {
    Double_t nMB = hNorm->GetBinContent(3);
    const Double_t vtxScale = (hNorm->GetBinContent(1)+hNorm->GetBinContent(2)+hNorm->GetBinContent(3)) /
      (hNorm->GetBinContent(2)+hNorm->GetBinContent(3));
    nMB *= vtxScale;
    hOmegaMB->Scale(1.0/nMB);
    if(isMC)
      hGenOmegaMB->Scale(1.0/nMB);
  }
  {
    Double_t nHM = hNorm->GetBinContent(4);
    hOmegaHM->Scale(1.0/nHM);
    if(isMC)
      hGenOmegaHM->Scale(1.0/nHM);
  }
  {
    Double_t nVHM = hNorm->GetBinContent(5);
    hOmegaVHM->Scale(1.0/nVHM);
    if(isMC)
      hGenOmegaVHM->Scale(1.0/nVHM);
  }

  // Normalize and correct histograms
  for(Int_t i = 0; i < nSpectra; i++) {
    NormalizeHistogram(hRaw[i]);
    hRaw[i]->Divide(hEff[i]);
    // Apply rapidity correction
    hRaw[i]->Divide(fRap);
    hRaw[i]->GetYaxis()->SetTitle("(#Omega+#bar{#Omega}): dN/dp_{T}");
    if(isMC){
      NormalizeHistogram(hGen[i]);
      hGen[i]->Divide(fRap);
      hGen[i]->GetYaxis()->SetTitle("(#Omega+#bar{#Omega}): dN/dp_{T}");
    }
  }
  WriteToFile(outFile, isMC);
}