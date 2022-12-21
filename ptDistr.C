// #include "helpFunc.C"
#include "ptDistr.h"

/*
  .L ptDistr.C+

  make_results("./outputAnal/mc_anal.root", "./outputPtHists")
  
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

void SignalExtractionPt(const Double_t *xPtBins, const Int_t nPtBins, TH3D *invMassHist, Int_t leftCentr, Int_t rightCentr, TH1D* inHist, TFile* outFile){
  TDirectory* dir = outFile->mkdir(Form("PtFitHists_%s", invMassHist->GetName()));
  dir->cd();
  invMassHist->GetZaxis()->SetRange(leftCentr, rightCentr);
  TH2D* hProfileInvMassZ = static_cast<TH2D*>(invMassHist->Project3D("xy"));
  for(Int_t i = 0; i < nPtBins - 1; ++i) {
    gROOT->SetBatch(kFALSE);
    TH1D* hProfileInvMassX = hProfileInvMassZ->ProjectionX("_px", i, i+1);
    hProfileInvMassX->SetStats(0);
    hProfileInvMassX->SetTitle(Form("M_{inv} in pt bins fit, %d bin",  i+1));
    TF1 *fitFcn = new TF1("fitFcn",fitFunctionG,-0.03, 0.03,6);
    fitFcn->SetParameters(1,1,1,20,0,0.001);
    fitFcn->SetNpx(1e4);
    gROOT->SetBatch(kTRUE);
    TCanvas *c1 = new TCanvas(Form("PtFit_%d", i+1),"PtHist",10,10,1200,900);
    c1->cd();
    TFitResultPtr FitResult = hProfileInvMassX->Fit("fitFcn","epS");
    Double_t par[6];
    fitFcn->GetParameters(par);
    // Fill histogram with fitted signal integral divided by the bin width
    Double_t FillValue = par[3]*TMath::Power(TMath::TwoPi(), 0.5)*par[5]/(0.06/nMinvBins);
    inHist->Fill((xPtBins[i] + xPtBins[i+1])/2., FillValue);
    inHist->SetBinError(i+1, FillValue*std::pow((std::pow(FitResult->ParError(5)/par[5], 2) + std::pow(FitResult->ParError(3)/par[3], 2)), 0.5));

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

    TLatex ptRange;
    ptRange.SetTextSize(0.03);
    ptRange.SetTextFont(42);
    ptRange.SetNDC();
    ptRange.DrawLatex(0.194491, 0.74, Form("%.2f < p_{T} < %.2f", xPtBins[i], xPtBins[i+1]));
    c1->Update();
    c1->Write();
  }
  dir->cd("/");
}

void WriteToFile(TFile* outFile){
  // Write to file
  outFile->cd();
  outFile->WriteObject(PtHistList, "ListOfPtHists");
  outFile->Close();
}

void make_results(const Char_t* fileNameData, const Char_t* outputDir, const Char_t* fileNameEff){
  // Get Efficiency
  fileEff = FindFileFresh(fileNameEff);
  hEffOmegaMB = (TH1D*)fileEff->Get("hEffOmegaMB");
  hEffOmegaHM = (TH1D*)fileEff->Get("hEffOmegaHM");
  hEffOmegaVHM = (TH1D*)fileEff->Get("hEffOmegaVHM");

  outFile = new TFile(Form("%s/PtHist.root", outputDir), "RECREATE");

  TFile* fileData = FindFileFresh(fileNameData);
  TList *ListOfHists = (TList*)fileData->Get("ListOfEInvMassHists");
  TH3D* hInvMassOmega = (TH3D*)ListOfHists->FindObject("hOmegaInvMassVsPt_Omega");
  TH3D* hInvMassOmegaBar = (TH3D*)ListOfHists->FindObject("hOmegaInvMassVsPt_OmegaBar");
  // Make summed hsitogram of Omega and OmegaBar
  hInvMassOmega->Add(hInvMassOmegaBar);

  hOmegaMB = new TH1D("hOmegaMB", "; p_{T} [GeV/c]", nPtBinsMB, xBinsMB);
  hOmegaMB->Sumw2();

  hOmegaHM = new TH1D("hGenOmegaHM", ";  p_{T} [GeV/c]", nPtBinsHM, xBinsHM);
  hOmegaHM->Sumw2();

  hOmegaVHM = new TH1D("hGenOmegaVHM", ";  p_{T} [GeV/c]", nPtBinsHM, xBinsHM);
  hOmegaVHM->Sumw2();

  SignalExtractionPt(xBinsMB, nPtBinsMB, hInvMassOmega, 1, 11, hOmegaMB, outFile); // 1 - 11 means 0 - 100 %
  SignalExtractionPt(xBinsHM, nPtBinsHM, hInvMassOmega, 1, 3, hOmegaHM, outFile); // 1 - 3 means 0 - 10 %
  SignalExtractionPt(xBinsHM, nPtBinsHM, hInvMassOmega, 1, 2, hOmegaVHM, outFile); // 1 - 2 means 0 - 5 %

  TF1* fRap = new TF1("fRap", rap_correction, 0.0, 50.0, 2);
  fRap->SetParameters(0.8, massOmega);

  // Normalize results
  TH1D* hNorm = (TH1D*)fileData->Get("hNorm");
  {
    Double_t nMB = hNorm->GetBinContent(3);
    const Double_t vtxScale = (hNorm->GetBinContent(1)+hNorm->GetBinContent(2)+hNorm->GetBinContent(3)) /
      (hNorm->GetBinContent(2)+hNorm->GetBinContent(3));
    nMB *= vtxScale;
    hXiMB->Scale(1.0/nMB);
  }
  {
    Double_t nHM = hNorm->GetBinContent(4);
    hXiHM->Scale(1.0/nHM);
  }
  {
    Double_t nVHM = hNorm->GetBinContent(5);
    hXiVHM->Scale(1.0/nVHM);
  }

  // Normalize and correct histograms
  for(Int_t i = 0; i < nSpectra; i++) {
    NormalizeHistogram(hRaw[i]);
    hRaw[i]->Divide(hEff[i]);
    // Apply rapidity correction
    hRaw[i]->Divide(fRap);
    hRaw[i]->GetYaxis()->SetTitle("(#Omega+#bar{#Omega}): dN/dp_{T}");
  }

  // hOmegaMB->Scale(1./nMB);
  // NormalizeHistogram(hOmegaMB);
  // Eff
  // hOmegaMB->Divide(fRap);
  PtHistList = new TList();
  PtHistList->Add(hOmegaMB);
  WriteToFile(outFile);
}