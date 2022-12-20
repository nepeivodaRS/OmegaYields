// #include "helpFunc.C"
#include "ptDistr.h"

/*
  .L ptDistr.C+

  make_results("./outputAnal/mc_anal.root", "./outputPtHists")
  
 */

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

void SignalExtractionPt(const Double_t *xPtBins, const Int_t nPtBins, TH3D *invMassHist, Int_t leftCentr, Int_t rightCentr, TH1D* inHist, const Char_t* outputDirName){
  TFile f(Form("%s/PtFitHists.root", outputDirName),"RECREATE");
  invMassHist->GetZaxis()->SetRange(leftCentr, rightCentr);
  TH2D* hProfileInvMassZ = static_cast<TH2D*>(invMassHist->Project3D("xy"));
  for(Int_t i = 0; i < nPtBins - 1; ++i) {
    gROOT->SetBatch(kFALSE);
    TH1D* hProfileInvMassX = hProfileInvMassZ->ProjectionX("_px", i, i+1);
    hProfileInvMassX->SetStats(0);
    hProfileInvMassX->SetTitle(Form("M_{inv} vs Pt Fit in %d bin",  i));
    TF1 *fitFcn = new TF1("fitFcn",fitFunctionG,-0.03, 0.03,6);
    fitFcn->SetParameters(1,1,1,20,0,0.001);
    fitFcn->SetNpx(1e4);
    gROOT->SetBatch(kTRUE);
    TCanvas *c1 = new TCanvas(Form("PtFit_%d", i),"PtHist",10,10,1200,900);
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
    legend->SetTextFont(60);
    legend->SetTextSize(0.03);
    legend->SetLineColorAlpha(0.,0.);
    legend->SetFillColorAlpha(0.,0.);
    legend->AddEntry(hProfileInvMassX,"pp data","lpe");
    legend->AddEntry(backFcn,"pol2 background fit","l");
    legend->AddEntry(signalFcn,"gauss signal fit","l");
    legend->AddEntry(fitFcn,"global Fit","l");
    legend->Draw();
    c1->Write();
  }
  f.Close();
}

void WriteToFile(TFile* outFile){
  // Write to file
  outFile->cd();
  outFile->WriteObject(PtHistList, "ListOfPtHists");
  outFile->Close();
}

void make_results(const Char_t* fileNameData, const Char_t* outputDir){
  PtHistList = new TList();
  outFile = new TFile(Form("%s/PtHist.root", outputDir), "RECREATE");
  TFile* fileData = FindFileFresh(fileNameData);

  TList *ListOfHists = (TList*)fileData->Get("ListOfEInvMassHists");
  TH3D* hInvMassOmega = (TH3D*)ListOfHists->FindObject("hOmegaInvMassVsPt_Omega");
  TH3D* hInvMassOmegaBar = (TH3D*)ListOfHists->FindObject("hOmegaInvMassVsPt_OmegaBar");
  hInvMassOmega->Add(hInvMassOmegaBar);
  TF1* fRap = new TF1("fRap", rap_correction, 0.0, 50.0, 2);
  fRap->SetParameters(0.8, massOmega);
  TH1D* hOmegaMB = new TH1D("hOmegaMB", "; p_{T} [GeV/c]", nPtBinsMB, xBinsMB);
  hOmegaMB->Sumw2();

  SignalExtractionPt(xBinsMB, nPtBinsMB, hInvMassOmega, 1, 11, hOmegaMB, outputDir); // 1 - 11
  // hOmegaMB->Scale(1./nMB);
  // NormalizeHistogram(hOmegaMB);
  // Eff
  // hOmegaMB->Divide(fRap);
  PtHistList->Add(hOmegaMB);
  WriteToFile(outFile);
}