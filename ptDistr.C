// #include "helpFunc.C"
#include "ptDistr.h"

/*
  .L ptDistr.C

  make_results("./outputAnal/mc_anal_27_MCclosureFixed.root", "./outputEff/mc_Eff_26.root", "./outputPtHists/PtHist_28.root", 1)

  make_results("./outputAnal/mc_anal_2feb_injected.root", "./outputEff/mc_Eff_2feb_injected.root", "./outputPtHists/PtHist_2feb_injected.root", 1)
  
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

Double_t GAUSSredifined(Double_t *x, Double_t *par){
return par[0]*(0.06/nMinvBins)/(par[2]*TMath::Power(TMath::TwoPi(), 0.5))*exp(-0.5*(((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2])));
}

Double_t fitFunctionG(Double_t *x, Double_t *par) {
  return background(x,par) + GAUSS(x,&par[3]);
}

Double_t fitFunctionRedifined(Double_t *x, Double_t *par) {
  return background(x,par) + GAUSSredifined(x,&par[3]);
}

Double_t fitFunctionCB(Double_t *x, Double_t *par) {
  return background(x,par) + DoubleSidedCB(x,&par[3]);
}

void CreateRatioPlot(TH1D *h1In,TH1D *h2In, TFile* outFile){
  gROOT->SetBatch(kTRUE);
  TH1D* h1 = (TH1D*)h1In->Clone();
  TH1D* h2 = (TH1D*)h2In->Clone();
  Double_t MaxHeightHist = (h1->GetMaximum() > h2->GetMaximum()) ? h1->GetMaximum() : h2->GetMaximum();

  TDirectory* dir = outFile->mkdir(Form("closureMC_%s", h1->GetName()));
  dir->cd();
  h1->Write();
  h2->Write();
  TCanvas *canvas= new TCanvas(h1->GetName(), "canvas", 10,10,1200,900);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  //pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  h1->SetStats(0);          // No statistics on upper plot
  h1->Draw();               // Draw h1
  h2->Draw("same");         // Draw h2 on top of h1

  TAxis *axis = h1->GetYaxis();
  //axis->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
  axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axis->SetLabelSize(19);
  axis->SetRange(0, MaxHeightHist);

  TLegend*  lLegend     = new TLegend(0.60,0.658,0.782,0.759);
  lLegend->SetTextFont(42);
  lLegend->SetLineColorAlpha(0.,0.);
  lLegend->SetFillColorAlpha(0.,0.);
  lLegend->SetBorderSize(0.);
  lLegend->SetTextSize(0.04);
  lLegend->AddEntry  (h1,h1->GetName(),"lpf");
  lLegend->AddEntry  (h2,h2->GetName(),"lpf");
  lLegend->Draw();

  // SqrtSnn ALICE label
  TLatex sqrtSnn;
  sqrtSnn.SetTextSize(0.04);
  sqrtSnn.SetTextFont(42);
  sqrtSnn.SetNDC();
  sqrtSnn.DrawLatex(0.60, 0.79, "ALICE pp #sqrt{#it{s}} = 13 TeV");

  // MC closure label
  TLatex mcClosure;
  mcClosure.SetTextSize(0.04);
  mcClosure.SetTextFont(42);
  mcClosure.SetNDC();
  mcClosure.DrawLatex(0.60, 0.59, "MC closure");

  // lower plot will be in pad
  canvas->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridy(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  // Define the ratio plot
  TH1F *h3 = (TH1F*)h1->Clone("h3");
  h3->SetLineColor(kBlack);
  h3->SetMinimum(0.75);  // Define Y ..
  h3->SetMaximum(1.25); // .. range
  h3->SetStats(0);      // No statistics on lower plot
  h3->Divide(h2);
  //h3->SetMarkerStyle(21);
  h3->Draw("ep");       // Draw the ratio plot

  // h1 settings
  h1->SetLineColor(kBlue+1);
  h1->SetLineWidth(2);

  // Y axis h1 plot settings
  h1->GetYaxis()->SetTitleSize(22);
  h1->GetYaxis()->SetTitleFont(43);
  h1->GetYaxis()->SetTitleOffset(1.55);

  // h2 settings
  h2->SetLineColor(kRed);
  h2->SetLineWidth(2);

  // Ratio plot (h3) settings
  h3->SetTitle(""); // Remove the ratio title

  // Y axis ratio plot settings
  h3->GetYaxis()->SetTitle("Rec/Gen");
  h3->GetYaxis()->SetNdivisions(5);
  h3->GetYaxis()->SetTitleSize(22);
  h3->GetYaxis()->SetTitleFont(43);
  h3->GetYaxis()->SetTitleOffset(1.55);
  h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h3->GetYaxis()->SetLabelSize(19);
  h3->GetYaxis()->CenterTitle(true);
  //TAxis *axisH3 = h3->GetYaxis();
  //axisH3->ChangeLabel(1, -1, -1, -1, -1, -1, " ");

  // X axis ratio plot settings
  h3->GetXaxis()->SetTitleSize(22);
  h3->GetXaxis()->SetTitleFont(43);
  h3->GetXaxis()->SetTitleOffset(3.5);
  h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h3->GetXaxis()->SetLabelSize(19);

  canvas->Update();
  canvas->Write();
  gROOT->SetBatch(kFALSE);
  dir->cd("/");
}

void SignalExtractionPt(const Double_t *xPtBins, const Int_t nPtBins, TH3D *inHist3D, Int_t leftCentr, Int_t rightCentr, TH1D* inHist, TFile* outFile){
  // Set stat box to show only mean and number of entries
  gStyle->SetOptStat("me");
  gStyle->SetOptFit(1);
  // Setup for Gaussian mean graph
  double bins[7] = {1, 2, 3, 4, 5, 6, 7};
  double mean[7] = {0};
  double sigma[7] = {0};
  double errBins[7] = {0};
  double errMean[7] = {0};
  double errSigma[7] = {0};
  // Setup for Signal graph
  double signal[7] = {0};
  double errSignal[7] = {0};
  // Create dir to store all the fitted hists in centrality region from 'leftCentr' to 'rightCentr' bin
  TDirectory* dir = outFile->mkdir(Form("PtFitHists_%s_from_%d_to_%d_centr", inHist3D->GetName(), leftCentr, rightCentr));
  dir->cd();
  // Clone hist not to change the original one
  TH3D* invMassHist = (TH3D*)inHist3D->Clone();
  invMassHist->GetZaxis()->SetRange(leftCentr, rightCentr);
  invMassHist->Write();
  TH2D* hProfileInvMassZ = static_cast<TH2D*>(invMassHist->Project3D("xy")); // doesn't work w/o static cast; projects in range that was set above
  hProfileInvMassZ->Write();
  // Loop over all PtBins and fit inv mass spectra
  for(Int_t i = 0; i < nPtBins; ++i) {
    gROOT->SetBatch(kTRUE);
    TH1D* hProfileInvMassX = hProfileInvMassZ->ProjectionX("_px", hProfileInvMassZ->GetYaxis()->FindBin(xPtBins[i] + 0.00001), hProfileInvMassZ->GetYaxis()->FindBin(xPtBins[i+1]-0.00001));
    hProfileInvMassX->SetTitle(Form("M_{inv} in pt bins fit, %d bin",  i+1));
    hProfileInvMassX->SetXTitle("#it{M}_{inv} - #it{M}_{#Omega^{-} (#bar{#Omega}^{+})} [GeV/#it{c}^{2}]");
    TF1 *fitFcn = new TF1("fitFcn",fitFunctionG,-0.03, 0.03,6);
    Double_t sigPickApprox = 20; // 0 guess for the first fit
    fitFcn->SetParameters(1,1,1,sigPickApprox,0,0.001);
    // fitFcn->SetNpx(1e5);
    // fitFcn->SetLineColor(kRed+1);
    TH1D* hProfileInvMassXFirstFit = (TH1D*)hProfileInvMassX->Clone();
    hProfileInvMassXFirstFit->Fit("fitFcn","eRL");
    Double_t par[6];
    fitFcn->GetParameters(par);
    par[5] = TMath::Abs(par[5]); // sometimes sigma is negative
    par[3] = par[3]*par[5]*TMath::Power(TMath::TwoPi(), 0.5)/(0.06/nMinvBins); // from initial guess of A to S
    // Make a fit once again with initial guess based on the previous fit
    TF1 *fitFcnRedefined = new TF1("fitFcn2",fitFunctionRedifined,-0.03, 0.03,6);
    fitFcnRedefined->SetParameters(par);
    fitFcnRedefined->SetNpx(1e5);
    fitFcnRedefined->SetLineColor(kRed+1);
    // Create canvas for fitted inv mass
    TCanvas *c1 = new TCanvas(Form("InvMassFit_bin_%d_from_%d_to_%d_centr", i+1, leftCentr, rightCentr),"PtHist",10,10,1200,900);
    c1->cd();
    TFitResultPtr FitResult = hProfileInvMassX->Fit("fitFcn2","eRSL");
    fitFcnRedefined->GetParameters(par);
    // Fill histogram with fitted signal: integral divided by the bin width. And eval the error of each bin
    Double_t FillValue = par[3];
    std::cout << "bin number: " << i+1 << " signal value: " << FillValue << std::endl;
    inHist->Fill((xPtBins[i] + xPtBins[i+1])/2., FillValue);
    inHist->SetBinError(i+1, FitResult->ParError(3));
    // Points for Gaussian mean evolution graph
    mean[i] = par[4];
    errMean[i] = FitResult->ParError(4);
    // Points for Gaussian sigma evolution graph
    sigma[i] = par[5];
    errSigma[i] = FitResult->ParError(5);
    // Points for signal graph
    signal[i] = par[3];
    errSignal[i] = FitResult->ParError(3);
    // Addition of fitted signal and background curves to our canvas
    TF1 *backFcn = new TF1("backFcn",background,-0.03,0.03,3);
    TF1 *signalFcn = new TF1("signalFcn",GAUSSredifined,-0.03,0.03,3);
    
    signalFcn->SetNpx(1e4);
    backFcn->SetParameters(par);
    backFcn->SetLineStyle(2);
    backFcn->SetLineColor(kCyan+2); //{kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
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
    legend->SetBorderSize(0.);
    legend->AddEntry(hProfileInvMassX,"Data (stat uncert.)","lpe");
    legend->AddEntry(fitFcnRedefined,"fit (signal + bkg.)","l");
    legend->AddEntry(backFcn,"estimated bkg. (pol2)","l");
    legend->AddEntry(signalFcn,"estimated signal (Gauss)","l");
    legend->Draw("same");

    // Settings of latex label of pt range
    TLatex ptRange;
    ptRange.SetTextSize(0.03);
    ptRange.SetTextFont(42);
    ptRange.SetNDC();
    ptRange.DrawLatex(0.155, 0.814, Form("%.2f GeV/#it{c} < #it{p}_{T} < %.2f GeV/#it{c}", xPtBins[i], xPtBins[i+1]));

    // Settings of latex label of SqrtSnn
    TLatex sqrtSnn;
    sqrtSnn.SetTextSize(0.03);
    sqrtSnn.SetTextFont(42);
    sqrtSnn.SetNDC();
    sqrtSnn.DrawLatex(0.155, 0.850, "ALICE pp #sqrt{#it{s}} = 13 TeV");

    // Settings of latex label of OmegaLabel
    TLatex omegaLabel;
    omegaLabel.SetTextSize(0.0411899);
    omegaLabel.SetTextFont(42);
    omegaLabel.SetNDC();
    omegaLabel.DrawLatex(0.1569, 0.574, "#Omega^{-}+#bar{#Omega}^{+}");

    // Settings of stat box
    gPad->Update(); // update to find 'stats' box
    TPaveStats *st = (TPaveStats*)hProfileInvMassX->FindObject("stats");
    //st->SetTextSize(0.03);
    st->SetX1NDC(0.146);
    st->SetX2NDC(0.363);
    st->SetY1NDC(0.617);
    st->SetY2NDC(0.801);
    st->SetBorderSize(0);

    // Write canvas of fitted inv mass distr
    gROOT->SetBatch(kFALSE);
    gPad->Update();
    c1->Update();
    c1->Write();
  }
  // Plot the Gaussian mean fit value evolution
  gROOT->SetBatch(kTRUE);
  TGraphErrors *MeanGaussFit = new TGraphErrors(7, bins, mean, errBins, errMean);
  MeanGaussFit->SetTitle(Form("#mu of the Gaussain fit for %s from %d to %d mult", inHist3D->GetName(), leftCentr, rightCentr));
  MeanGaussFit->GetXaxis()->SetTitle("Bin number");
  MeanGaussFit->GetYaxis()->SetTitle("#mu");
  MeanGaussFit->SetMarkerStyle(20);
  TCanvas *c2 = new TCanvas(Form("#mu of the Gaussain fit for %s from %d to %d mult", inHist3D->GetName(), leftCentr, rightCentr),"PtHist",10,10,1200,900);
  c2->cd();
  MeanGaussFit->Draw("AP");
  c2->Update();
  gROOT->SetBatch(kFALSE);
  c2->Write();
  // Plot the Gaussian sigma fit value evolution
  gROOT->SetBatch(kTRUE);
  TGraphErrors *SigmaGaussFit = new TGraphErrors(7, bins, sigma, errBins, errSigma);
  SigmaGaussFit->SetTitle(Form("#sigma of the Gaussain fit for %s from %d to %d mult", inHist3D->GetName(), leftCentr, rightCentr));
  SigmaGaussFit->GetXaxis()->SetTitle("Bin number");
  SigmaGaussFit->GetYaxis()->SetTitle("#sigma");
  SigmaGaussFit->SetMarkerStyle(20);
  TCanvas *c4 = new TCanvas(Form("#sigma of the Gaussain fit for %s from %d to %d mult", inHist3D->GetName(), leftCentr, rightCentr),"PtHist",10,10,1200,900);
  c4->cd();
  SigmaGaussFit->Draw("AP");
  c4->Update();
  gROOT->SetBatch(kFALSE);
  c4->Write();
  // Plot the signal evolution
  gROOT->SetBatch(kTRUE);
  TGraphErrors *NumberOfCascades = new TGraphErrors(7, bins, signal, errBins, errSignal);
  NumberOfCascades->SetTitle(Form("Number_of_cascades_for_%s_from_%d_to_%d_mult", inHist3D->GetName(), leftCentr, rightCentr));
  NumberOfCascades->GetXaxis()->SetTitle("Bin number");
  NumberOfCascades->GetYaxis()->SetTitle("Number of cascades");
  NumberOfCascades->SetMarkerStyle(20);
  TCanvas *c3 = new TCanvas(Form("Number_of_cascades_for_%s_from_%d_to_%d_mult", inHist3D->GetName(), leftCentr, rightCentr),"PtHist",10,10,1200,900);
  c3->cd();
  NumberOfCascades->Draw("AP");
  c3->Update();
  gROOT->SetBatch(kFALSE);
  c3->Write();
  // Go back to default dir in output file
  dir->cd("/");
}

void SignalMC(const Double_t *xPtBins, const Int_t nPtBins, TH3D *inHist3D, Int_t leftCentr, Int_t rightCentr, TH1D* inHist, TFile* outFile){
  // Set stat box to show only mean and number of entries
  gStyle->SetOptStat("me");
  // Setup for Signal graph
  double bins[7] = {1, 2, 3, 4, 5, 6, 7};
  double errBins[7] = {0};
  double signal[7] = {0};
  double errSignal[7] = {0};
  // Create dir to store all the fitted hists in centrality region from 'leftCentr' to 'rightCentr' bin
  TDirectory* dir = outFile->mkdir(Form("PtFitHists_%s_from_%d_to_%d_centr", inHist3D->GetName(), leftCentr, rightCentr));
  dir->cd();
  // Clone hist not to change the original one
  TH3D* invMassHist = (TH3D*)inHist3D->Clone();
  invMassHist->GetZaxis()->SetRange(leftCentr, rightCentr);
  invMassHist->Write();
  TH2D* hProfileInvMassZ = static_cast<TH2D*>(invMassHist->Project3D("xy")); // doesn't work w/o static cast; projects in range that was set above
  hProfileInvMassZ->Write();
  // Loop over all PtBins and fit inv mass spectra
  for(Int_t i = 0; i < nPtBins; ++i) {
    gROOT->SetBatch(kTRUE);
    TH1D* hProfileInvMassX = hProfileInvMassZ->ProjectionX("_px", hProfileInvMassZ->GetYaxis()->FindBin(xPtBins[i] + 0.00001), hProfileInvMassZ->GetYaxis()->FindBin(xPtBins[i+1]-0.00001));
    hProfileInvMassX->SetTitle(Form("M_{inv} in pt bins fit, %d bin",  i+1));
    hProfileInvMassX->SetXTitle("#it{M}_{inv} - #it{M}_{#Omega^{-} (#bar{#Omega}^{+})} [GeV/#it{c}^{2}]");

    TCanvas *c1 = new TCanvas(Form("InvMass_bin_%d_from_%d_to_%d", i+1, leftCentr, rightCentr),"PtHist",10,10,1200,900);
    c1->cd();
    // Points for signal graph
    signal[i] = hProfileInvMassX->GetEntries();
    errSignal[i] = 0;

    hProfileInvMassX->Draw();

    // Settings of Legend
    TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->SetLineColorAlpha(0.,0.);
    legend->SetFillColorAlpha(0.,0.);
    legend->SetBorderSize(0.);
    legend->AddEntry(hProfileInvMassX,"Real signal","lpe");
    legend->Draw("same");

    // Settings of latex label of pt range
    TLatex ptRange;
    ptRange.SetTextSize(0.03);
    ptRange.SetTextFont(42);
    ptRange.SetNDC();
    ptRange.DrawLatex(0.155, 0.814, Form("%.2f GeV/#it{c} < #it{p}_{T} < %.2f GeV/#it{c}", xPtBins[i], xPtBins[i+1]));

    // Settings of latex label of SqrtSnn
    TLatex sqrtSnn;
    sqrtSnn.SetTextSize(0.03);
    sqrtSnn.SetTextFont(42);
    sqrtSnn.SetNDC();
    sqrtSnn.DrawLatex(0.155, 0.850, "ALICE pp #sqrt{#it{s}} = 13 TeV");

    // Settings of latex label of OmegaLabel
    TLatex omegaLabel;
    omegaLabel.SetTextSize(0.0411899);
    omegaLabel.SetTextFont(42);
    omegaLabel.SetNDC();
    omegaLabel.DrawLatex(0.1569, 0.574, "#Omega^{-}+#bar{#Omega}^{+}");

    // Settings of stat box
    gPad->Update(); // update to find 'stats' box
    TPaveStats *st = (TPaveStats*)hProfileInvMassX->FindObject("stats");
    //st->SetTextSize(0.03);
    st->SetX1NDC(0.146);
    st->SetX2NDC(0.363);
    st->SetY1NDC(0.617);
    st->SetY2NDC(0.801);
    st->SetBorderSize(0);

    // Write canvas of MC inv mass distr
    gROOT->SetBatch(kFALSE);
    gPad->Update();
    c1->Update();
    c1->Write();
  }
  // Plot the signal evolution
  gROOT->SetBatch(kTRUE);
  TGraphErrors *NumberOfCascades = new TGraphErrors(7, bins, signal, errBins, errSignal);
  NumberOfCascades->SetTitle(Form("Number of cascades for %s from %d to %d mult", inHist3D->GetName(), leftCentr, rightCentr));
  NumberOfCascades->GetXaxis()->SetTitle("Bin number");
  NumberOfCascades->GetYaxis()->SetTitle("Number of cascades");
  NumberOfCascades->SetMarkerStyle(20);
  TCanvas *c3 = new TCanvas(Form("Number_of_cascades_for_%s_from_%d_to_%d_mult", inHist3D->GetName(), leftCentr, rightCentr),"PtHist",10,10,1200,900);
  c3->cd();
  NumberOfCascades->Draw("AP");
  c3->Update();
  gROOT->SetBatch(kFALSE);
  c3->Write();
  // Go back to default dir in output file
  dir->cd("/");
}

void WriteToFile(TFile* outFile, Bool_t isMC = kFALSE){
  // Write to file
  outFile->cd();
  // Fill reconstructed pt hists
  TDirectory* dirOut = outFile->mkdir("PtHists");
  dirOut->cd();
  hOmegaMB->Write();
  hOmegaHM->Write();
  hOmegaVHM->Write();

  // Fill MC generated pt hists
  if(isMC){
    TDirectory* dirOutGen = outFile->mkdir("PtHistsGen");
    dirOutGen->cd();
    hGenOmegaMB->Write();
    hGenOmegaHM->Write();
    hGenOmegaVHM->Write();
  }

  // Fill efficiency pt hists
  TDirectory* dirOutEff = outFile->mkdir("PtHistsEff");
  dirOutEff->cd();
  hEffOmegaMB->Write();
  hEffOmegaHM->Write();
  hEffOmegaVHM->Write();

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

    TH3D* hInvMassOmegaMC = (TH3D*)ListOfHists->FindObject("hOmegaInvMassVsPtTrue_Omega");
    TH3D* hInvMassOmegaBarMC = (TH3D*)ListOfHists->FindObject("hOmegaInvMassVsPtTrue_OmegaBar");
    hInvMassSumMC = (TH3D*)hInvMassOmegaMC->Clone();
    TH3D* histToAddMC  = (TH3D*)hInvMassOmegaBarMC->Clone();
    hInvMassSumMC->Add(histToAddMC);
  }

  // Setup hists
  hOmegaMB = new TH1D("hOmegaMB", "; #it{p}_{T} (GeV/c)", nPtBinsMB, xBinsMB);
  hOmegaMB->Sumw2();

  hOmegaMBMC = new TH1D("hOmegaMBMC", "; #it{p}_{T} (GeV/c)", nPtBinsMB, xBinsMB);
  hOmegaMBMC->Sumw2();

  hOmegaHM = new TH1D("hOmegaHM", ";  #it{p}_{T} (GeV/c)", nPtBinsHM, xBinsHM);
  hOmegaHM->Sumw2();

  hOmegaVHM = new TH1D("hOmegaVHM", ";  #it{p}_{T} (GeV/c)", nPtBinsHM, xBinsHM);
  hOmegaVHM->Sumw2();

  // Create output file
  outFile = new TFile(outputFileName, "RECREATE");
  // Extract the signal from 3d inv mass hist
  SignalExtractionPt(xBinsMB, nPtBinsMB, hInvMassSum, 1, 11, hOmegaMB, outFile); // 1 - 11 means 0 - 100 %
  SignalExtractionPt(xBinsHM, nPtBinsHM, hInvMassSum, 1, 3, hOmegaHM, outFile); // 1 - 3 means 0 - 10 %
  SignalExtractionPt(xBinsHM, nPtBinsHM, hInvMassSum, 1, 1, hOmegaVHM, outFile); // 1 - 1 means 0 - 1 %

  if(isMC){
    SignalMC(xBinsMB, nPtBinsMB, hInvMassSumMC, 1, 11, hOmegaMBMC, outFile); // 1 - 11 means 0 - 100 %
  }

  // Setup rapidity correction TF1
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
    hRaw[i]->GetYaxis()->SetTitle("(#Omega^{-}+#bar{#Omega}^{+}):  d^{2}#it{N}/d#it{p}_{T} ((GeV/#it{c})^{-1})");
    if(isMC){
      NormalizeHistogram(hGen[i]);
      hGen[i]->Divide(fRap);
      hGen[i]->GetYaxis()->SetTitle("(#Omega^{-}+#bar{#Omega}^{+}): d^{2}#it{N}/d#it{p}_{T} ((GeV/#it{c})^{-1})");
      CreateRatioPlot(hRaw[i], hGen[i], outFile);
    }
  }
  WriteToFile(outFile, isMC);
}

