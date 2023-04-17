#include "closureAndDiffRecComp.h"

/*
  .L closureAndDiffRecComp.C

  closure_inv_mass("./outputPtHists/PtHist_6april_data.root", "/Users/rnepeiv/workLund/PhD_work/OmegaYields/published_data/correctedspectrumfittedgfcorrphyseff_spectra23.root", "./outputClosure/SignalClosure_6april_data.root", 0)
  closure_inv_mass("./outputPtHists/PtHist_6april_mc.root", "/Users/rnepeiv/workLund/PhD_work/OmegaYields/published_data/correctedspectrumfittedgfcorrphyseff_spectra23.root", "./outputClosure/SignalClosure_6april_mc.root", 1)
  
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

TCanvas* createRatioPlotRecOverGenMC(TGraphErrors* graph1, TGraphErrors* graph2) {
    gROOT->SetBatch(kTRUE);
    TGraphErrors* ratio = new TGraphErrors();
    TCanvas *canvas3 = new TCanvas("RecOverGenPlot","RecOverGenPlot",10,10,1200,900);
    int n = graph1->GetN();
    for (int i = 0; i < n; i++) {
        double x1, y1, x2, y2, xerr1, yerr1, xerr2, yerr2;
        graph1->GetPoint(i, x1, y1);
        graph2->GetPoint(i, x2, y2);
        xerr1 = graph1->GetErrorX(i);
        yerr1 = graph1->GetErrorY(i);
        xerr2 = graph2->GetErrorX(i);
        yerr2 = graph2->GetErrorY(i);
        if (y2 != 0) {
            double ratio_y = y1 / y2;
            //double ratio_yerr = ratio_y * sqrt(pow(yerr1 / y1, 2) + pow(yerr2 / y2, 2));
            ratio->SetPoint(i, x1, ratio_y);
            //ratio->SetPointError(i, xerr1, ratio_yerr);
        }
    }
    ratio->Draw("AP");
    ratio->SetMarkerStyle(20);
    ratio->GetXaxis()->SetTitle("Bin number");
    ratio->GetYaxis()->SetTitle("Rec/Gen");
    gROOT->SetBatch(kFALSE);
    return canvas3;
    // f->cd();
    // canvas3->Update();
    // canvas3->Write();
}

void createRatioPlotDiffRecMethods(TH1D *h1In,TH1D *h2In, TFile* File, const Char_t* CanvName, const Char_t* h1Name, const Char_t* h2Name){
  gROOT->SetBatch(kTRUE);
  TH1D* h1 = (TH1D*)h1In->Clone();
  TH1D* h2 = (TH1D*)h2In->Clone();
  Double_t MaxHeightHist = (h1->GetMaximum() > h2->GetMaximum()) ? h1->GetMaximum() : h2->GetMaximum();

  TDirectory* dir = File->mkdir(CanvName);
  dir->cd();
  h1->Write();
  h2->Write();
  TCanvas *canvas4= new TCanvas(CanvName, "canvas", 10,10,1200,900);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  //pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  h1->SetStats(0);          // No statistics on upper plot
  h2->SetStats(0);          // No statistics on upper plot

  TAxis *axis = h1->GetYaxis();
  //axis->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
  axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axis->SetLabelSize(19);
  axis->SetRange(0, MaxHeightHist);

  h2->Draw();         // Draw h2
  h1->Draw("same");               //  Draw h1 on top of h2

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
  mcClosure.DrawLatex(0.60, 0.59, "Rec. Comparison");

  // lower plot will be in pad
  canvas4->cd();          // Go back to the main canvas before defining pad2
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
  h3->GetYaxis()->SetTitle(Form("%s / %s", h1Name, h2Name));
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

  canvas4->Update();
  canvas4->Write();
  gROOT->SetBatch(kFALSE);
  dir->cd("/");
}

TCanvas* findCanvas(TDirectory* dir, const Char_t* canvasName){
  TKey* key;
  TCanvas *canvasReturn;
  TIter nextkey(dir->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
    //std::cout << key << std::endl;
    const char * classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname);
    if (!cl) continue;
    TObject *obj = key->ReadObj();
    if (obj->InheritsFrom("TCanvas")){
      TCanvas * can = (TCanvas*)obj;
      //std::cout << can->GetName() << std::endl;
      if(strcmp(can->GetName(), canvasName) == 0){
        canvasReturn = can;
      }
    }
  }
  return canvasReturn;
}

void WriteToFile(TFile* outFile, Bool_t isMC = kFALSE){
  outFile->cd();
  if(isMC){
    cRatioPlot->Write();
  }
}

void closure_inv_mass(const Char_t* fileNameData, const Char_t* publishedData, const Char_t* outputFileName, Bool_t isMC = kFALSE){
  // Create output file
  outFile = new TFile(outputFileName, "RECREATE");
  TFile* publishedFile = FindFileFresh(publishedData);
  // Get Data
  TFile* f = FindFileFresh(fileNameData);
  f->cd();
  //f->ls();
  if(isMC){
    gROOT->SetBatch(kTRUE);
    TCanvas *cRecNumberOfCasc = new TCanvas();
    TCanvas *cGenNumberOfcasc = new TCanvas();
    gROOT->SetBatch(kFALSE);
    // Search for dir with specific rec method
    //f->cd("DefPtFitHists_hOmegaInvMassVsPt_Omega_from_1_to_11_centr");
    //f->cd("FixedBGPtFitHists_hOmegaInvMassVsPt_Omega_from_1_to_11_centr");
    f->cd("CombinedPtFitHists_hOmegaInvMassVsPt_Omega_from_1_to_11_centr");

    // Get canvas of number of rec cascades
    cRecNumberOfCasc = findCanvas(gDirectory, "Number_of_cascades_Subtr_for_hOmegaInvMassVsPt_Omega_from_1_to_11_mult");
    //cRecNumberOfCasc->ls();
    //cRecNumberOfCasc->Draw();
    // Get canvas of number of gen cascades
    f->cd("PtFitHists_hOmegaInvMassVsPtTrueEffCorr_Omega_from_1_to_11_centr");
    cGenNumberOfcasc = findCanvas(gDirectory, "Number_of_cascades_for_hOmegaInvMassVsPtTrueEffCorr_Omega_from_1_to_11_mult");

    // f->cd("CombinedPtFitHists_hOmegaInvMassVsPt_Omega_from_1_to_11_centr");
    // cGenNumberOfcasc = findCanvas(gDirectory, "Number_of_cascades_for_hOmegaInvMassVsPt_Omega_from_1_to_11_mult");

    // Get graphs from canvases
    //std::cout << (TGraphErrors*)canvas1->GetListOfPrimitives() << std::endl;
    TGraphErrors* graphRec = (TGraphErrors*)cRecNumberOfCasc->GetListOfPrimitives()->FindObject("Graph");
    TGraphErrors* graphGen = (TGraphErrors*)cGenNumberOfcasc->GetListOfPrimitives()->FindObject("Graph");
    // Create ratio plot
    cRatioPlot = createRatioPlotRecOverGenMC(graphRec, graphGen);
  }
  //gDirectory->ls();
  // Setup hists
  //hist->Draw();
  hOmegaMBdef = (TH1D*)f->Get("PtHists/hOmegaMBdef");
  hOmegaMBSideBand = (TH1D*)f->Get("PtHists/hOmegaMBSideBand");
  hOmegaMBBGfix = (TH1D*)f->Get("PtHists/hOmegaMBBGfix");
  hOmegaMBCombined = (TH1D*)f->Get("PtHists/hOmegaMBCombined");

  hOmegaMBPublished = (TH1D*)publishedFile->Get("hptcorr_norm_5");
  //hOmegaMBdef->Draw();
  createRatioPlotDiffRecMethods(hOmegaMBdef, hOmegaMBSideBand, outFile, "DefToSideBand", "Def", "SideBand");
  createRatioPlotDiffRecMethods(hOmegaMBdef, hOmegaMBBGfix, outFile, "DefToFixBG", "Def", "FixBG");
  createRatioPlotDiffRecMethods(hOmegaMBSideBand, hOmegaMBBGfix, outFile, "SideBandToFixBG", "SideBand", "FixBG");

  createRatioPlotDiffRecMethods(hOmegaMBCombined, hOmegaMBPublished, outFile, "CombinedToPublished", "This", "Publ.");

  createRatioPlotDiffRecMethods(hOmegaMBBGfix, hOmegaMBPublished, outFile, "FixBGToPublished", "FixBG", "Publ.");
  createRatioPlotDiffRecMethods(hOmegaMBdef, hOmegaMBPublished, outFile, "DefToPublished", "Def", "Publ.");
  createRatioPlotDiffRecMethods(hOmegaMBSideBand, hOmegaMBPublished, outFile, "SideBandToPublished", "SideBand", "Publ.");

  WriteToFile(outFile, isMC);
  //TDirectory* dir = outFile->mkdir("RatioSignalClosure");
}
