#include "closure_inv_mass.h"

/*
  .L closure_inv_mass.C

  closure_inv_mass("./outputPtHists/PtHist_2march_mc.root", "./outputClosure/SignalClosure_2march_mc_def", 1)
  
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

void createRatioPlotMC(TGraphErrors* graph1, TGraphErrors* graph2) {
    TGraphErrors* ratio = new TGraphErrors();
    TCanvas *canvas3 = new TCanvas();
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
    ratio->GetYaxis()->SetTitle("Fit/Gen");
    // f->cd();
    // canvas3->Update();
    // canvas3->Write();
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

void closure_inv_mass(const Char_t* fileNameData, const Char_t* outputFileName, Bool_t isMC = kFALSE){
  // Get Data
  TFile* f = FindFileFresh(fileNameData);
  f->cd();
  outFile = new TFile(outputFileName, "RECREATE");
  //f->cd("SideBandPtFitHists_hOmegaInvMassVsPt_Omega_from_1_to_11_centr");
  f->cd("DefPtFitHists_hOmegaInvMassVsPt_Omega_from_1_to_11_centr");
  //f->ls();
  gROOT->SetBatch(kTRUE);
  TCanvas *canvas1 = new TCanvas();
  TCanvas *canvas2 = new TCanvas();
  gROOT->SetBatch(kFALSE);
  // if(graph1 && graph2)
  canvas1 = findCanvas(gDirectory, "Number_of_cascades_for_hOmegaInvMassVsPt_Omega_from_1_to_11_mult");
  // canvas1->Draw();
  f->cd("PtFitHists_hOmegaInvMassVsPtTrue_Omega_from_1_to_11_centr");
  canvas2 = findCanvas(gDirectory, "Number_of_cascades_for_hOmegaInvMassVsPtTrue_Omega_from_1_to_11_mult");
  canvas1->ls();
  std::cout << (TGraphErrors*)canvas1->GetListOfPrimitives() << std::endl;
  TGraphErrors* graph1 = (TGraphErrors*)canvas1->GetListOfPrimitives()->FindObject("Graph");
  //graph1->Draw();
  TGraphErrors* graph2 = (TGraphErrors*)canvas2->GetListOfPrimitives()->FindObject("Graph");
  createRatioPlotMC(graph1, graph2);
  // canvas2->Draw();
}
