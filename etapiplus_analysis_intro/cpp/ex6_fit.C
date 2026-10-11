// =====================================================================
// Exercise 6: Fit the a0(980) peak.
//
// What this macro does:
//   1. Makes the eta pi+ mass histogram with the cuts from exercise 3,
//      weighted by accidweight (accidental-subtracted).
//   2. Fits the peak near 1 GeV with a Gaussian (the peak) plus a
//      straight line (the background underneath it).
//   3. Prints the fitted mass, width, fit quality and signal yield.
//
// How to run it (in a terminal):
//   cd tutorial/cpp
//   root -l -b -q ex6_fit.C
// To change the fit range (here 0.88 to 1.10 GeV), pass it in:
//   root -l -b -q 'ex6_fit.C("../../flat_pipneta__B4_v6.root", 0.85, 1.15)'
//
// What you should see: ROOT's fit printout, then a short summary with
// the mean, sigma, chi2/ndf and yield (see README), and a new picture:
// plots/cpp_ex6_fit.png
// =====================================================================

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TMath.h"
#include <iostream>

// This function has THREE arguments. The last two (fitLo, fitHi) are
// the edges of the fit range in GeV; they also have default values.
void ex6_fit(const char* filename = "../../flat_pipneta__B4_v6.root",
             double fitLo = 0.88, double fitHi = 1.10)
{
   gSystem->mkdir("plots", true);

   TFile* f = TFile::Open(filename);
   if (f == nullptr) { std::cout << "Could not open " << filename << std::endl; return; }
   TTree* tree = (TTree*)f->Get("pipneta__B4");
   tree->SetAlias("M_etapi",
      "sqrt(pow(PiPlus_E+Eta_E,2) - pow(PiPlus_Px+Eta_Px,2)"
      " - pow(PiPlus_Py+Eta_Py,2) - pow(PiPlus_Pz+Eta_Pz,2))");

   // ---- Make the histogram ----------------------------------------
   // 10 MeV bins this time, so the peak has enough bins to fit.
   TH1D* h = new TH1D("h",
         "a_{0}(980) #rightarrow #eta#pi^{+} fit;M(#eta#pi^{+}) (GeV/c^{2});Weighted entries / 10 MeV",
         180, 0.7, 2.5);
   h->Sumw2();
   tree->Draw("M_etapi>>h",
              "accidweight * (Chi2_NDF < 10 && -Mandelstam_t < 1.0 && abs(MissingMassSquared) < 0.05)",
              "goff");

   // ---- Define the model ------------------------------------------
   // A TF1 is a function of x that ROOT can fit. "gaus(0) + pol1(3)" means
   //   [0]*exp(-0.5*((x-[1])/[2])^2)   a Gaussian using parameters 0,1,2
   // + [3] + [4]*x                      a straight line using parameters 3,4
   // The [n] are the numbers the fit will adjust ("parameters").
   TF1* model = new TF1("model", "gaus(0) + pol1(3)", fitLo, fitHi);
   model->SetParNames("Amplitude", "Mean", "Sigma", "Bkg const", "Bkg slope");

   // A fit needs starting values close to the answer, or it can wander
   // off. Look at the plot from exercise 3 and guess:
   // peak height ~500 above background, at ~0.99 GeV, ~0.03 GeV wide,
   // background ~450 per 10 MeV and roughly flat.
   model->SetParameters(500, 0.99, 0.03, 450, 0);
   model->SetParLimits(2, 0.005, 0.2);   // keep sigma positive and sensible

   gStyle->SetOptFit(1111);   // show the fit results in a box on the plot
   gStyle->SetOptStat(0);
   TCanvas* c = new TCanvas("c", "c", 800, 600);

   // ---- Fit -------------------------------------------------------
   // "R" = only use the bins inside [fitLo, fitHi]
   // "S" = give back a result object so we can read the fit quality
   TFitResultPtr result = h->Fit(model, "RS");

   double mean  = model->GetParameter(1);
   double sigma = model->GetParameter(2);
   double chi2  = model->GetChisquare();
   int    ndf   = model->GetNDF();
   std::cout << "\n===== Summary =====" << std::endl;
   std::cout << "Mean     = " << 1000 * mean  << " +/- " << 1000 * model->GetParError(1) << " MeV" << std::endl;
   std::cout << "Sigma    = " << 1000 * sigma << " +/- " << 1000 * model->GetParError(2) << " MeV" << std::endl;
   // TMath::Prob gives the p-value: the chance of a chi2 this big or
   // bigger if the model were exactly right.
   std::cout << "chi2/ndf = " << chi2 << " / " << ndf
             << "  (p-value " << TMath::Prob(chi2, ndf) << ")" << std::endl;

   // Signal yield = area under the Gaussian / bin width.
   // The area of a Gaussian is amplitude * sigma * sqrt(2*pi).
   double area  = model->GetParameter(0) * sigma * TMath::Sqrt(2 * TMath::Pi());
   double yield = area / h->GetBinWidth(1);
   std::cout << "Signal yield ~ " << yield << " (weighted entries)" << std::endl;

   // ---- Draw ------------------------------------------------------
   // Draw the background part on its own, so we can see what is under the peak.
   TF1* bkg = new TF1("bkg", "pol1", fitLo, fitHi);
   bkg->SetParameters(model->GetParameter(3), model->GetParameter(4));
   bkg->SetLineStyle(2);
   bkg->SetLineColor(kGray + 2);

   h->GetXaxis()->SetRangeUser(0.7, 1.6);   // zoom in on the region of interest
   h->SetMinimum(0);
   h->SetMaximum(1.6 * h->GetMaximum());     // leave room for the fit box
   h->SetMarkerStyle(20);
   h->SetMarkerSize(0.7);
   h->Draw("e");
   bkg->Draw("same");

   TLegend* leg = new TLegend(0.12, 0.72, 0.45, 0.88);
   leg->AddEntry(h, "Data", "lep");
   leg->AddEntry(model, "Gaussian + line", "l");
   leg->AddEntry(bkg, "Background only", "l");
   leg->Draw();

   c->SaveAs("plots/cpp_ex6_fit.png");
   f->Close();
}
