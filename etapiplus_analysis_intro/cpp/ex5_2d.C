// =====================================================================
// Exercise 5: A 2D plot and its projections.
//
// What this macro does:
//   left:  a 2D histogram of the eta pi+ mass (x) against -t (y), the
//          momentum transfer to the neutron. Colour = number of entries.
//   right: two "projections" of that 2D plot onto the mass axis, one
//          for small -t and one for larger -t, drawn on top of each other.
//
// How to run it (in a terminal):
//   cd tutorial/cpp
//   root -l -b -q ex5_2d.C
//
// What you should see: the weighted sum in each slice printed (see
// README) and a new picture: plots/cpp_ex5_2d.png
// =====================================================================

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include <iostream>

void ex5_2d(const char* filename = "../../flat_pipneta__B4_v6.root")
{
   gSystem->mkdir("plots", true);

   TFile* f = TFile::Open(filename);
   if (f == nullptr) { std::cout << "Could not open " << filename << std::endl; return; }
   TTree* tree = (TTree*)f->Get("pipneta__B4");

   // The same alias as in exercise 4: a nickname for the mass formula.
   tree->SetAlias("M_etapi",
      "sqrt(pow(PiPlus_E+Eta_E,2) - pow(PiPlus_Px+Eta_Px,2)"
      " - pow(PiPlus_Py+Eta_Py,2) - pow(PiPlus_Pz+Eta_Pz,2))");

   // ---- Book a 2D histogram ---------------------------------------
   // new TH2D(name, "title;x label;y label",
   //          x bins, x low, x high,   y bins, y low, y high)
   TH2D* h2 = new TH2D("h2",
         "#eta#pi^{+} mass vs -t;M(#eta#pi^{+}) (GeV/c^{2});-t (GeV^{2})",
         90, 0.7, 2.5,    40, 0.0, 2.0);
   h2->Sumw2();

   // ---- Fill it ---------------------------------------------------
   // For a 2D plot, Draw takes "y:x" (y first, then x!).
   // -Mandelstam_t is -t, which is a positive number.
   // Weight = accidweight, cuts = the same ones as in exercise 3.
   tree->Draw("-Mandelstam_t:M_etapi>>h2",
              "accidweight * (Chi2_NDF < 10 && abs(MissingMassSquared) < 0.05)", "goff");

   // ---- Projections -----------------------------------------------
   // ProjectionX(name, first y bin, last y bin) adds up all the y bins
   // in that range and gives a normal 1D histogram of x.
   // FindBin(value) tells us which bin a value falls in.
   // Small -t slice: 0 to 0.5 GeV^2. The "- 0.001" makes sure we take
   // the bin BELOW 0.5, not the one that starts at 0.5.
   int binLo1 = h2->GetYaxis()->FindBin(0.0);
   int binHi1 = h2->GetYaxis()->FindBin(0.5 - 0.001);
   int binLo2 = h2->GetYaxis()->FindBin(0.5);
   int binHi2 = h2->GetYaxis()->FindBin(1.0 - 0.001);
   TH1D* pLow  = h2->ProjectionX("pLow",  binLo1, binHi1);
   TH1D* pHigh = h2->ProjectionX("pHigh", binLo2, binHi2);

   std::cout << "Weighted entries in the 2D plot:     " << h2->Integral()    << std::endl;
   std::cout << "Weighted entries with 0 < -t < 0.5:   " << pLow->Integral()  << std::endl;
   std::cout << "Weighted entries with 0.5 < -t < 1.0: " << pHigh->Integral() << std::endl;

   // ---- Draw ------------------------------------------------------
   gStyle->SetOptStat(0);
   TCanvas* c = new TCanvas("c", "c", 1200, 500);
   c->Divide(2, 1);

   c->cd(1);
   gPad->SetRightMargin(0.14);   // room for the colour scale
   h2->SetMinimum(0);            // do not colour bins below 0 (see README)
   h2->Draw("colz");             // "colz" = colours plus a colour scale (z axis)

   c->cd(2);
   pLow->SetTitle("Projections onto the mass axis;M(#eta#pi^{+}) (GeV/c^{2});Weighted entries / 20 MeV");
   pLow->SetLineColor(kBlue + 1);  pLow->SetLineWidth(2);
   pHigh->SetLineColor(kRed + 1);  pHigh->SetLineWidth(2);
   pLow->SetMinimum(0);
   pLow->SetMaximum(1.3 * TMath::Max(pLow->GetMaximum(), pHigh->GetMaximum()));
   pLow->Draw("hist");
   pHigh->Draw("hist same");
   TLegend* leg = new TLegend(0.50, 0.75, 0.88, 0.88);
   leg->SetTextSize(0.04);
   leg->AddEntry(pLow,  "0 < -t < 0.5 GeV^{2}", "l");
   leg->AddEntry(pHigh, "0.5 < -t < 1.0 GeV^{2}", "l");
   leg->Draw();

   c->SaveAs("plots/cpp_ex5_2d.png");
   f->Close();
}
