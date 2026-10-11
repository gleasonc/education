// =====================================================================
// Exercise 4: What does a cut do? And what are accidentals?
//
// What this macro does (two plots side by side):
//   left:  the eta pi+ mass with NO cuts and WITH the chi2/ndf cut,
//          both scaled to the same area so you can compare their SHAPES.
//   right: after the chi2/ndf cut, the mass for
//            - in-time ("prompt") beam photons,
//            - the estimate of accidental photons hiding in that sample,
//            - the difference = accidental-subtracted signal.
//
// New ideas: an "alias" (a nickname for a formula) and a weight in
// TTree::Draw.
//
// How to run it (in a terminal):
//   cd etapiplus_analysis_intro/cpp
//   root -l -b -q ex4_compare.C
//
// What you should see: printed sums for each curve (see README) and a
// new picture: plots/cpp_ex4_compare.png
// =====================================================================

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include <iostream>

void ex4_compare(const char* filename = "../../flat_pipneta__B4_v6.root")
{
   gSystem->mkdir("plots", true);

   TFile* f = TFile::Open(filename);
   if (f == nullptr) { std::cout << "Could not open " << filename << std::endl; return; }
   TTree* tree = (TTree*)f->Get("pipneta__B4");

   // ---- An alias: a nickname for a formula ------------------------
   // The tree has no "mass" branch, but it has the four-momenta. The
   // invariant mass of the eta pi+ pair is
   //    M = sqrt( (E1+E2)^2 - (px1+px2)^2 - (py1+py2)^2 - (pz1+pz2)^2 )
   // SetAlias lets us call that whole formula "M_etapi" from now on.
   tree->SetAlias("M_etapi",
      "sqrt(pow(PiPlus_E+Eta_E,2) - pow(PiPlus_Px+Eta_Px,2)"
      " - pow(PiPlus_Py+Eta_Py,2) - pow(PiPlus_Pz+Eta_Pz,2))");

   // Book the histograms: 90 bins from 0.7 to 2.5 GeV (20 MeV per bin).
   // A small helper "function" saves us typing the same line 5 times.
   // [&] and "auto" are C++ shorthand; you can read the next line as
   // "book(name, title) makes a TH1D with our usual binning".
   auto book = [&](const char* name, const char* title) {
      TH1D* h = new TH1D(name, title, 90, 0.7, 2.5);
      h->Sumw2();   // correct error bars when filling with weights
      return h;
   };
   const char* axes = ";M(#eta#pi^{+}) (GeV/c^{2})";
   TH1D* hNoCut   = book("hNoCut",   Form("Effect of the #chi^{2}/ndf cut%s;Fraction / 20 MeV", axes));
   TH1D* hCut     = book("hCut",     Form("Effect of the #chi^{2}/ndf cut%s;Fraction / 20 MeV", axes));
   TH1D* hPrompt  = book("hPrompt",  Form("Accidental subtraction (#chi^{2}/ndf < 10)%s;Entries / 20 MeV", axes));
   TH1D* hAccid   = book("hAccid",   "");
   TH1D* hSub     = book("hSub",     "");

   // ---- Fill with weights -----------------------------------------
   // The second argument of Draw can be a WEIGHT as well as a cut.
   // A true/false condition counts as 1 or 0, so
   //    "accidweight * (Chi2_NDF < 10)"
   // means "fill with weight accidweight, but only if Chi2_NDF < 10".
   tree->Draw("M_etapi>>hNoCut",  "accidweight", "goff");
   tree->Draw("M_etapi>>hCut",    "accidweight * (Chi2_NDF < 10)", "goff");

   // In-time photons only: accidweight == 1, each one counts as 1.
   tree->Draw("M_etapi>>hPrompt", "(accidweight == 1) * (Chi2_NDF < 10)", "goff");
   // Out-of-time photons: their weight is negative. Multiplying by -1
   // gives the (positive) number of accidentals expected under the prompt peak.
   tree->Draw("M_etapi>>hAccid",  "-accidweight * (accidweight < 0) * (Chi2_NDF < 10)", "goff");
   // Everything together with its own weight = prompt minus accidentals.
   tree->Draw("M_etapi>>hSub",    "accidweight * (Chi2_NDF < 10)", "goff");

   // Integral() adds up the bin contents (weights) inside the axis range.
   std::cout << "Weighted sum, no cut:              " << hNoCut->Integral()  << std::endl;
   std::cout << "Weighted sum, chi2/ndf < 10:       " << hCut->Integral()    << std::endl;
   std::cout << "Fraction kept by the cut:          " << hCut->Integral() / hNoCut->Integral() << std::endl;
   std::cout << "Prompt entries (chi2/ndf < 10):    " << hPrompt->Integral() << std::endl;
   std::cout << "Accidental estimate:               " << hAccid->Integral()  << std::endl;
   std::cout << "Prompt - accidental:               " << hSub->Integral()    << std::endl;

   // ---- Normalize: scale each histogram so its area is 1 ----------
   // Then shapes can be compared even though one has far more entries.
   hNoCut->Scale(1.0 / hNoCut->Integral());
   hCut->Scale(1.0 / hCut->Integral());

   // ---- Draw ------------------------------------------------------
   gStyle->SetOptStat(0);
   TCanvas* c = new TCanvas("c", "c", 1200, 500);
   c->Divide(2, 1);

   c->cd(1);
   gPad->SetLeftMargin(0.13);    // room for the y-axis label
   hNoCut->SetLineColor(kBlack);  hNoCut->SetLineWidth(2);
   hCut->SetLineColor(kRed + 1);  hCut->SetLineWidth(2);
   // Make sure the y-axis is tall enough for the taller of the two.
   hNoCut->SetMaximum(1.15 * TMath::Max(hNoCut->GetMaximum(), hCut->GetMaximum()));
   hNoCut->SetMinimum(0);
   hNoCut->Draw("hist");
   hCut->Draw("hist same");
   TLegend* leg1 = new TLegend(0.50, 0.15, 0.88, 0.30);
   leg1->SetTextSize(0.04);
   leg1->AddEntry(hNoCut, "No cut", "l");
   leg1->AddEntry(hCut, "#chi^{2}/ndf < 10", "l");
   leg1->Draw();

   c->cd(2);
   hPrompt->SetLineColor(kBlue + 1); hPrompt->SetLineWidth(2);
   hAccid->SetLineColor(kGray + 2);  hAccid->SetFillColor(kGray);
   hSub->SetLineColor(kRed + 1);     hSub->SetLineWidth(2);
   hPrompt->SetMinimum(0);
   hPrompt->SetMaximum(1.4 * hPrompt->GetMaximum());   // leave room for the legend
   hPrompt->Draw("hist");
   hAccid->Draw("hist same");
   hSub->Draw("hist same");
   TLegend* leg2 = new TLegend(0.40, 0.70, 0.88, 0.88);
   leg2->SetTextSize(0.04);
   leg2->AddEntry(hPrompt, "In-time (prompt)", "l");
   leg2->AddEntry(hAccid, "Accidental estimate", "f");
   leg2->AddEntry(hSub, "Prompt #minus accidental", "l");
   leg2->Draw();

   c->SaveAs("plots/cpp_ex4_compare.png");
   f->Close();
}
