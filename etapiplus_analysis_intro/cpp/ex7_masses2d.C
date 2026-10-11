// =====================================================================
// Exercise 7: Two-dimensional mass plots (and a Dalitz plot).
//
// What this macro does:
//   For every entry that passes the standard cuts (chi2/ndf < 10,
//   -t < 1 GeV^2, |MM^2| < 0.05 GeV^2) it computes the three two-body
//   masses M(eta pi+), M(eta n), M(pi+ n) and fills, weighted by accidweight:
//     pad 1: M(eta pi+) vs M(eta n)
//     pad 2: M(eta pi+) vs M(pi+ n)
//     pad 3: M(eta n)   vs M(pi+ n)
//     pad 4: Dalitz plot, M^2(eta pi+) vs M^2(pi+ n)
//   and a second picture comparing M(eta pi+) with and without a
//   "baryon veto" (M(pi+ n) > 2 GeV and M(eta n) > 2 GeV).
//
// How to read the 2D plots:
//   A meson decaying to eta pi+ (a0(980), a2(1320)) has a fixed
//   M(eta pi+), so it makes a HORIZONTAL band (M(eta pi+) is on the y-axis).
//   A baryon decaying to pi+ n (Delta(1232)) or eta n (N*) makes a
//   VERTICAL band. Where bands cross, the two can interfere.
//
// How to run it (in a terminal):
//   cd etapiplus_analysis_intro/cpp
//   root -l -b -q ex7_masses2d.C
//
// What you should see: printed weighted sums (see README) and two
// pictures: plots/cpp_ex7_masses2d.png and plots/cpp_ex7_veto.png
// =====================================================================

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <iostream>

void ex7_masses2d(const char* filename = "../../flat_pipneta__B4_v6.root")
{
   gSystem->mkdir("plots", true);

   TFile* f = TFile::Open(filename);
   if (f == nullptr) { std::cout << "Could not open " << filename << std::endl; return; }
   TTree* tree = (TTree*)f->Get("pipneta__B4");

   // ---- Connect the branches we need to variables (see exercise 3) ----
   double piPx, piPy, piPz, piE, etaPx, etaPy, etaPz, etaE;
   double nPx, nPy, nPz, nE;
   double chi2ndf, t, mm2, weight;
   tree->SetBranchAddress("PiPlus_Px", &piPx);   tree->SetBranchAddress("PiPlus_Py", &piPy);
   tree->SetBranchAddress("PiPlus_Pz", &piPz);   tree->SetBranchAddress("PiPlus_E",  &piE);
   tree->SetBranchAddress("Eta_Px", &etaPx);     tree->SetBranchAddress("Eta_Py", &etaPy);
   tree->SetBranchAddress("Eta_Pz", &etaPz);     tree->SetBranchAddress("Eta_E",  &etaE);
   tree->SetBranchAddress("Neutron_Px", &nPx);   tree->SetBranchAddress("Neutron_Py", &nPy);
   tree->SetBranchAddress("Neutron_Pz", &nPz);   tree->SetBranchAddress("Neutron_E",  &nE);
   tree->SetBranchAddress("Chi2_NDF", &chi2ndf);
   tree->SetBranchAddress("Mandelstam_t", &t);
   tree->SetBranchAddress("MissingMassSquared", &mm2);
   tree->SetBranchAddress("accidweight", &weight);

   // ---- Book histograms ----
   // TH2D(name, "title;x label;y label", x bins, x low, x high, y bins, y low, y high)
   TH2D* hEtaPi_vs_EtaN = new TH2D("hEtaPi_vs_EtaN",
         "M(#eta#pi^{+}) vs M(#eta n);M(#eta n) (GeV/c^{2});M(#eta#pi^{+}) (GeV/c^{2})",
         100, 1.4, 4.4, 90, 0.7, 2.5);
   TH2D* hEtaPi_vs_PiN = new TH2D("hEtaPi_vs_PiN",
         "M(#eta#pi^{+}) vs M(#pi^{+}n);M(#pi^{+}n) (GeV/c^{2});M(#eta#pi^{+}) (GeV/c^{2})",
         100, 1.0, 4.0, 90, 0.7, 2.5);
   TH2D* hEtaN_vs_PiN = new TH2D("hEtaN_vs_PiN",
         "M(#eta n) vs M(#pi^{+}n);M(#pi^{+}n) (GeV/c^{2});M(#eta n) (GeV/c^{2})",
         100, 1.0, 4.0, 100, 1.4, 4.4);
   TH2D* hDalitz = new TH2D("hDalitz",
         "Dalitz plot;M^{2}(#pi^{+}n) (GeV^{2}/c^{4});M^{2}(#eta#pi^{+}) (GeV^{2}/c^{4})",
         100, 1.0, 16.0, 100, 0.0, 9.0);
   TH1D* hMass     = new TH1D("hMass", "Baryon veto;M(#eta#pi^{+}) (GeV/c^{2});Weighted entries / 20 MeV", 90, 0.7, 2.5);
   TH1D* hMassVeto = new TH1D("hMassVeto", "", 90, 0.7, 2.5);
   hMass->Sumw2();
   hMassVeto->Sumw2();

   // ---- Event loop ----
   long nEntries = tree->GetEntries();
   std::cout << "Looping over " << nEntries << " entries..." << std::endl;
   for (long i = 0; i < nEntries; i++) {
      tree->GetEntry(i);

      // Standard cuts (exercise 3). Failing entries are skipped.
      if (chi2ndf >= 10) continue;
      if (-t >= 1.0) continue;
      if (TMath::Abs(mm2) >= 0.05) continue;

      TLorentzVector pip(piPx, piPy, piPz, piE);
      TLorentzVector eta(etaPx, etaPy, etaPz, etaE);
      TLorentzVector neu(nPx, nPy, nPz, nE);
      double mEtaPi = (eta + pip).M();
      double mEtaN  = (eta + neu).M();
      double mPiN   = (pip + neu).M();

      // A 2D histogram is filled with Fill(x, y, weight).
      hEtaPi_vs_EtaN->Fill(mEtaN, mEtaPi, weight);
      hEtaPi_vs_PiN->Fill(mPiN, mEtaPi, weight);
      hEtaN_vs_PiN->Fill(mPiN, mEtaN, weight);
      // M2() is the mass squared.
      hDalitz->Fill((pip + neu).M2(), (eta + pip).M2(), weight);

      hMass->Fill(mEtaPi, weight);
      // Baryon veto: keep only entries far from the low-mass pi+ n and eta n regions.
      if (mPiN > 2.0 && mEtaN > 2.0) hMassVeto->Fill(mEtaPi, weight);
   }

   std::cout << "Weighted entries in the Dalitz plot:          " << hDalitz->Integral() << std::endl;
   std::cout << "M(eta pi+) 0.7-2.5 GeV, standard cuts:        " << hMass->Integral() << std::endl;
   std::cout << "M(eta pi+) 0.7-2.5 GeV, plus baryon veto:     " << hMassVeto->Integral() << std::endl;
   std::cout << "Fraction kept by the veto:                    " << hMassVeto->Integral() / hMass->Integral() << std::endl;

   // ---- Draw the four 2D plots ----
   gStyle->SetOptStat(0);
   TCanvas* c = new TCanvas("c", "c", 1200, 1000);
   c->Divide(2, 2);
   TH2D* all2d[4] = {hEtaPi_vs_EtaN, hEtaPi_vs_PiN, hEtaN_vs_PiN, hDalitz};
   for (int p = 0; p < 4; p++) {
      c->cd(p + 1);
      gPad->SetRightMargin(0.14);
      all2d[p]->SetMinimum(0.1);    // white = empty or negative after subtraction
      all2d[p]->Draw("colz");
   }
   c->SaveAs("plots/cpp_ex7_masses2d.png");

   // ---- Draw the effect of the baryon veto ----
   TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
   hMass->SetLineColor(kBlack);     hMass->SetLineWidth(2);
   hMassVeto->SetLineColor(kRed + 1); hMassVeto->SetLineWidth(2);
   hMass->SetMinimum(0);
   hMass->SetMaximum(1.3 * hMass->GetMaximum());
   hMass->Draw("hist");
   hMassVeto->Draw("hist same");
   TLegend* leg = new TLegend(0.15, 0.75, 0.6, 0.88);
   leg->AddEntry(hMass, "Standard cuts", "l");
   leg->AddEntry(hMassVeto, "+ M(#pi^{+}n) > 2 and M(#eta n) > 2 GeV", "l");
   leg->Draw();
   c2->SaveAs("plots/cpp_ex7_veto.png");
   f->Close();
}
