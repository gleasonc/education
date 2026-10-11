// =====================================================================
// Exercise 8: Decay angles in the Gottfried-Jackson (GJ) frame.
//
// Why: a meson's spin shows up in how its decay products are spread in
// angle. A spin-0 a0(980) decays the same way in every direction (flat
// in cos theta); a spin-2 a2(1320) does not. To see this we look at the
// decay in the meson's own rest frame, with axes fixed by the reaction.
//
// The method (the same steps as the group's MakeFlat.C):
//   1. Lab -> centre-of-mass (CM) frame of beam + target. In the CM frame
//      the total momentum is zero, so the beam and the target come in
//      head-on and the eta pi+ system and the neutron go out back-to-back.
//   2. CM -> eta pi+ rest frame (the GJ frame). In this frame the eta
//      pi+ system is at rest, and the eta and pi+ fly apart back-to-back.
//   A "boost" changes the frame you look from: Boost(b) moves to a frame
//   travelling with velocity -b (in units of c). To go to the rest frame
//   of something with momentum p and energy E, use b = -p/E.
//   3. Axes of the GJ frame:
//        z = direction of the beam, seen in the eta pi+ rest frame
//        y = perpendicular to the production plane: beam x (eta pi+), in CM
//        x = y x z   (completes a right-handed set)
//   4. Take the eta direction in the GJ frame:
//        cos(theta_GJ) = its component along z
//        phi_GJ        = its angle around z, measured from x (degrees)
//
// Plots: cos(theta_GJ), phi_GJ, cos(theta_GJ) vs M(eta pi+), and
// cos(theta_GJ) in the a0(980) and a2(1320) mass windows.
// Cross-check: cos(theta_GJ) is also computed a second way, boosting
// straight from the lab to the eta pi+ rest frame and taking the angle
// between eta and beam. The two must agree for every entry.
//
// How to run it (in a terminal):
//   cd tutorial/cpp
//   root -l -b -q ex8_gj_angles.C
//
// What you should see: printed checks and sums (see README) and
// plots/cpp_ex8_gj_angles.png
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
#include "TVector3.h"
#include "TMath.h"
#include <iostream>

void ex8_gj_angles(const char* filename = "../../flat_pipneta__B4_v6.root")
{
   gSystem->mkdir("plots", true);

   TFile* f = TFile::Open(filename);
   if (f == nullptr) { std::cout << "Could not open " << filename << std::endl; return; }
   TTree* tree = (TTree*)f->Get("pipneta__B4");

   double piPx, piPy, piPz, piE, etaPx, etaPy, etaPz, etaE;
   double bPx, bPy, bPz, bE;
   double chi2ndf, t, mm2, weight;
   tree->SetBranchAddress("PiPlus_Px", &piPx);   tree->SetBranchAddress("PiPlus_Py", &piPy);
   tree->SetBranchAddress("PiPlus_Pz", &piPz);   tree->SetBranchAddress("PiPlus_E",  &piE);
   tree->SetBranchAddress("Eta_Px", &etaPx);     tree->SetBranchAddress("Eta_Py", &etaPy);
   tree->SetBranchAddress("Eta_Pz", &etaPz);     tree->SetBranchAddress("Eta_E",  &etaE);
   tree->SetBranchAddress("Beam_Px", &bPx);      tree->SetBranchAddress("Beam_Py", &bPy);
   tree->SetBranchAddress("Beam_Pz", &bPz);      tree->SetBranchAddress("Beam_E",  &bE);
   tree->SetBranchAddress("Chi2_NDF", &chi2ndf);
   tree->SetBranchAddress("Mandelstam_t", &t);
   tree->SetBranchAddress("MissingMassSquared", &mm2);
   tree->SetBranchAddress("accidweight", &weight);

   TH1D* hCos  = new TH1D("hCos", "cos#theta_{GJ} (all masses);cos#theta_{GJ};Weighted entries", 50, -1, 1);
   TH1D* hPhi  = new TH1D("hPhi", "#phi_{GJ} (all masses);#phi_{GJ} (deg);Weighted entries", 60, -180, 180);
   TH2D* hCosM = new TH2D("hCosM", "cos#theta_{GJ} vs M(#eta#pi^{+});M(#eta#pi^{+}) (GeV/c^{2});cos#theta_{GJ}",
                          45, 0.7, 2.5, 40, -1, 1);
   TH1D* hCosA0 = new TH1D("hCosA0", "cos#theta_{GJ} in mass windows;cos#theta_{GJ};Weighted entries", 25, -1, 1);
   TH1D* hCosA2 = new TH1D("hCosA2", "", 25, -1, 1);
   hCos->Sumw2(); hPhi->Sumw2(); hCosA0->Sumw2(); hCosA2->Sumw2();

   // Mass windows (GeV)
   double a0Lo = 0.92, a0Hi = 1.04;
   double a2Lo = 1.24, a2Hi = 1.40;

   // The target proton is at rest in the lab: (px, py, pz, E) = (0, 0, 0, m_p)
   TLorentzVector target(0, 0, 0, 0.938272);

   double cosMin = 2, cosMax = -2, phiMin = 999, phiMax = -999, maxDiff = 0;
   long nUsed = 0;

   long nEntries = tree->GetEntries();
   std::cout << "Looping over " << nEntries << " entries..." << std::endl;
   for (long i = 0; i < nEntries; i++) {
      tree->GetEntry(i);
      if (chi2ndf >= 10) continue;
      if (-t >= 1.0) continue;
      if (TMath::Abs(mm2) >= 0.05) continue;

      TLorentzVector beam(bPx, bPy, bPz, bE);
      TLorentzVector pip(piPx, piPy, piPz, piE);
      TLorentzVector eta(etaPx, etaPy, etaPz, etaE);
      TLorentzVector etaPi = eta + pip;
      double mEtaPi = etaPi.M();

      // ---- Step 1: lab -> CM frame ----
      TLorentzVector cm = beam + target;
      TVector3 boostCM = -(cm.Vect()) * (1.0 / cm.E());   // -p/E of the whole system
      TLorentzVector beamCM = beam;   beamCM.Boost(boostCM);
      TLorentzVector etaCM  = eta;    etaCM.Boost(boostCM);
      TLorentzVector etaPiCM = etaPi; etaPiCM.Boost(boostCM);

      // ---- Step 2: CM -> eta pi+ rest frame ----
      TVector3 boostGJ = -(etaPiCM.Vect()) * (1.0 / etaPiCM.E());
      TLorentzVector beamGJ = beamCM; beamGJ.Boost(boostGJ);
      TLorentzVector etaGJ  = etaCM;  etaGJ.Boost(boostGJ);

      // ---- Step 3: the GJ axes ----
      TVector3 zGJ = beamGJ.Vect().Unit();
      TVector3 yGJ = beamCM.Vect().Cross(etaPiCM.Vect()).Unit();
      TVector3 xGJ = yGJ.Cross(zGJ);

      // ---- Step 4: eta direction in those axes ----
      // a.Dot(b) (or a*b) is the component of a along b when b has length 1.
      TVector3 v(etaGJ.Vect().Dot(xGJ), etaGJ.Vect().Dot(yGJ), etaGJ.Vect().Dot(zGJ));
      double cosTheta = v.CosTheta();
      double phi = v.Phi() * 180.0 / TMath::Pi();

      // ---- Cross-check: lab -> eta pi+ rest frame in one boost ----
      // BoostVector() is p/E, so -BoostVector() takes us to the rest frame.
      TLorentzVector beam2 = beam; beam2.Boost(-etaPi.BoostVector());
      TLorentzVector eta2  = eta;  eta2.Boost(-etaPi.BoostVector());
      double cosCheck = TMath::Cos(eta2.Angle(beam2.Vect()));   // Angle() = angle between two 3-vectors
      maxDiff = TMath::Max(maxDiff, TMath::Abs(cosTheta - cosCheck));

      cosMin = TMath::Min(cosMin, cosTheta);  cosMax = TMath::Max(cosMax, cosTheta);
      phiMin = TMath::Min(phiMin, phi);       phiMax = TMath::Max(phiMax, phi);
      nUsed++;

      hCos->Fill(cosTheta, weight);
      hPhi->Fill(phi, weight);
      hCosM->Fill(mEtaPi, cosTheta, weight);
      if (mEtaPi > a0Lo && mEtaPi < a0Hi) hCosA0->Fill(cosTheta, weight);
      if (mEtaPi > a2Lo && mEtaPi < a2Hi) hCosA2->Fill(cosTheta, weight);
   }

   std::cout << "Entries used:                    " << nUsed << std::endl;
   std::cout << "cos(theta_GJ) range:             " << cosMin << " to " << cosMax << std::endl;
   std::cout << "phi_GJ range (deg):              " << phiMin << " to " << phiMax << std::endl;
   std::cout << "Max |cos - cross-check|:         " << maxDiff << std::endl;
   std::cout << "Weighted entries, all masses:    " << hCos->Integral() << std::endl;
   std::cout << "Weighted entries, a0 window:     " << hCosA0->Integral()
             << "  (mean cos = " << hCosA0->GetMean() << ")" << std::endl;
   std::cout << "Weighted entries, a2 window:     " << hCosA2->Integral()
             << "  (mean cos = " << hCosA2->GetMean() << ")" << std::endl;

   // ---- Draw ----
   gStyle->SetOptStat(0);
   TCanvas* c = new TCanvas("c", "c", 1200, 1000);
   c->Divide(2, 2);

   c->cd(1);
   hCos->SetMinimum(0);
   hCos->SetLineWidth(2);
   hCos->Draw("hist");

   c->cd(2);
   hPhi->SetMinimum(0);
   hPhi->SetLineWidth(2);
   hPhi->Draw("hist");

   c->cd(3);
   gPad->SetRightMargin(0.14);
   hCosM->SetMinimum(0.1);
   hCosM->Draw("colz");

   c->cd(4);
   // Scale both to area 1 to compare shapes.
   hCosA0->Scale(1.0 / hCosA0->Integral());
   hCosA2->Scale(1.0 / hCosA2->Integral());
   hCosA0->GetYaxis()->SetTitle("Fraction / 0.08");
   hCosA0->SetLineColor(kBlue + 1); hCosA0->SetLineWidth(2);
   hCosA2->SetLineColor(kRed + 1);  hCosA2->SetLineWidth(2);
   hCosA0->SetMinimum(0);
   hCosA0->SetMaximum(1.5 * TMath::Max(hCosA0->GetMaximum(), hCosA2->GetMaximum()));
   hCosA0->Draw("hist");
   hCosA2->Draw("hist same");
   TLegend* leg = new TLegend(0.15, 0.75, 0.6, 0.88);
   leg->AddEntry(hCosA0, "a_{0}(980): 0.92-1.04 GeV", "l");
   leg->AddEntry(hCosA2, "a_{2}(1320): 1.24-1.40 GeV", "l");
   leg->Draw();

   c->SaveAs("plots/cpp_ex8_gj_angles.png");
   f->Close();
}
