// =====================================================================
// Exercise 3: Loop over the events, apply your own cuts, make a cut flow.
//
// What this macro does:
//   1. Connects each branch we need to a normal C++ variable.
//   2. Goes through the tree one entry at a time (an "event loop").
//   3. For each entry: works out the eta pi+ invariant mass from the
//      four-momenta, then applies three cuts one after another.
//   4. Counts how many entries survive each cut (the "cut flow").
//   5. Fills the eta pi+ mass histogram for entries that pass all cuts,
//      weighted by "accidweight" to subtract accidental beam photons.
//
// How to run it (in a terminal):
//   cd tutorial/cpp
//   root -l -b -q ex3_cutflow.C
// It reads all 6 million entries, so give it 10 to 30 seconds.
//
// What you should see:
//   A cut-flow table in the terminal (see README for the numbers) and
//   a new picture: plots/cpp_ex3_cutflow.png
//
// Try it: change the cut values in the "YOUR CUTS" block below, run
// again, and see how the table and the mass plot change.
// =====================================================================

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <iostream>
#include <iomanip>   // lets us line up columns when printing (std::setw)

void ex3_cutflow(const char* filename = "../../flat_pipneta__B4_v6.root")
{
   gSystem->mkdir("plots", true);

   // ---- YOUR CUTS: change these numbers and run again --------------
   // "double" is a variable type that holds a decimal number.
   double maxChi2NDF = 10.0;   // keep entries with Chi2_NDF below this
   double maxMinusT  = 1.0;    // keep entries with -t below this (GeV^2)
   double maxAbsMM2  = 0.05;   // keep entries with |missing mass^2| below this (GeV^2)

   TFile* f = TFile::Open(filename);
   if (f == nullptr) { std::cout << "Could not open " << filename << std::endl; return; }
   TTree* tree = (TTree*)f->Get("pipneta__B4");

   // ---- Connect branches to variables -----------------------------
   // First make ordinary variables, one per branch we want to read.
   // Their types must match the branch types printed in exercise 1
   // ("/D" = double).
   double piPx, piPy, piPz, piE;       // pi+ four-momentum (GeV)
   double etaPx, etaPy, etaPz, etaE;   // eta four-momentum (GeV)
   double chi2ndf, t, mm2, weight;

   // SetBranchAddress("branch name", &variable) tells the tree:
   // "every time you read an entry, put the value of this branch into
   // this variable". The "&" means "the address of" the variable, i.e.
   // where in memory to write the number.
   tree->SetBranchAddress("PiPlus_Px", &piPx);
   tree->SetBranchAddress("PiPlus_Py", &piPy);
   tree->SetBranchAddress("PiPlus_Pz", &piPz);
   tree->SetBranchAddress("PiPlus_E",  &piE);
   tree->SetBranchAddress("Eta_Px", &etaPx);
   tree->SetBranchAddress("Eta_Py", &etaPy);
   tree->SetBranchAddress("Eta_Pz", &etaPz);
   tree->SetBranchAddress("Eta_E",  &etaE);
   tree->SetBranchAddress("Chi2_NDF", &chi2ndf);
   tree->SetBranchAddress("Mandelstam_t", &t);
   tree->SetBranchAddress("MissingMassSquared", &mm2);
   tree->SetBranchAddress("accidweight", &weight);

   // Reading is faster if we switch off the branches we do not use:
   // "*" means "all branches", 0 means off, 1 means on.
   tree->SetBranchStatus("*", 0);
   const char* used[12] = {"PiPlus_Px", "PiPlus_Py", "PiPlus_Pz", "PiPlus_E",
                           "Eta_Px", "Eta_Py", "Eta_Pz", "Eta_E",
                           "Chi2_NDF", "Mandelstam_t", "MissingMassSquared", "accidweight"};
   for (int i = 0; i < 12; i++) tree->SetBranchStatus(used[i], 1);

   // ---- Book the histograms ---------------------------------------
   // The cut flow: one bin per step. Bin labels say which step.
   const int nSteps = 4;
   const char* stepNames[nSteps] = {"All entries",
                                    Form("#chi^{2}/ndf < %g", maxChi2NDF),
                                    Form("-t < %g GeV^{2}", maxMinusT),
                                    Form("|MM^{2}| < %g GeV^{2}", maxAbsMM2)};
   TH1D* hFlow = new TH1D("hFlow", "Cut flow (accidental-subtracted);;Weighted entries",
                          nSteps, 0, nSteps);

   // The eta pi+ mass after all cuts. 90 bins from 0.7 to 2.5 GeV = 20 MeV per bin.
   TH1D* hMass = new TH1D("hMass",
         "#eta#pi^{+} mass after all cuts;M(#eta#pi^{+}) (GeV/c^{2});Weighted entries / 20 MeV",
         90, 0.7, 2.5);
   // Sumw2() tells ROOT to keep track of the error bars correctly when
   // we fill with weights instead of just counting 1 per entry.
   hFlow->Sumw2();
   hMass->Sumw2();

   // Raw counters (how many entries, ignoring the weights).
   // "long" is a whole number that can be very big.
   long rawCount[nSteps] = {0, 0, 0, 0};

   // ---- The event loop --------------------------------------------
   long nEntries = tree->GetEntries();
   std::cout << "Looping over " << nEntries << " entries..." << std::endl;

   for (long i = 0; i < nEntries; i++) {
      // Read entry number i. After this line, all the variables we
      // connected above hold the values for this entry.
      tree->GetEntry(i);

      // Print a progress message every million entries.
      // "%" gives the remainder of a division: i % 1000000 is 0 only
      // when i is 0, 1000000, 2000000, ...
      if (i % 1000000 == 0) std::cout << "  entry " << i << std::endl;

      // Build four-vectors (px, py, pz, E). TLorentzVector knows how to
      // add them and how to work out a mass: M^2 = E^2 - |p|^2.
      TLorentzVector pip(piPx, piPy, piPz, piE);
      TLorentzVector eta(etaPx, etaPy, etaPz, etaE);
      double massEtaPi = (pip + eta).M();

      // Step 0: every entry counts.
      rawCount[0]++;                // "++" adds 1
      hFlow->Fill(0.5, weight);     // add "weight" to the bin that contains x = 0.5

      // Each cut: if the entry FAILS, "continue" skips the rest of the
      // loop and goes straight to the next entry.
      if (chi2ndf >= maxChi2NDF) continue;
      rawCount[1]++;
      hFlow->Fill(1.5, weight);

      if (-t >= maxMinusT) continue;   // t is negative, so -t is positive
      rawCount[2]++;
      hFlow->Fill(2.5, weight);

      if (TMath::Abs(mm2) >= maxAbsMM2) continue;   // Abs = absolute value
      rawCount[3]++;
      hFlow->Fill(3.5, weight);

      // Only entries that passed every cut get this far.
      hMass->Fill(massEtaPi, weight);
   }

   // ---- Print the cut-flow table ----------------------------------
   // std::setw(n) makes the next item take n characters, so columns line up.
   std::cout << "\n" << std::left << std::setw(26) << "Step"
             << std::right << std::setw(12) << "Entries"
             << std::setw(10) << "% kept"
             << std::setw(14) << "Weighted" << std::endl;
   // Plain names for the terminal (the "#chi" codes only work on plots).
   const char* plain[nSteps] = {"All entries", "Chi2_NDF cut", "-t cut", "|MM2| cut"};
   for (int s = 0; s < nSteps; s++) {
      hFlow->GetXaxis()->SetBinLabel(s + 1, stepNames[s]);   // label for the plot
      std::cout << std::left << std::setw(26) << plain[s]
                << std::right << std::setw(12) << rawCount[s]
                << std::setw(9) << std::fixed << std::setprecision(1)
                << 100.0 * rawCount[s] / rawCount[0] << "%"
                << std::setw(14) << std::setprecision(0) << hFlow->GetBinContent(s + 1)
                << std::endl;
   }
   std::cout << "\nWeighted entries in the mass plot: "
             << std::setprecision(0) << hMass->Integral() << std::endl;

   // ---- Draw ------------------------------------------------------
   gStyle->SetOptStat(0);
   TCanvas* c = new TCanvas("c", "c", 1200, 500);
   c->Divide(2, 1);

   c->cd(1);
   gPad->SetLogy();
   gPad->SetBottomMargin(0.15);
   hFlow->SetFillColor(kAzure - 9);
   hFlow->SetMinimum(1e4);
   hFlow->GetXaxis()->SetLabelSize(0.05);
   hFlow->Draw("hist");

   c->cd(2);
   hMass->SetMarkerStyle(20);
   hMass->SetMarkerSize(0.7);
   hMass->Draw("e");   // "e" = draw points with error bars

   c->SaveAs("plots/cpp_ex3_cutflow.png");
   f->Close();
}
