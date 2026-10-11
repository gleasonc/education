// =====================================================================
// Exercise 2: Draw labelled histograms straight from the tree.
//
// What this macro does:
//   Uses TTree::Draw, the quickest way to histogram one branch, to make
//   two plots side by side:
//     left:  the photon beam energy (Beam_E)
//     right: the kinematic-fit quality (Chi2_NDF), all entries and
//            only the entries with Chi2_NDF < 10
//   Every plot gets a title and labelled axes with units.
//
// How to run it (in a terminal):
//   cd etapiplus_analysis_intro/cpp
//   root -l -b -q ex2_draw.C
//
// What you should see:
//   A few lines of printout with the number of entries in each histogram,
//   and a new picture: plots/cpp_ex2_draw.png
//   (Leave out "-b -q" if you want the plot to pop up on screen.)
// =====================================================================

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include <iostream>

void ex2_draw(const char* filename = "../../flat_pipneta__B4_v6.root")
{
   // Make a folder called "plots" to save pictures in. "true" means
   // "it is fine if the folder already exists".
   gSystem->mkdir("plots", true);

   // Open the file and get the tree (see exercise 1 for what each part
   // means: "*" = pointer, "->" = use what the pointer points to).
   TFile* f = TFile::Open(filename);
   if (f == nullptr) { std::cout << "Could not open " << filename << std::endl; return; }
   TTree* tree = (TTree*)f->Get("pipneta__B4");

   // ---- Make the empty histograms first ---------------------------
   // new TH1D(name, "title;x-axis label;y-axis label", number of bins, lowest x, highest x)
   // The text between the ";" signs sets the title and the axis labels.
   // "new" creates the object and gives back a pointer to it.
   TH1D* hBeam = new TH1D("hBeam",
         "Photon beam energy;E_{#gamma} (GeV);Entries / 20 MeV", 200, 8.0, 12.0);
   TH1D* hChi2All = new TH1D("hChi2All",
         "Kinematic-fit #chi^{2}/ndf;#chi^{2}/ndf;Entries / 0.5", 100, 0.0, 50.0);
   TH1D* hChi2Cut = new TH1D("hChi2Cut",
         "Kinematic-fit #chi^{2}/ndf;#chi^{2}/ndf;Entries / 0.5", 100, 0.0, 50.0);

   // ---- Fill them with TTree::Draw --------------------------------
   // tree->Draw("what to plot >> which histogram", "cut", "options")
   //   "Beam_E>>hBeam"  : take the branch Beam_E and fill histogram hBeam
   //   ""               : no cut, use every entry
   //   "goff"           : "graphics off" - fill only, do not draw yet
   tree->Draw("Beam_E>>hBeam", "", "goff");
   tree->Draw("Chi2_NDF>>hChi2All", "", "goff");
   // The second argument is a CUT: only entries where it is true are used.
   tree->Draw("Chi2_NDF>>hChi2Cut", "Chi2_NDF < 10", "goff");

   // "(long)" turns the number into a whole number so it prints as 6091343, not 6.09e+06.
   std::cout << "Entries in hBeam:    " << (long)hBeam->GetEntries() << std::endl;
   std::cout << "Entries in hChi2All: " << (long)hChi2All->GetEntries()
             << "  (inside 0-50: " << (long)hChi2All->Integral() << ")" << std::endl;
   std::cout << "Entries in hChi2Cut: " << (long)hChi2Cut->GetEntries() << std::endl;

   // ---- Draw them -------------------------------------------------
   gStyle->SetOptStat(0);        // hide the statistics box on the plots

   // A canvas is the window/picture we draw on. Divide(2,1) splits it
   // into 2 columns and 1 row: two "pads" side by side.
   TCanvas* c = new TCanvas("c", "c", 1200, 500);
   c->Divide(2, 1);

   // Left pad: beam energy
   c->cd(1);                     // "cd" = go to pad number 1
   hBeam->SetLineColor(kBlue + 1);
   hBeam->SetLineWidth(2);
   hBeam->Draw("hist");          // "hist" = draw as a histogram outline

   // Right pad: chi2/ndf, both versions on top of each other
   c->cd(2);
   gPad->SetLogy();              // logarithmic y-axis: shows big and small counts together
   hChi2All->SetLineColor(kBlack);
   hChi2All->SetLineWidth(2);
   hChi2Cut->SetLineColor(kRed + 1);
   hChi2Cut->SetFillColorAlpha(kRed - 9, 0.5);
   hChi2All->Draw("hist");
   hChi2Cut->Draw("hist same");  // "same" = draw on top of what is already there

   // A legend says which line is which.
   TLegend* leg = new TLegend(0.50, 0.75, 0.88, 0.88);
   leg->SetTextSize(0.04);
   leg->AddEntry(hChi2All, "All entries", "l");
   leg->AddEntry(hChi2Cut, "#chi^{2}/ndf < 10", "f");
   leg->Draw();

   c->SaveAs("plots/cpp_ex2_draw.png");
}
