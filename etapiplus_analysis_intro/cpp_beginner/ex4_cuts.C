// =====================================================================
// Beginner exercise 4: Apply cuts and compare before and after.
//
// Same loop as exercise 3, plus three cuts. We count how many entries
// pass each cut, and fill the eta pi+ mass twice: once before the cuts
// and once after. Both are drawn on one plot.
//
// HOW TO RUN:
//     cd etapiplus_analysis_intro/cpp_beginner
//     root -l ex4_cuts.C
//
// WHAT YOU SHOULD SEE:
//     All entries:   6091343
//     Pass chi2 cut: 699087
//     Pass -t cut:   433403
//     Pass MM2 cut:  362187
// and ex4_cuts.pdf (and .png): black = before cuts, red = after cuts.
// After the cuts there are clear peaks near 1.0 GeV (a0(980)) and
// 1.3 GeV (a2(1320)).
//
// TRY: change the numbers in the "if" lines and run again.
// =====================================================================

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLorentzVector.h"

void ex4_cuts()
{
    // ---- 1. Open the file and get the tree ----
    TFile *fileIn=new TFile("../../flat_pipneta__B4_v6.root");
    TTree *tree=(TTree*)fileIn->Get("pipneta__B4");

    // ---- 2. Variables for the branches we need ----
    double PiPlus_Px, PiPlus_Py, PiPlus_Pz, PiPlus_E;
    double Eta_Px, Eta_Py, Eta_Pz, Eta_E;
    double Chi2_NDF, Mandelstam_t, MissingMassSquared, accidweight;

    tree->SetBranchAddress("PiPlus_Px", &PiPlus_Px );
    tree->SetBranchAddress("PiPlus_Py", &PiPlus_Py );
    tree->SetBranchAddress("PiPlus_Pz", &PiPlus_Pz );
    tree->SetBranchAddress("PiPlus_E", &PiPlus_E );

    tree->SetBranchAddress("Eta_Px", &Eta_Px );
    tree->SetBranchAddress("Eta_Py", &Eta_Py );
    tree->SetBranchAddress("Eta_Pz", &Eta_Pz );
    tree->SetBranchAddress("Eta_E", &Eta_E );

    tree->SetBranchAddress("Chi2_NDF", &Chi2_NDF );
    tree->SetBranchAddress("Mandelstam_t", &Mandelstam_t );
    tree->SetBranchAddress("MissingMassSquared", &MissingMassSquared );
    tree->SetBranchAddress("accidweight", &accidweight );

    // ---- 3. Empty histograms ----
    TH1F* h_EtaPiPlusMass_NoCuts=new TH1F("h_EtaPiPlusMass_NoCuts",";M(#eta#pi^{+}) (GeV);Entries",90,0.7,2.5);
    TH1F* h_EtaPiPlusMass_Cuts=new TH1F("h_EtaPiPlusMass_Cuts",";M(#eta#pi^{+}) (GeV);Entries",90,0.7,2.5);

    // Counters: whole numbers that start at 0.
    int nAll=0;
    int nPassChi2=0;
    int nPassT=0;
    int nPassMM2=0;

    // ---- 4. Loop over every entry ----
    Long64_t NEntries=tree->GetEntries();
    cout<<"There are "<<NEntries<<" entries."<<endl;

    for (Long64_t i_entry=0; i_entry<NEntries; i_entry++) {
        tree->GetEntry(i_entry);
        if (i_entry%1000000==0)
            cout<<i_entry<<" events processed"<<endl;

        TLorentzVector pipP4;
        pipP4.SetPxPyPzE(PiPlus_Px, PiPlus_Py, PiPlus_Pz, PiPlus_E);
        TLorentzVector etaP4;
        etaP4.SetPxPyPzE(Eta_Px, Eta_Py, Eta_Pz, Eta_E);

        TLorentzVector EtaPiPlusP4=etaP4+pipP4;
        double EtaPiPlusMass=EtaPiPlusP4.M();

        // Before any cuts
        nAll=nAll+1;
        h_EtaPiPlusMass_NoCuts->Fill(EtaPiPlusMass, accidweight);

        // ---- Cuts ----
        // "continue" means: skip the rest of the loop for this entry and
        // go on to the next one. So each "if" throws some entries away.

        // Cut 1: a good kinematic fit
        if (Chi2_NDF>10)
            continue;
        nPassChi2=nPassChi2+1;

        // Cut 2: small momentum transfer. t is negative, so we use -t.
        if (-Mandelstam_t>1)
            continue;
        nPassT=nPassT+1;

        // Cut 3: missing mass squared close to zero, on either side
        if (MissingMassSquared>0.05)
            continue;
        if (MissingMassSquared<-0.05)
            continue;
        nPassMM2=nPassMM2+1;

        // Only entries that passed every cut get here.
        h_EtaPiPlusMass_Cuts->Fill(EtaPiPlusMass, accidweight);
    }

    // ---- Print the counts ----
    cout<<"All entries:   "<<nAll<<endl;
    cout<<"Pass chi2 cut: "<<nPassChi2<<endl;
    cout<<"Pass -t cut:   "<<nPassT<<endl;
    cout<<"Pass MM2 cut:  "<<nPassMM2<<endl;

    // Hide the statistics box (it would cover the plots).
    gStyle->SetOptStat(0);

    // ---- 5. Draw both on one plot and save ----
    TCanvas *c1=new TCanvas("c1","",800,600);

    h_EtaPiPlusMass_NoCuts->SetMarkerStyle(20);
    h_EtaPiPlusMass_NoCuts->SetMarkerColor(kBlack);
    h_EtaPiPlusMass_NoCuts->SetLineColor(kBlack);
    h_EtaPiPlusMass_NoCuts->SetMarkerSize(1);
    h_EtaPiPlusMass_NoCuts->SetMinimum(0);
    h_EtaPiPlusMass_NoCuts->SetMaximum(13000);    // leave room for the legend
    h_EtaPiPlusMass_NoCuts->Draw();

    h_EtaPiPlusMass_Cuts->SetMarkerStyle(20);
    h_EtaPiPlusMass_Cuts->SetMarkerColor(kRed);
    h_EtaPiPlusMass_Cuts->SetLineColor(kRed);
    h_EtaPiPlusMass_Cuts->SetMarkerSize(1);
    h_EtaPiPlusMass_Cuts->Draw("same");      // "same" = draw on top

    TLegend *leg=new TLegend(0.15,0.75,0.5,0.88);
    leg->AddEntry(h_EtaPiPlusMass_NoCuts,"Before cuts","p");
    leg->AddEntry(h_EtaPiPlusMass_Cuts,"After cuts","p");
    leg->Draw();

    c1->SaveAs("ex4_cuts.pdf");
    c1->SaveAs("ex4_cuts.png");
}
