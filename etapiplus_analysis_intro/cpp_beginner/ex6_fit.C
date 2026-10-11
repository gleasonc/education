// =====================================================================
// Beginner exercise 6: Fit the a0(980) peak.
//
// We fill the eta pi+ mass with the same cuts as exercise 3, then fit
// the peak near 1 GeV with a Gaussian (the peak) plus a straight line
// (the background under it).
//
// HOW TO RUN:
//     cd tutorial/cpp_beginner
//     root -l ex6_fit.C
//
// WHAT YOU SHOULD SEE: a table of fit results from ROOT. The ones to
// read are
//     Mean  = about 0.986 GeV
//     Sigma = about 0.026 GeV
// and a plot saved as ex6_fit.pdf (and .png).
//
// TRY: change the fit range (fitLow, fitHigh) and see how much the
// results move.
// =====================================================================

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLorentzVector.h"

void ex6_fit()
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

    // ---- 3. Empty histogram (10 MeV bins) ----
    TH1F* h_EtaPiPlusMass=new TH1F("h_EtaPiPlusMass","a_{0}(980) fit;M(#eta#pi^{+}) (GeV);Entries / 10 MeV",180,0.7,2.5);

    // Keep track of the error bars correctly when filling with weights.
    h_EtaPiPlusMass->Sumw2();

    // ---- 4. Loop over every entry (same cuts as exercise 3) ----
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
        double EtaPiPlusMass=(etaP4+pipP4).M();

        if (Chi2_NDF>10)
            continue;
        if (-Mandelstam_t>1)
            continue;
        if (MissingMassSquared>0.05)
            continue;
        if (MissingMassSquared<-0.05)
            continue;

        h_EtaPiPlusMass->Fill(EtaPiPlusMass, accidweight);
    }

    // ---- The fit ----
    double fitLow=0.88;      // fit range in GeV
    double fitHigh=1.10;

    // The fit function. "gaus(0)" is a Gaussian that uses parameters 0, 1, 2
    // (height, mean, width). "pol1(3)" is a straight line that uses
    // parameters 3 and 4 (constant, slope).
    TF1 *fitFunc=new TF1("fitFunc","gaus(0)+pol1(3)",fitLow,fitHigh);

    // Starting guesses, read off the plot from exercise 3:
    // height 500, mean 0.99 GeV, width 0.03 GeV, background 450, slope 0.
    fitFunc->SetParameters(500, 0.99, 0.03, 450, 0);
    fitFunc->SetParNames("Height","Mean","Sigma","Bkg const","Bkg slope");

    gStyle->SetOptFit(1111);    // show the fit results on the plot

    TCanvas *c1=new TCanvas("c1","",800,600);

    // "R" = only fit inside fitLow to fitHigh.
    h_EtaPiPlusMass->Fit(fitFunc,"R");

    cout<<"Mean  = "<<fitFunc->GetParameter(1)<<" +/- "<<fitFunc->GetParError(1)<<" GeV"<<endl;
    cout<<"Sigma = "<<fitFunc->GetParameter(2)<<" +/- "<<fitFunc->GetParError(2)<<" GeV"<<endl;
    cout<<"chi2/ndf = "<<fitFunc->GetChisquare()<<" / "<<fitFunc->GetNDF()<<endl;

    // ---- 5. Draw and save ----
    h_EtaPiPlusMass->GetXaxis()->SetRangeUser(0.7,1.6);   // zoom in
    h_EtaPiPlusMass->SetMinimum(0);
    h_EtaPiPlusMass->SetMaximum(1600);
    h_EtaPiPlusMass->Draw();
    c1->SaveAs("ex6_fit.pdf");
    c1->SaveAs("ex6_fit.png");
}
