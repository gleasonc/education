// =====================================================================
// Beginner exercise 5: 2D histograms.
//
// Three plots, all after the chi2 and missing-mass cuts:
//     1. M(eta pi+) against M(pi+ n): a "Dalitz-style" plot. A meson
//        decaying to eta pi+ makes a vertical band; a Delta or N* decaying
//        to pi+ n would make a horizontal band.
//     2. M(eta pi+) against -t.
//     3. Two slices of plot 2 projected onto the mass axis:
//        blue = 0 < -t < 0.5 GeV^2, red = 0.5 < -t < 1.0 GeV^2.
//
// HOW TO RUN:
//     cd etapiplus_analysis_intro/cpp_beginner
//     root -l ex5_2d.C
//
// WHAT YOU SHOULD SEE:
//     Entries with 0 < -t < 0.5:   32873.1
//     Entries with 0.5 < -t < 1.0: 43234.7
// and ex5_2d.pdf (and .png). Bright vertical stripes near 1.0 and 1.3 GeV
// are the a0(980) and a2(1320).
// =====================================================================

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLorentzVector.h"

void ex5_2d()
{
    // ---- 1. Open the file and get the tree ----
    TFile *fileIn=new TFile("../../flat_pipneta__B4_v6.root");
    TTree *tree=(TTree*)fileIn->Get("pipneta__B4");

    // ---- 2. Variables for the branches we need ----
    double PiPlus_Px, PiPlus_Py, PiPlus_Pz, PiPlus_E;
    double Eta_Px, Eta_Py, Eta_Pz, Eta_E;
    double Neutron_Px, Neutron_Py, Neutron_Pz, Neutron_E;
    double Chi2_NDF, Mandelstam_t, MissingMassSquared, accidweight;

    tree->SetBranchAddress("PiPlus_Px", &PiPlus_Px );
    tree->SetBranchAddress("PiPlus_Py", &PiPlus_Py );
    tree->SetBranchAddress("PiPlus_Pz", &PiPlus_Pz );
    tree->SetBranchAddress("PiPlus_E", &PiPlus_E );

    tree->SetBranchAddress("Eta_Px", &Eta_Px );
    tree->SetBranchAddress("Eta_Py", &Eta_Py );
    tree->SetBranchAddress("Eta_Pz", &Eta_Pz );
    tree->SetBranchAddress("Eta_E", &Eta_E );

    tree->SetBranchAddress("Neutron_Px", &Neutron_Px );
    tree->SetBranchAddress("Neutron_Py", &Neutron_Py );
    tree->SetBranchAddress("Neutron_Pz", &Neutron_Pz );
    tree->SetBranchAddress("Neutron_E", &Neutron_E );

    tree->SetBranchAddress("Chi2_NDF", &Chi2_NDF );
    tree->SetBranchAddress("Mandelstam_t", &Mandelstam_t );
    tree->SetBranchAddress("MissingMassSquared", &MissingMassSquared );
    tree->SetBranchAddress("accidweight", &accidweight );

    // ---- 3. Empty 2D histograms ----
    // new TH2F("name", ";x label;y label", x bins, x low, x high, y bins, y low, y high)
    TH2F* h_PiPlusNeutronMass_vs_EtaPiPlusMass=new TH2F("h_PiPlusNeutronMass_vs_EtaPiPlusMass",";M(#eta#pi^{+}) (GeV);M(#pi^{+}n) (GeV)",90,0.7,2.5,100,1.0,4.0);
    TH2F* h_t_vs_EtaPiPlusMass=new TH2F("h_t_vs_EtaPiPlusMass",";M(#eta#pi^{+}) (GeV);-t (GeV^{2})",90,0.7,2.5,40,0,2);

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
        TLorentzVector neutronP4;
        neutronP4.SetPxPyPzE(Neutron_Px, Neutron_Py, Neutron_Pz, Neutron_E);

        TLorentzVector EtaPiPlusP4=etaP4+pipP4;
        double EtaPiPlusMass=EtaPiPlusP4.M();

        TLorentzVector PiPlusNeutronP4=pipP4+neutronP4;
        double PiPlusNeutronMass=PiPlusNeutronP4.M();

        // Cuts (no -t cut here, because we want to see -t)
        if (Chi2_NDF>10)
            continue;
        if (MissingMassSquared>0.05)
            continue;
        if (MissingMassSquared<-0.05)
            continue;

        // A 2D histogram is filled with (x, y, weight).
        h_PiPlusNeutronMass_vs_EtaPiPlusMass->Fill(EtaPiPlusMass, PiPlusNeutronMass, accidweight);
        h_t_vs_EtaPiPlusMass->Fill(EtaPiPlusMass, -Mandelstam_t, accidweight);
    }

    // ---- Projections ----
    // ProjectionX adds up a range of y bins and gives a 1D histogram.
    // The -t axis has 40 bins from 0 to 2, so each bin is 0.05 wide:
    //     bins 1 to 10  ->  0   < -t < 0.5
    //     bins 11 to 20 ->  0.5 < -t < 1.0
    TH1D* h_Mass_LowT=h_t_vs_EtaPiPlusMass->ProjectionX("h_Mass_LowT",1,10);
    TH1D* h_Mass_HighT=h_t_vs_EtaPiPlusMass->ProjectionX("h_Mass_HighT",11,20);

    cout<<"Entries with 0 < -t < 0.5:   "<<h_Mass_LowT->Integral()<<endl;
    cout<<"Entries with 0.5 < -t < 1.0: "<<h_Mass_HighT->Integral()<<endl;

    // Hide the statistics box (it would cover the plots).
    gStyle->SetOptStat(0);

    // ---- 5. Draw and save ----
    TCanvas *c1=new TCanvas("c1","",1500,500);
    c1->Divide(3,1);

    c1->cd(1);
    h_PiPlusNeutronMass_vs_EtaPiPlusMass->SetMinimum(0.1);
    h_PiPlusNeutronMass_vs_EtaPiPlusMass->Draw("colz");   // "colz" = colours with a colour scale

    c1->cd(2);
    h_t_vs_EtaPiPlusMass->SetMinimum(0.1);
    h_t_vs_EtaPiPlusMass->Draw("colz");

    c1->cd(3);
    h_Mass_LowT->SetLineColor(kBlue);
    h_Mass_LowT->SetLineWidth(2);
    h_Mass_HighT->SetLineColor(kRed);
    h_Mass_HighT->SetLineWidth(2);
    h_Mass_LowT->SetTitle(";M(#eta#pi^{+}) (GeV);Entries");
    h_Mass_LowT->SetMinimum(0);
    h_Mass_LowT->SetMaximum(1500);
    h_Mass_LowT->Draw("hist");
    h_Mass_HighT->Draw("hist same");

    TLegend *leg=new TLegend(0.45,0.75,0.88,0.88);
    leg->AddEntry(h_Mass_LowT,"0 < -t < 0.5 GeV^{2}","l");
    leg->AddEntry(h_Mass_HighT,"0.5 < -t < 1.0 GeV^{2}","l");
    leg->Draw();

    c1->SaveAs("ex5_2d.pdf");
    c1->SaveAs("ex5_2d.png");
}
