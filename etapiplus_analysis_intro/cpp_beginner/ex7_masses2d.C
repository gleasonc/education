// =====================================================================
// Beginner exercise 7: Two-dimensional mass plots.
//
// Three particles in the final state (eta, pi+, neutron) give three
// pairs, so three two-body masses: M(eta pi+), M(eta n), M(pi+ n).
// Plotting one against another shows where particles were made:
//     - a meson decaying to eta pi+ (a0(980), a2(1320)) has a fixed
//       M(eta pi+), so it makes a HORIZONTAL band (M(eta pi+) is on y);
//     - a baryon decaying to pi+ n (Delta(1232)) or eta n (N*) makes
//       a VERTICAL band.
// A "Dalitz plot" is the same idea using the masses SQUARED.
//
// At the end we cut away the low M(pi+ n) and low M(eta n) regions
// (a "baryon veto", like in MakeFlat.C) and see what it does to M(eta pi+).
//
// HOW TO RUN:
//     cd etapiplus_analysis_intro/cpp_beginner
//     root -l ex7_masses2d.C
//
// WHAT YOU SHOULD SEE:
//     M(eta pi+) entries, standard cuts:  76107.8
//     M(eta pi+) entries, + baryon veto:  40348.1
// and two pictures: ex7_masses2d.pdf and ex7_veto.pdf (and .png).
// =====================================================================

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLorentzVector.h"

void ex7_masses2d()
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

    // ---- 3. Empty histograms ----
    // TH2F: ";x label;y label", x bins, x low, x high, y bins, y low, y high
    TH2F* h_EtaPiPlusMass_vs_EtaNeutronMass=new TH2F("h_EtaPiPlusMass_vs_EtaNeutronMass",";M(#eta n) (GeV);M(#eta#pi^{+}) (GeV)",100,1.4,4.4,90,0.7,2.5);
    TH2F* h_EtaPiPlusMass_vs_PiPlusNeutronMass=new TH2F("h_EtaPiPlusMass_vs_PiPlusNeutronMass",";M(#pi^{+}n) (GeV);M(#eta#pi^{+}) (GeV)",100,1.0,4.0,90,0.7,2.5);
    TH2F* h_EtaNeutronMass_vs_PiPlusNeutronMass=new TH2F("h_EtaNeutronMass_vs_PiPlusNeutronMass",";M(#pi^{+}n) (GeV);M(#eta n) (GeV)",100,1.0,4.0,100,1.4,4.4);
    TH2F* h_Dalitz=new TH2F("h_Dalitz",";M^{2}(#pi^{+}n) (GeV^{2});M^{2}(#eta#pi^{+}) (GeV^{2})",100,1.0,16.0,100,0.0,9.0);

    TH1F* h_EtaPiPlusMass=new TH1F("h_EtaPiPlusMass",";M(#eta#pi^{+}) (GeV);Entries",90,0.7,2.5);
    TH1F* h_EtaPiPlusMass_BaryonVeto=new TH1F("h_EtaPiPlusMass_BaryonVeto",";M(#eta#pi^{+}) (GeV);Entries",90,0.7,2.5);

    // ---- 4. Loop over every entry ----
    Long64_t NEntries=tree->GetEntries();
    cout<<"There are "<<NEntries<<" entries."<<endl;

    for (Long64_t i_entry=0; i_entry<NEntries; i_entry++) {
        tree->GetEntry(i_entry);
        if (i_entry%1000000==0)
            cout<<i_entry<<" events processed"<<endl;

        // Standard cuts (exercise 4)
        if (Chi2_NDF>10)
            continue;
        if (-Mandelstam_t>1)
            continue;
        if (MissingMassSquared>0.05)
            continue;
        if (MissingMassSquared<-0.05)
            continue;

        // Four-vectors
        TLorentzVector pipP4;
        pipP4.SetPxPyPzE(PiPlus_Px, PiPlus_Py, PiPlus_Pz, PiPlus_E);
        TLorentzVector etaP4;
        etaP4.SetPxPyPzE(Eta_Px, Eta_Py, Eta_Pz, Eta_E);
        TLorentzVector neutronP4;
        neutronP4.SetPxPyPzE(Neutron_Px, Neutron_Py, Neutron_Pz, Neutron_E);

        // The three pairs. .M() is the mass, .M2() is the mass squared.
        TLorentzVector EtaPiPlusP4=etaP4+pipP4;
        TLorentzVector EtaNeutronP4=etaP4+neutronP4;
        TLorentzVector PiPlusNeutronP4=pipP4+neutronP4;

        double EtaPiPlusMass=EtaPiPlusP4.M();
        double EtaNeutronMass=EtaNeutronP4.M();
        double PiPlusNeutronMass=PiPlusNeutronP4.M();

        // Fill the 2D plots: Fill(x, y, weight)
        h_EtaPiPlusMass_vs_EtaNeutronMass->Fill(EtaNeutronMass, EtaPiPlusMass, accidweight);
        h_EtaPiPlusMass_vs_PiPlusNeutronMass->Fill(PiPlusNeutronMass, EtaPiPlusMass, accidweight);
        h_EtaNeutronMass_vs_PiPlusNeutronMass->Fill(PiPlusNeutronMass, EtaNeutronMass, accidweight);
        h_Dalitz->Fill(PiPlusNeutronP4.M2(), EtaPiPlusP4.M2(), accidweight);

        h_EtaPiPlusMass->Fill(EtaPiPlusMass, accidweight);

        // Baryon veto: only keep entries with BOTH baryon masses above 2 GeV.
        // && means "and".
        if (PiPlusNeutronMass>2 && EtaNeutronMass>2)
            h_EtaPiPlusMass_BaryonVeto->Fill(EtaPiPlusMass, accidweight);
    }

    cout<<"M(eta pi+) entries, standard cuts:  "<<h_EtaPiPlusMass->Integral()<<endl;
    cout<<"M(eta pi+) entries, + baryon veto:  "<<h_EtaPiPlusMass_BaryonVeto->Integral()<<endl;

    // ---- 5. Draw and save ----
    gStyle->SetOptStat(0);      // hide the statistics box

    TCanvas *c1=new TCanvas("c1","",1200,1000);
    c1->Divide(2,2);
    c1->cd(1);
    h_EtaPiPlusMass_vs_EtaNeutronMass->SetMinimum(0.1);
    h_EtaPiPlusMass_vs_EtaNeutronMass->Draw("colz");
    c1->cd(2);
    h_EtaPiPlusMass_vs_PiPlusNeutronMass->SetMinimum(0.1);
    h_EtaPiPlusMass_vs_PiPlusNeutronMass->Draw("colz");
    c1->cd(3);
    h_EtaNeutronMass_vs_PiPlusNeutronMass->SetMinimum(0.1);
    h_EtaNeutronMass_vs_PiPlusNeutronMass->Draw("colz");
    c1->cd(4);
    h_Dalitz->SetMinimum(0.1);
    h_Dalitz->Draw("colz");
    c1->SaveAs("ex7_masses2d.pdf");
    c1->SaveAs("ex7_masses2d.png");

    TCanvas *c2=new TCanvas("c2","",800,600);
    h_EtaPiPlusMass->SetLineColor(kBlack);
    h_EtaPiPlusMass->SetLineWidth(2);
    h_EtaPiPlusMass->SetMinimum(0);
    h_EtaPiPlusMass->SetMaximum(2500);
    h_EtaPiPlusMass->Draw("hist");
    h_EtaPiPlusMass_BaryonVeto->SetLineColor(kRed);
    h_EtaPiPlusMass_BaryonVeto->SetLineWidth(2);
    h_EtaPiPlusMass_BaryonVeto->Draw("hist same");

    TLegend *leg=new TLegend(0.15,0.75,0.6,0.88);
    leg->AddEntry(h_EtaPiPlusMass,"Standard cuts","l");
    leg->AddEntry(h_EtaPiPlusMass_BaryonVeto,"+ M(#pi^{+}n) > 2 and M(#eta n) > 2 GeV","l");
    leg->Draw();
    c2->SaveAs("ex7_veto.pdf");
    c2->SaveAs("ex7_veto.png");
}
