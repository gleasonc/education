// =====================================================================
// Beginner exercise 3: Build four-vectors and plot invariant masses.
//
// The tree stores (Px, Py, Pz, E) for the eta, the pi+ and the neutron.
// Adding the four-vectors of two particles and taking the mass tells us
// whether they came from the decay of one heavier particle: if so, we
// see a peak at that particle's mass.
//
// We plot three pairs (no cuts yet):
//     M(eta pi+)    peaks: a0(980), a2(1320)
//     M(eta n)      peaks: N* resonances?
//     M(pi+ n)      peaks: Delta(1232)?
//
// HOW TO RUN:
//     cd tutorial/cpp_beginner
//     root -l ex3_masses.C
//
// WHAT YOU SHOULD SEE: three mass plots side by side, saved as
// ex3_masses.pdf (and .png). Without cuts the peaks are weak; exercise 4
// adds cuts.
// =====================================================================

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

void ex3_masses()
{
    // ---- 1. Open the file and get the tree ----
    TFile *fileIn=new TFile("../../flat_pipneta__B4_v6.root");
    TTree *tree=(TTree*)fileIn->Get("pipneta__B4");

    // ---- 2. Variables for the branches we need ----
    double PiPlus_Px, PiPlus_Py, PiPlus_Pz, PiPlus_E;
    double Eta_Px, Eta_Py, Eta_Pz, Eta_E;
    double Neutron_Px, Neutron_Py, Neutron_Pz, Neutron_E;
    double accidweight;

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

    tree->SetBranchAddress("accidweight", &accidweight );

    // ---- 3. Empty histograms ----
    TH1F* h_EtaPiPlusMass=new TH1F("h_EtaPiPlusMass",";M(#eta#pi^{+}) (GeV);Entries",100,0.6,3.6);
    TH1F* h_EtaNeutronMass=new TH1F("h_EtaNeutronMass",";M(#eta n) (GeV);Entries",100,1.4,4.4);
    TH1F* h_PiPlusNeutronMass=new TH1F("h_PiPlusNeutronMass",";M(#pi^{+}n) (GeV);Entries",100,1.0,4.0);

    // ---- 4. Loop over every entry ----
    Long64_t NEntries=tree->GetEntries();
    cout<<"There are "<<NEntries<<" entries."<<endl;

    for (Long64_t i_entry=0; i_entry<NEntries; i_entry++) {
        tree->GetEntry(i_entry);
        if (i_entry%1000000==0)
            cout<<i_entry<<" events processed"<<endl;

        // One four-vector per particle: (Px, Py, Pz, E)
        TLorentzVector pipP4;
        pipP4.SetPxPyPzE(PiPlus_Px, PiPlus_Py, PiPlus_Pz, PiPlus_E);
        TLorentzVector etaP4;
        etaP4.SetPxPyPzE(Eta_Px, Eta_Py, Eta_Pz, Eta_E);
        TLorentzVector neutronP4;
        neutronP4.SetPxPyPzE(Neutron_Px, Neutron_Py, Neutron_Pz, Neutron_E);

        // Add two four-vectors, then .M() gives the invariant mass.
        TLorentzVector EtaPiPlusP4=etaP4+pipP4;
        double EtaPiPlusMass=EtaPiPlusP4.M();

        TLorentzVector EtaNeutronP4=etaP4+neutronP4;
        double EtaNeutronMass=EtaNeutronP4.M();

        TLorentzVector PiPlusNeutronP4=pipP4+neutronP4;
        double PiPlusNeutronMass=PiPlusNeutronP4.M();

        // Fill with accidweight: +1 for in-time beam photons, negative for
        // out-of-time ones, so accidental photons are subtracted.
        h_EtaPiPlusMass->Fill(EtaPiPlusMass, accidweight);
        h_EtaNeutronMass->Fill(EtaNeutronMass, accidweight);
        h_PiPlusNeutronMass->Fill(PiPlusNeutronMass, accidweight);
    }

    // ---- 5. Draw and save ----
    TCanvas *c1=new TCanvas("c1","",1500,500);
    c1->Divide(3,1);
    c1->cd(1);
    h_EtaPiPlusMass->Draw("hist");
    c1->cd(2);
    h_EtaNeutronMass->Draw("hist");
    c1->cd(3);
    h_PiPlusNeutronMass->Draw("hist");
    c1->SaveAs("ex3_masses.pdf");
    c1->SaveAs("ex3_masses.png");
}
