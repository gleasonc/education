// =====================================================================
// Beginner exercise 8: Angles in the Gottfried-Jackson (GJ) frame.
//
// WHY: the spin of a meson shows up in the ANGLE at which it decays.
// A spin-0 meson like the a0(980) decays equally in all directions;
// a spin-2 meson like the a2(1320) does not. To measure that angle we
// must look at the decay from the meson's own point of view: its
// "rest frame", where it is not moving.
//
// WHAT A BOOST IS: the same event looks different to observers moving
// at different speeds. A "boost" re-writes all the four-vectors as seen
// by an observer moving with some velocity. TLorentzVector does this
// for us with .Boost(velocity).
//
// THE STEPS (the same ones as in MakeFlat.C):
//   1. Boost to the centre-of-mass (CM) frame of beam + target, where
//      the total momentum is zero.
//   2. Boost from the CM frame to the eta pi+ rest frame.
//   3. Choose the axes:
//        z = the beam direction (seen in the eta pi+ rest frame)
//        y = perpendicular to the plane made by the beam and the eta pi+
//        x = perpendicular to both y and z
//   4. Measure the eta direction:
//        theta_GJ = angle between the eta and the z axis (we plot cos theta)
//        phi_GJ   = angle around the z axis, starting from the x axis
//
// HOW TO RUN:
//     cd etapiplus_analysis_intro/cpp_beginner
//     root -l ex8_gj_angles.C
//
// WHAT YOU SHOULD SEE: four plots saved as ex8_gj_angles.pdf (and .png):
//   cos(theta_GJ), phi_GJ, cos(theta_GJ) vs M(eta pi+), and cos(theta_GJ)
//   for the a0(980) and a2(1320) mass regions. Also printed:
//     a0 window entries: 8314.83
//     a2 window entries: 8536
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
#include "TVector3.h"
#include "TMath.h"

void ex8_gj_angles()
{
    // ---- 1. Open the file and get the tree ----
    TFile *fileIn=new TFile("../../flat_pipneta__B4_v6.root");
    TTree *tree=(TTree*)fileIn->Get("pipneta__B4");

    // ---- 2. Variables for the branches we need ----
    double PiPlus_Px, PiPlus_Py, PiPlus_Pz, PiPlus_E;
    double Eta_Px, Eta_Py, Eta_Pz, Eta_E;
    double Beam_Px, Beam_Py, Beam_Pz, Beam_E;
    double Chi2_NDF, Mandelstam_t, MissingMassSquared, accidweight;

    tree->SetBranchAddress("PiPlus_Px", &PiPlus_Px );
    tree->SetBranchAddress("PiPlus_Py", &PiPlus_Py );
    tree->SetBranchAddress("PiPlus_Pz", &PiPlus_Pz );
    tree->SetBranchAddress("PiPlus_E", &PiPlus_E );

    tree->SetBranchAddress("Eta_Px", &Eta_Px );
    tree->SetBranchAddress("Eta_Py", &Eta_Py );
    tree->SetBranchAddress("Eta_Pz", &Eta_Pz );
    tree->SetBranchAddress("Eta_E", &Eta_E );

    tree->SetBranchAddress("Beam_Px", &Beam_Px );
    tree->SetBranchAddress("Beam_Py", &Beam_Py );
    tree->SetBranchAddress("Beam_Pz", &Beam_Pz );
    tree->SetBranchAddress("Beam_E", &Beam_E );

    tree->SetBranchAddress("Chi2_NDF", &Chi2_NDF );
    tree->SetBranchAddress("Mandelstam_t", &Mandelstam_t );
    tree->SetBranchAddress("MissingMassSquared", &MissingMassSquared );
    tree->SetBranchAddress("accidweight", &accidweight );

    // ---- 3. Empty histograms ----
    TH1F* h_CosTheta=new TH1F("h_CosTheta",";cos#theta_{GJ};Entries",50,-1.0,1.0);
    TH1F* h_Phi=new TH1F("h_Phi",";#phi_{GJ} (degrees);Entries",60,-180,180);
    TH2F* h_CosTheta_vs_EtaPiPlusMass=new TH2F("h_CosTheta_vs_EtaPiPlusMass",";M(#eta#pi^{+}) (GeV);cos#theta_{GJ}",45,0.7,2.5,40,-1,1);
    TH1F* h_CosTheta_a0=new TH1F("h_CosTheta_a0",";cos#theta_{GJ};Fraction",25,-1.0,1.0);
    TH1F* h_CosTheta_a2=new TH1F("h_CosTheta_a2",";cos#theta_{GJ};Fraction",25,-1.0,1.0);

    // The target proton sits still in the lab: momentum 0, energy = its mass
    TLorentzVector TargetP4;
    TargetP4.SetPxPyPzE(0,0,0,0.938272);

    // ---- 4. Loop over every entry ----
    Long64_t NEntries=tree->GetEntries();
    cout<<"There are "<<NEntries<<" entries."<<endl;

    for (Long64_t i_entry=0; i_entry<NEntries; i_entry++) {
        tree->GetEntry(i_entry);
        if (i_entry%1000000==0)
            cout<<i_entry<<" events processed"<<endl;

        // Standard cuts
        if (Chi2_NDF>10)
            continue;
        if (-Mandelstam_t>1)
            continue;
        if (MissingMassSquared>0.05)
            continue;
        if (MissingMassSquared<-0.05)
            continue;

        // Four-vectors in the lab
        TLorentzVector BeamP4;
        BeamP4.SetPxPyPzE(Beam_Px, Beam_Py, Beam_Pz, Beam_E);
        TLorentzVector pipP4;
        pipP4.SetPxPyPzE(PiPlus_Px, PiPlus_Py, PiPlus_Pz, PiPlus_E);
        TLorentzVector etaP4;
        etaP4.SetPxPyPzE(Eta_Px, Eta_Py, Eta_Pz, Eta_E);

        TLorentzVector EtaPiPlusP4=etaP4+pipP4;
        double EtaPiPlusMass=EtaPiPlusP4.M();

        // ---- Step 1: boost to the CM frame ----
        // The velocity of the CM frame is (total momentum)/(total energy).
        // Boosting by MINUS that velocity takes us into the CM frame.
        TLorentzVector locCoMP4=BeamP4+TargetP4;
        TVector3 boostCoM=-(locCoMP4.Vect())*(1.0/locCoMP4.E());

        // Make copies, then boost the copies (the lab versions stay as they were).
        TLorentzVector locBeamP4_CM=BeamP4;
        TLorentzVector locEtaP4_CM=etaP4;
        TLorentzVector locEtaPiPlusP4_CM=EtaPiPlusP4;
        locBeamP4_CM.Boost(boostCoM);
        locEtaP4_CM.Boost(boostCoM);
        locEtaPiPlusP4_CM.Boost(boostCoM);

        // ---- Step 2: boost from the CM frame to the eta pi+ rest frame ----
        TVector3 boostGJ=-(locEtaPiPlusP4_CM.Vect())*(1.0/locEtaPiPlusP4_CM.E());

        TLorentzVector locBeamP4_GJ=locBeamP4_CM;
        TLorentzVector locEtaP4_GJ=locEtaP4_CM;
        locBeamP4_GJ.Boost(boostGJ);
        locEtaP4_GJ.Boost(boostGJ);

        // ---- Step 3: the GJ axes ----
        // .Unit() makes a vector of length 1 pointing the same way.
        // a.Cross(b) is a vector perpendicular to both a and b.
        TVector3 z_hat_GJ=locBeamP4_GJ.Vect().Unit();
        TVector3 y_hat_GJ=locBeamP4_CM.Vect().Cross(locEtaPiPlusP4_CM.Vect()).Unit();
        TVector3 x_hat_GJ=y_hat_GJ.Cross(z_hat_GJ);

        // ---- Step 4: the eta direction in those axes ----
        // a.Dot(b) with b of length 1 = how far a points along b.
        TVector3 etaVec=locEtaP4_GJ.Vect();
        TVector3 v(etaVec.Dot(x_hat_GJ), etaVec.Dot(y_hat_GJ), etaVec.Dot(z_hat_GJ));

        double cosTheta=v.CosTheta();
        double phi=v.Phi()*180./TMath::Pi();     // Phi() is in radians; change to degrees

        // ---- Fill ----
        h_CosTheta->Fill(cosTheta, accidweight);
        h_Phi->Fill(phi, accidweight);
        h_CosTheta_vs_EtaPiPlusMass->Fill(EtaPiPlusMass, cosTheta, accidweight);

        if (EtaPiPlusMass>0.92 && EtaPiPlusMass<1.04)      // a0(980) region
            h_CosTheta_a0->Fill(cosTheta, accidweight);
        if (EtaPiPlusMass>1.24 && EtaPiPlusMass<1.40)      // a2(1320) region
            h_CosTheta_a2->Fill(cosTheta, accidweight);
    }

    cout<<"a0 window entries: "<<h_CosTheta_a0->Integral()<<endl;
    cout<<"a2 window entries: "<<h_CosTheta_a2->Integral()<<endl;

    // Scale the two mass-window plots to area 1 so their shapes can be compared.
    h_CosTheta_a0->Scale(1.0/h_CosTheta_a0->Integral());
    h_CosTheta_a2->Scale(1.0/h_CosTheta_a2->Integral());

    // ---- 5. Draw and save ----
    gStyle->SetOptStat(0);

    TCanvas *c1=new TCanvas("c1","",1200,1000);
    c1->Divide(2,2);
    c1->cd(1);
    h_CosTheta->SetMinimum(0);
    h_CosTheta->Draw("hist");
    c1->cd(2);
    h_Phi->SetMinimum(0);
    h_Phi->Draw("hist");
    c1->cd(3);
    h_CosTheta_vs_EtaPiPlusMass->SetMinimum(0.1);
    h_CosTheta_vs_EtaPiPlusMass->Draw("colz");
    c1->cd(4);
    h_CosTheta_a0->SetLineColor(kBlue);
    h_CosTheta_a0->SetLineWidth(2);
    h_CosTheta_a0->SetMinimum(0);
    h_CosTheta_a0->SetMaximum(0.08);
    h_CosTheta_a0->Draw("hist");
    h_CosTheta_a2->SetLineColor(kRed);
    h_CosTheta_a2->SetLineWidth(2);
    h_CosTheta_a2->Draw("hist same");

    TLegend *leg=new TLegend(0.15,0.75,0.6,0.88);
    leg->AddEntry(h_CosTheta_a0,"a_{0}(980): 0.92-1.04 GeV","l");
    leg->AddEntry(h_CosTheta_a2,"a_{2}(1320): 1.24-1.40 GeV","l");
    leg->Draw();

    c1->SaveAs("ex8_gj_angles.pdf");
    c1->SaveAs("ex8_gj_angles.png");
}
