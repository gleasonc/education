// =====================================================================
// All-in-one macro (optional): exercises 2 to 5 in one file.
//
// This is laid out like a typical analysis macro you will meet in the
// group (for example MakeFlat.C): read the tree once, fill many
// histograms in one loop, then draw them all at the end.
//
// HOW TO RUN:
//     cd tutorial/cpp_beginner
//     root -l -b -q MakePlots_EtaPiPlus.C
//
// WHAT YOU SHOULD SEE: progress messages, then five PDF files in this
// folder: ChiSqNDF.pdf, NeutronBeta.pdf, Masses.pdf,
// EtaPiPlusMass_Cuts.pdf and Masses2D.pdf.
// =====================================================================

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLorentzVector.h"

void MakePlots_EtaPiPlus()
{
    TFile *fileIn=new TFile("../../flat_pipneta__B4_v6.root");
    TTree *tree=(TTree*)fileIn->Get("pipneta__B4");

    // One variable per branch
    double MissingMassSquared, accidweight, PiPlus_Px, PiPlus_Py, PiPlus_Pz, PiPlus_E, Eta_Px, Eta_Py, Eta_Pz, Eta_E, Neutron_Px, Neutron_Py, Neutron_Pz, Neutron_E, Beam_Px, Beam_Py, Beam_Pz, Beam_E, Neutron_Beta, Mandelstam_t, Chi2_NDF;
    int pol;

    tree->SetBranchAddress("accidweight", &accidweight );
    tree->SetBranchAddress("Mandelstam_t", &Mandelstam_t );
    tree->SetBranchAddress("MissingMassSquared", &MissingMassSquared );

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

    tree->SetBranchAddress("Neutron_Px", &Neutron_Px );
    tree->SetBranchAddress("Neutron_Py", &Neutron_Py );
    tree->SetBranchAddress("Neutron_Pz", &Neutron_Pz );
    tree->SetBranchAddress("Neutron_E", &Neutron_E );
    tree->SetBranchAddress("Neutron_Beta", &Neutron_Beta );
    tree->SetBranchAddress("Chi2_NDF", &Chi2_NDF );

    tree->SetBranchAddress("pol", &pol );

    // Histograms filled before any cuts
    TH1F* h_ChiSqNDF=new TH1F("h_ChiSqNDF",";#chi^{2}/NDF",100,0,15);
    TH1F* h_MissingMassSquared_NoCuts=new TH1F("h_MissingMassSquared_NoCuts",";MM^{2} (GeV^{2})",100,-0.1,0.1);
    TH2F* h_ChiSqNDF_vs_MissingMassSq=new TH2F("h_ChiSqNDF_vs_MissingMassSq",";MM^{2} (GeV^{2});#chi^{2}/NDF",100,-0.1,0.1,100,0,20);
    TH1F* h_NeutronBeta=new TH1F("h_NeutronBeta",";#beta_{n}",100,0,1.1);
    TH1F* h_EtaPiPlusMass_NoCuts=new TH1F("h_EtaPiPlusMass_NoCuts",";M(#eta#pi^{+}) (GeV)",90,0.7,2.5);

    // Histograms filled after the cuts
    TH1F* h_Mandelstam_t=new TH1F("h_Mandelstam_t",";-t (GeV^{2})",100,0,2);
    TH1F* h_EtaPiPlusMass=new TH1F("h_EtaPiPlusMass",";M(#eta#pi^{+}) (GeV)",90,0.7,2.5);
    TH1F* h_EtaNeutronMass=new TH1F("h_EtaNeutronMass",";M(#eta n) (GeV)",100,1.4,4.4);
    TH1F* h_PiPlusNeutronMass=new TH1F("h_PiPlusNeutronMass",";M(#pi^{+}n) (GeV)",100,1.0,4.0);
    TH2F* h_EtaPiPlusMass_vs_PiPlusNeutronMass=new TH2F("h_EtaPiPlusMass_vs_PiPlusNeutronMass",";M(#pi^{+}n) (GeV);M(#eta#pi^{+}) (GeV)",100,1.0,4.0,90,0.7,2.5);
    TH2F* h_EtaPiPlusMass_vs_EtaNeutronMass=new TH2F("h_EtaPiPlusMass_vs_EtaNeutronMass",";M(#eta n) (GeV);M(#eta#pi^{+}) (GeV)",100,1.4,4.4,90,0.7,2.5);

    // Same mass, with one extra cut each, to see what the cut does
    TH1F* h_EtaPiPlusMass_CutPiPlusN=new TH1F("h_EtaPiPlusMass_CutPiPlusN",";M(#eta#pi^{+}) (GeV)",90,0.7,2.5);
    TH1F* h_EtaPiPlusMass_CutEtaN=new TH1F("h_EtaPiPlusMass_CutEtaN",";M(#eta#pi^{+}) (GeV)",90,0.7,2.5);

    Long64_t NEntries=tree->GetEntries();
    cout<<"There are "<<NEntries<<" entries."<<endl;

    for (Long64_t i_entry=0; i_entry<NEntries; i_entry++) {
        tree->GetEntry(i_entry);
        if (i_entry%1000000==0)
            cout<<i_entry<<" events processed"<<endl;

        // Four-vectors
        TLorentzVector neutronP4;
        neutronP4.SetPxPyPzE(Neutron_Px, Neutron_Py, Neutron_Pz, Neutron_E);
        TLorentzVector etaP4;
        etaP4.SetPxPyPzE(Eta_Px, Eta_Py, Eta_Pz, Eta_E);
        TLorentzVector pipP4;
        pipP4.SetPxPyPzE(PiPlus_Px, PiPlus_Py, PiPlus_Pz, PiPlus_E);

        // Invariant masses of each pair
        TLorentzVector EtaPiPlusP4=etaP4+pipP4;
        double EtaPiPlusMass=EtaPiPlusP4.M();

        TLorentzVector EtaNeutronP4=etaP4+neutronP4;
        double EtaNeutronMass=EtaNeutronP4.M();

        TLorentzVector PiPlusNeutronP4=pipP4+neutronP4;
        double PiPlusNeutronMass=PiPlusNeutronP4.M();

        // Before cuts
        h_ChiSqNDF->Fill(Chi2_NDF, accidweight);
        h_ChiSqNDF_vs_MissingMassSq->Fill(MissingMassSquared,Chi2_NDF,accidweight);
        h_MissingMassSquared_NoCuts->Fill(MissingMassSquared, accidweight);
        h_NeutronBeta->Fill(Neutron_Beta, accidweight);
        h_EtaPiPlusMass_NoCuts->Fill(EtaPiPlusMass, accidweight);

        // Cuts
        if (Chi2_NDF>10)
            continue;
        if (Neutron_Beta>0.9)          // already applied when the tree was made, so this removes nothing
            continue;
        if (MissingMassSquared>0.05 || MissingMassSquared<-0.05)    // || means "or"
            continue;

        h_Mandelstam_t->Fill(-Mandelstam_t, accidweight);

        if (-Mandelstam_t>1)
            continue;

        // After cuts
        h_EtaPiPlusMass->Fill(EtaPiPlusMass, accidweight);
        h_EtaNeutronMass->Fill(EtaNeutronMass, accidweight);
        h_PiPlusNeutronMass->Fill(PiPlusNeutronMass, accidweight);
        h_EtaPiPlusMass_vs_PiPlusNeutronMass->Fill(PiPlusNeutronMass, EtaPiPlusMass, accidweight);
        h_EtaPiPlusMass_vs_EtaNeutronMass->Fill(EtaNeutronMass, EtaPiPlusMass, accidweight);

        // Extra cuts: remove low-mass pi+ n (Delta region) or low-mass eta n (N* region)
        if (PiPlusNeutronMass>1.4)
            h_EtaPiPlusMass_CutPiPlusN->Fill(EtaPiPlusMass, accidweight);
        if (EtaNeutronMass>2.0)
            h_EtaPiPlusMass_CutEtaN->Fill(EtaPiPlusMass, accidweight);
    }

    // Drawing
    TCanvas *c1=new TCanvas("c1","",800,600);
    c1->Divide(2,2);
    c1->cd(1);
    h_ChiSqNDF->Draw("hist");
    c1->cd(2);
    h_ChiSqNDF_vs_MissingMassSq->SetMinimum(0.1);
    h_ChiSqNDF_vs_MissingMassSq->Draw("colz");
    c1->cd(3);
    h_MissingMassSquared_NoCuts->Draw("hist");
    c1->cd(4);
    h_Mandelstam_t->Draw("hist");
    c1->SaveAs("ChiSqNDF.pdf");

    TCanvas *c2=new TCanvas("c2","",800,600);
    h_NeutronBeta->Draw("hist");
    c2->SaveAs("NeutronBeta.pdf");

    TCanvas *c3=new TCanvas("c3","",1200,400);
    c3->Divide(3,1);
    c3->cd(1);
    h_EtaPiPlusMass->Draw("hist");
    c3->cd(2);
    h_EtaNeutronMass->Draw("hist");
    c3->cd(3);
    h_PiPlusNeutronMass->Draw("hist");
    c3->SaveAs("Masses.pdf");

    TCanvas *c4=new TCanvas("c4","",800,600);
    h_EtaPiPlusMass_NoCuts->SetMarkerStyle(20);
    h_EtaPiPlusMass_NoCuts->SetMarkerColor(kBlack);
    h_EtaPiPlusMass_NoCuts->SetLineColor(kBlack);
    h_EtaPiPlusMass_NoCuts->SetMinimum(0);
    h_EtaPiPlusMass_NoCuts->SetMaximum(14000);    // leave room for the legend
    h_EtaPiPlusMass_NoCuts->Draw();

    h_EtaPiPlusMass->SetMarkerStyle(20);
    h_EtaPiPlusMass->SetMarkerColor(kRed);
    h_EtaPiPlusMass->SetLineColor(kRed);
    h_EtaPiPlusMass->Draw("same");

    h_EtaPiPlusMass_CutPiPlusN->SetMarkerStyle(20);
    h_EtaPiPlusMass_CutPiPlusN->SetMarkerColor(kBlue);
    h_EtaPiPlusMass_CutPiPlusN->SetLineColor(kBlue);
    h_EtaPiPlusMass_CutPiPlusN->Draw("same");

    h_EtaPiPlusMass_CutEtaN->SetMarkerStyle(20);
    h_EtaPiPlusMass_CutEtaN->SetMarkerColor(kGreen+2);
    h_EtaPiPlusMass_CutEtaN->SetLineColor(kGreen+2);
    h_EtaPiPlusMass_CutEtaN->Draw("same");

    TLegend *leg=new TLegend(0.15,0.68,0.5,0.88);
    leg->AddEntry(h_EtaPiPlusMass_NoCuts,"No cuts","p");
    leg->AddEntry(h_EtaPiPlusMass,"Cuts","p");
    leg->AddEntry(h_EtaPiPlusMass_CutPiPlusN,"Cuts + M(#pi^{+}n) > 1.4","p");
    leg->AddEntry(h_EtaPiPlusMass_CutEtaN,"Cuts + M(#eta n) > 2.0","p");
    leg->Draw();
    c4->SaveAs("EtaPiPlusMass_Cuts.pdf");

    TCanvas *c5=new TCanvas("c5","",1200,500);
    c5->Divide(2,1);
    c5->cd(1);
    h_EtaPiPlusMass_vs_PiPlusNeutronMass->SetMinimum(0.1);
    h_EtaPiPlusMass_vs_PiPlusNeutronMass->Draw("colz");
    c5->cd(2);
    h_EtaPiPlusMass_vs_EtaNeutronMass->SetMinimum(0.1);
    h_EtaPiPlusMass_vs_EtaNeutronMass->Draw("colz");
    c5->SaveAs("Masses2D.pdf");

    cout<<"Entries (weighted) after cuts: "<<h_EtaPiPlusMass->Integral()<<" in 0.7-2.5 GeV"<<endl;
}
