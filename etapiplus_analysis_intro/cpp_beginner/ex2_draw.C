// =====================================================================
// Beginner exercise 2: Fill and draw histograms.
//
// Every macro from here on has the same 5 parts:
//     1. open the file and get the tree
//     2. make a variable for each branch and connect it (SetBranchAddress)
//     3. make empty histograms
//     4. loop over every entry and fill the histograms
//     5. draw the histograms and save a picture
//
// HOW TO RUN:
//     cd tutorial/cpp_beginner
//     root -l ex2_draw.C
//
// WHAT YOU SHOULD SEE:
//     "There are 6091343 entries." and progress messages, then a
//     window with two plots, also saved as ex2_draw.pdf (and .png) in this folder.
//     Left: the beam energy. Right: the kinematic-fit chi2/ndf.
// =====================================================================

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

void ex2_draw()
{
    // ---- 1. Open the file and get the tree ----
    TFile *fileIn=new TFile("../../flat_pipneta__B4_v6.root");
    TTree *tree=(TTree*)fileIn->Get("pipneta__B4");

    // ---- 2. Variables for the branches we need ----
    double Beam_E, Chi2_NDF;

    tree->SetBranchAddress("Beam_E", &Beam_E );
    tree->SetBranchAddress("Chi2_NDF", &Chi2_NDF );

    // ---- 3. Empty histograms ----
    // new TH1F("name", "title;x-axis label;y-axis label", number of bins, lowest x, highest x)
    TH1F* h_BeamE=new TH1F("h_BeamE","Beam energy;E_{#gamma} (GeV);Entries",200,8.0,12.0);
    TH1F* h_ChiSqNDF=new TH1F("h_ChiSqNDF","Kinematic fit;#chi^{2}/NDF;Entries",100,0,50);

    // ---- 4. Loop over every entry ----
    Long64_t NEntries=tree->GetEntries();
    cout<<"There are "<<NEntries<<" entries."<<endl;

    for (Long64_t i_entry=0; i_entry<NEntries; i_entry++) {
        tree->GetEntry(i_entry);

        // Print a message every 1,000,000 entries so we know it is working.
        if (i_entry%1000000==0)
            cout<<i_entry<<" events processed"<<endl;

        // Add this entry's values to the histograms.
        h_BeamE->Fill(Beam_E);
        h_ChiSqNDF->Fill(Chi2_NDF);
    }

    // ---- 5. Draw and save ----
    TCanvas *c1=new TCanvas("c1","",1200,500);
    c1->Divide(2,1);          // split the canvas into 2 pads, side by side

    c1->cd(1);                // go to the left pad
    h_BeamE->Draw();

    c1->cd(2);                // go to the right pad
    gPad->SetLogy();          // log scale on the y-axis of this pad
    h_ChiSqNDF->Draw();

    c1->SaveAs("ex2_draw.pdf");
    c1->SaveAs("ex2_draw.png");
}
