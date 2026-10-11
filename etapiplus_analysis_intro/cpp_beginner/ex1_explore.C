// =====================================================================
// Beginner exercise 1: What is inside the file?
//
// HOW TO RUN (in a terminal):
//     cd etapiplus_analysis_intro/cpp_beginner
//     root -l ex1_explore.C
// Type .q and press Enter to quit ROOT when it finishes.
//
// WHAT YOU SHOULD SEE:
//     A list of the 23 branches in the tree, the number of entries
//     (6091343), the values of a few branches for the first 5 entries,
//     and a plot of the beam energy saved as ex1_BeamE.pdf (and .png).
// =====================================================================

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

// ROOT runs the function with the same name as the file.
void ex1_explore()
{
    // Open the data file. The "../../" means "go up two folders".
    // If your file is somewhere else, change the text in quotes.
    TFile *fileIn=new TFile("../../flat_pipneta__B4_v6.root");

    // Inside the file is a TTree (a table) called "pipneta__B4".
    // Each row is one event. Each column is called a "branch".
    TTree *tree=(TTree*)fileIn->Get("pipneta__B4");

    // Print the list of branches and their types.
    // /D means a decimal number (double), /I means a whole number (int).
    tree->Print();

    // How many rows (entries) are in the tree?
    Long64_t NEntries=tree->GetEntries();
    cout<<"There are "<<NEntries<<" entries."<<endl;

    // ---------------------------------------------------------------
    // Read a few values ourselves.
    // Step 1: make one variable for each branch we want to read.
    //         The name does not have to match the branch, but it is
    //         easier to remember if it does.
    // ---------------------------------------------------------------
    double Beam_E, Chi2_NDF, Mandelstam_t, accidweight;
    int pol;

    // Step 2: tell the tree where to put each branch's value.
    //         The & means "put the value into this variable".
    tree->SetBranchAddress("Beam_E", &Beam_E );
    tree->SetBranchAddress("Chi2_NDF", &Chi2_NDF );
    tree->SetBranchAddress("Mandelstam_t", &Mandelstam_t );
    tree->SetBranchAddress("accidweight", &accidweight );
    tree->SetBranchAddress("pol", &pol );

    // Step 3: loop over the first 5 entries and print the values.
    // This "for" line means: start i_entry at 0, keep going while it is
    // less than 5, and add 1 to it each time around.
    for (Long64_t i_entry=0; i_entry<5; i_entry++) {

        // Read entry number i_entry. This fills all the variables above.
        tree->GetEntry(i_entry);

        cout<<"Entry "<<i_entry
            <<"   Beam_E = "<<Beam_E
            <<"   Chi2_NDF = "<<Chi2_NDF
            <<"   Mandelstam_t = "<<Mandelstam_t
            <<"   accidweight = "<<accidweight
            <<"   pol = "<<pol<<endl;
    }

    // ---------------------------------------------------------------
    // The quickest way to look at one branch: tree->Draw("branch name").
    // ROOT makes a histogram for you and draws it on a canvas.
    // ---------------------------------------------------------------
    TCanvas *c1=new TCanvas("c1","",800,600);
    tree->Draw("Beam_E");
    c1->SaveAs("ex1_BeamE.pdf");     // pdf for printing
    c1->SaveAs("ex1_BeamE.png");     // png for a quick look
}
