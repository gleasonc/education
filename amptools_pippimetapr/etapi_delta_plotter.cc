#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "TROOT.h"
#include "TH1.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"

#include "AMPTOOLS_DATAIO/EtaPiDeltaPlotGenerator.h"
#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_AMPS/TwoPiAngles.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/Zlm.h"

typedef EtaPiDeltaPlotGenerator PlotGen;

void atiSetup(){

  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface::registerAmplitude( Zlm() );
  AmpToolsInterface::registerDataReader( ROOTDataReader() );
}

//  THE USER SHOULD NOT HAVE TO CHANGE ANYTHING BELOW THIS LINE
// *************************************************************

using namespace std;

void usage( const char* exe ){
    cout << "Usage:" << endl << endl;
    cout << "  " << exe << " <results file name> [options]" << endl << endl;
    cout << "\t -o <file>\t output ROOT file (default etapi_plot.root); the PDF gets the same name with .pdf" << endl;
    cout << "\t -p <file>\t output PDF file (overrides the default above)" << endl;
    cout << "\t -t <title>\t title printed on the PDF (e.g. the mass / t bin)" << endl;
    cout << "\t -s \t\t total only: skip the histograms for each individual sum (much faster)" << endl;
    cout << "\t -G \t\t skip generated MC histograms (faster)" << endl << endl;
}

int main( int argc, char* argv[] ){

    gROOT->SetBatch( kTRUE );   // never open windows

    // ************************
    // usage
    // ************************

    cout << endl << " *** Writing fit-quality histograms (ROOT + PDF) *** " << endl << endl;

    if (argc < 2){
        usage( argv[0] );
        return 0;
    }

    string outName = "etapi_plot.root";
    string pdfName = "";
    string title = "";
    bool perSum = true;     // -s: total only
    bool doGen = true;      // -G: skip generated MC
    string resultsName(argv[1]);
    for (int i = 2; i < argc; i++){

        string arg(argv[i]);

        if (arg == "-o" && i + 1 < argc) outName = argv[++i];
        else if (arg == "-p" && i + 1 < argc) pdfName = argv[++i];
        else if (arg == "-t" && i + 1 < argc) title = argv[++i];
        else if (arg == "-s") perSum = false;
        else if (arg == "-G") doGen = false;
        else {
            usage( argv[0] );
            exit(1);
        }
    }
    if (pdfName.empty()){
        pdfName = outName;
        size_t dot = pdfName.rfind(".root");
        if (dot != string::npos) pdfName.erase(dot);
        pdfName += ".pdf";
    }


    // ************************
    // parse the command line parameters
    // ************************
    cout << "Fit results file name    = " << resultsName << endl;
    cout << "Output ROOT file name    = " << outName << endl;
    cout << "Output PDF file name     = " << pdfName << endl << endl;
    FitResults results( resultsName );
    if( !results.valid() ){

        cout << "Invalid fit results in file:  " << resultsName << endl;
        exit( 1 );
    }

    // ************************
    // set up the plot generator
    // ************************

    atiSetup();
    PlotGen plotGen( results );

    // ************************
    // set up an output ROOT file to store histograms
    // ************************


    TFile* plotfile = new TFile( outName.c_str(), "recreate");
    TH1::AddDirectory(kFALSE);

    string reactionName = results.reactionList()[0];
    plotGen.enableReaction( reactionName );
    vector<string> sums = plotGen.uniqueSums();

    // histograms of the full fit (all sums on), kept for drawing
    map<string, TH1*> totalHists;

    // loop over sum configurations (one for each of the individual contributions, and the combined sum of all)
    for (unsigned int isum = 0; isum <= sums.size(); isum++){

        if (!perSum && isum < sums.size()) continue;

        // turn on all sums by default
        for (unsigned int i = 0; i < sums.size(); i++){
            plotGen.enableSum(i);
        }

        // for individual contributions turn off all sums but the one of interest
        if (isum < sums.size()){
            for (unsigned int i = 0; i < sums.size(); i++){
                if (i != isum) plotGen.disableSum(i);
            }
        }


        // loop over data, accMC, and genMC and kBkgnd
        for (unsigned int iplot = 0; iplot < PlotGenerator::kNumTypes; iplot++){
            if (isum < sums.size() && iplot == PlotGenerator::kData) continue; // only plot data once
            if (!doGen && iplot == PlotGenerator::kGenMC) continue;

            // loop over different variables
            for (unsigned int ivar  = 0; ivar  < EtaPiDeltaPlotGenerator::kNumHists; ivar++){

                // set unique histogram name for each plot (could put in directories...)
                string histname =  "";
                if (ivar == EtaPiDeltaPlotGenerator::kEtaPiMass)  histname += "EtaPiMass";
                else if (ivar == EtaPiDeltaPlotGenerator::kDeltaPPMass)  histname += "DeltaPPMass";
                else if (ivar == EtaPiDeltaPlotGenerator::kEtaCosTheta)  histname += "cosTheta";
                else if (ivar == EtaPiDeltaPlotGenerator::kPhi)  histname += "Phi";
                else if (ivar == EtaPiDeltaPlotGenerator::kt)  histname += "t";
                else continue;

                if (iplot == PlotGenerator::kData) histname += "dat";
                if (iplot == PlotGenerator::kAccMC) histname += "acc";
                if (iplot == PlotGenerator::kGenMC) histname += "gen";
                //if (iplot == PlotGenerator::kBkgnd) histname += "bkgnd";

                if (isum < sums.size()){
                    // get name of sum for naming histogram
                    string sumName = sums[isum];
                    histname += "_";
                    histname += sumName;
                }

                Histogram* hist = plotGen.projection(ivar, reactionName, iplot);
                TH1* thist = hist->toRoot();
                thist->SetName(histname.c_str());
                plotfile->cd();
                thist->Write();

                if (isum == sums.size()) totalHists[histname] = thist;
                else delete thist;
            }
        }
    }


    // ************************
    // draw data vs. fit (accepted MC, all sums) and save to PDF + ROOT file
    // ************************

    gStyle->SetOptStat(0);

    const char* vars[] = { "EtaPiMass", "DeltaPPMass", "cosTheta", "Phi", "t" };
    const char* xTitles[] = { "M(#eta'#pi^{-}) [GeV]", "M(p#pi^{+}) [GeV]",
                              "cos#theta", "#phi [rad]", "-t [GeV^{2}]" };
    const int nVars = 5;

    TCanvas canvas( "fit_quality", title.c_str(), 1500, 900 );
    canvas.Divide( 3, 2 );
    TLegend legend( 0.1, 0.4, 0.9, 0.85 );
    legend.SetBorderSize( 0 );

    for (int v = 0; v < nVars; v++){

        string datName = string(vars[v]) + "dat";
        string accName = string(vars[v]) + "acc";
        if (!totalHists.count( datName ) || !totalHists.count( accName )) continue;
        TH1* dat = totalHists[datName];
        TH1* acc = totalHists[accName];

        canvas.cd( v + 1 );

        int mb = dat->GetMaximumBin();
        double yMax = max( dat->GetBinContent(mb) + dat->GetBinError(mb), acc->GetMaximum() );

        acc->SetTitle( (string(";") + xTitles[v] + ";events").c_str() );
        acc->SetMinimum( 0 );
        acc->SetMaximum( 1.15 * yMax );
        acc->SetFillColor( kAzure - 9 );
        acc->SetLineColor( kAzure - 9 );
        acc->Draw( "hist" );

        dat->SetMarkerStyle( 20 );
        dat->SetMarkerSize( 0.6 );
        dat->SetLineColor( kBlack );
        dat->Draw( "E1 same" );

        if (legend.GetNRows() == 0){
            legend.AddEntry( dat, "data", "pe" );
            legend.AddEntry( acc, "fit (acc. MC)", "f" );
        }
    }

    // last pad: legend and bin / fit information
    canvas.cd( 6 );
    if (!title.empty()) legend.SetHeader( title.c_str() );
    legend.AddEntry( (TObject*)0, ("likelihood = " + to_string( results.likelihood() )).c_str(), "" );
    legend.Draw();

    canvas.SaveAs( pdfName.c_str() );
    plotfile->cd();
    canvas.Write();

    plotfile->Close();
    cout << "Wrote " << outName << " and " << pdfName << endl;

    return 0;

}
