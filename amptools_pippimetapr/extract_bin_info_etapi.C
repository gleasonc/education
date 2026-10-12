// extract_bin_info_etapi.C
// ROOT macro: summarize ROOTDataReader-format trees (one csv row per file)
//
// Columns: file, entries, weighted_events, mass_mean, mass_rms, mass_min, mass_max
// The mass is computed from the final-state particles listed in 'indices'
// (default "1,2" = EtaPrime, Pi- for "reaction EtaPi Beam Proton EtaPrime Pi- Pi+",
//  where the proton is final-state index 0).
//
// Usage: root -l -b -q 'extract_bin_info_etapi.C("list.txt","data.csv","kin","1,2")'

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

void extract_bin_info_etapi( const char* listName, const char* outName,
                             const char* treeName = "kin", const char* indices = "1,2" )
{
    std::vector<int> idx;
    { std::istringstream ss( indices ); std::string t;
      while( std::getline( ss, t, ',' ) ) if( !t.empty() ) idx.push_back( std::stoi( t ) ); }

    std::vector<std::string> files;
    { std::ifstream in( listName ); std::string l;
      while( std::getline( in, l ) ) if( !l.empty() ) files.push_back( l ); }

    std::ofstream out( outName );
    out << "file,entries,weighted_events,mass_mean,mass_rms,mass_min,mass_max\n";
    out.precision( 10 );

    const int kMax = 32;
    for( auto& fname : files ){
        TFile f( fname.c_str() );
        TTree* t = f.IsOpen() ? (TTree*)f.Get( treeName ) : nullptr;
        if( !t ){ std::cerr << "Warning: no tree '" << treeName << "' in " << fname << std::endl; continue; }

        Int_t nfs = 0; Float_t E[kMax], Px[kMax], Py[kMax], Pz[kMax], w = 1;
        t->SetBranchAddress( "NumFinalState", &nfs );
        t->SetBranchAddress( "E_FinalState",  E  );
        t->SetBranchAddress( "Px_FinalState", Px );
        t->SetBranchAddress( "Py_FinalState", Py );
        t->SetBranchAddress( "Pz_FinalState", Pz );
        bool hasW = t->GetBranch( "Weight" ) != nullptr;
        if( hasW ) t->SetBranchAddress( "Weight", &w );

        double sw = 0, swm = 0, swm2 = 0, mmin = 1e30, mmax = -1e30;
        Long64_t n = t->GetEntries();
        for( Long64_t i = 0; i < n; ++i ){
            t->GetEntry( i );
            TLorentzVector p;
            for( int k : idx ) if( k < nfs ) p += TLorentzVector( Px[k], Py[k], Pz[k], E[k] );
            double m = p.M();
            sw += w; swm += w * m; swm2 += w * m * m;
            mmin = std::min( mmin, m ); mmax = std::max( mmax, m );
        }
        double mean = sw != 0 ? swm / sw : 0;
        double rms  = sw != 0 ? std::sqrt( std::max( 0., swm2 / sw - mean * mean ) ) : 0;
        out << fname << "," << n << "," << sw << "," << mean << "," << rms << ","
            << ( n ? mmin : 0 ) << "," << ( n ? mmax : 0 ) << "\n";
        std::cout << "Processed " << fname << " (" << n << " entries)" << std::endl;
    }
    std::cout << "Wrote " << outName << std::endl;
}
