// split_bins.C
// ROOT macro: split a ROOTDataReader-format tree into mass x (-t) bins in one pass,
// optionally keeping only one beam-polarization orientation.
// SDME version of amptools_pippimetapr/split_bins.C, for
// "reaction Pi0Pim Beam Proton Pi- Pi0" trees from the v2 DSelector.
//
//   mass = M of the final-state particles listed in massIdx
//          (default "1,2" = Pi-, Pi0; the proton is final-state index 0)
//   -t   = -(beam - those particles)^2, i.e. -t to the pi- pi0 system
//          (does not depend on the target mass)
//   pol  = with polBranch set (e.g. "pol"), only events whose Int_t branch
//          polBranch equals polValue are kept (0/45/90/135, -1 = AMO)
//
// Writes <outBase>_m<i>_t<j>.root for every bin (empty bins give an empty tree),
// i = 0..nMass-1, j = 0..nT-1 (j = 0 only when nT <= 0, i.e. no t binning).
// Every branch is copied unchanged (Weight etc. included).  Bins are
// low <= x < high; events outside the ranges are dropped.
//
// Usage:
//   root -l -b -q 'split_bins.C("in.root","out/data",0.6,0.95,1, 0.1,1.1,10)'
//   root -l -b -q 'split_bins.C("in.root","out/data_45",0.6,0.95,1, 0.1,1.1,10, false,"kin","1,2","pol",45)'
//   arguments: in, outBase, massLow, massHigh, nMass, tLow, tHigh, nT,
//              tLog (logarithmic t bins), tree ("kin"), massIdx ("1,2"),
//              polBranch ("" = no selection), polValue

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TBranch.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TSystem.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TString.h"

void split_bins( const char* inName, const char* outBase,
                 double mLow, double mHigh, int nMass,
                 double tLow = 0, double tHigh = 0, int nT = 0, bool tLog = false,
                 const char* treeName = "kin", const char* massIdx = "1,2",
                 const char* polBranch = "", int polValue = 0 )
{
    std::vector<int> idx;
    { std::istringstream ss( massIdx ); std::string s;
      while( std::getline( ss, s, ',' ) ) if( !s.empty() ) idx.push_back( std::stoi( s ) ); }

    bool useT = nT > 0;
    if( !useT ) nT = 1;
    if( useT && tLog && tLow <= 0 ){ std::cerr << "split_bins: log t bins need tLow > 0" << std::endl; return; }

    TFile in( inName );
    TTree* t = in.IsOpen() ? (TTree*)in.Get( treeName ) : nullptr;
    if( !t ){ std::cerr << "split_bins: no tree '" << treeName << "' in " << inName << std::endl; return; }

    // the branches ROOTDataReader needs, with the types it expects
    const char* need[]  = { "NumFinalState", "E_FinalState", "Px_FinalState", "Py_FinalState",
                            "Pz_FinalState", "E_Beam", "Px_Beam", "Py_Beam", "Pz_Beam" };
    const char* types[] = { "Int_t", "Float_t", "Float_t", "Float_t", "Float_t",
                            "Float_t", "Float_t", "Float_t", "Float_t" };
    for( int k = 0; k < 9; ++k ){
        TBranch* br = t->GetBranch( need[k] );
        TLeaf* lf = br ? br->GetLeaf( need[k] ) : nullptr;
        if( !lf ){
            std::cerr << "split_bins: ERROR: " << inName << " has no branch " << need[k]
                      << " in tree '" << treeName << "' (not ROOTDataReader format)" << std::endl;
            gSystem->Exit( 1 );
        }
        if( std::string( lf->GetTypeName() ) != types[k] ){
            std::cerr << "split_bins: ERROR: branch " << need[k] << " in " << inName << " is "
                      << lf->GetTypeName() << ", expected " << types[k] << std::endl;
            gSystem->Exit( 1 );
        }
    }

    bool usePol = polBranch && polBranch[0];
    Int_t pol = 0;
    if( usePol ){
        TBranch* br = t->GetBranch( polBranch );
        TLeaf* lf = br ? br->GetLeaf( polBranch ) : nullptr;
        if( !lf || std::string( lf->GetTypeName() ) != "Int_t" ){
            std::cerr << "split_bins: ERROR: " << inName << " has no Int_t branch " << polBranch
                      << " to select the polarization orientation" << std::endl;
            gSystem->Exit( 1 );
        }
        t->SetBranchAddress( polBranch, &pol );
    }

    const int kMax = 32;
    Int_t nfs = 0;
    Float_t E[kMax], Px[kMax], Py[kMax], Pz[kMax], eB, pxB, pyB, pzB;
    t->SetBranchAddress( "NumFinalState", &nfs );
    t->SetBranchAddress( "E_FinalState",  E  );
    t->SetBranchAddress( "Px_FinalState", Px );
    t->SetBranchAddress( "Py_FinalState", Py );
    t->SetBranchAddress( "Pz_FinalState", Pz );
    t->SetBranchAddress( "E_Beam",  &eB  );
    t->SetBranchAddress( "Px_Beam", &pxB );
    t->SetBranchAddress( "Py_Beam", &pyB );
    t->SetBranchAddress( "Pz_Beam", &pzB );

    std::vector<TFile*> files;
    std::vector<TTree*> trees;
    for( int i = 0; i < nMass; ++i ){
        for( int j = 0; j < nT; ++j ){
            TFile* f = new TFile( Form( "%s_m%d_t%d.root", outBase, i, j ), "RECREATE" );
            TTree* c = t->CloneTree( 0 );
            c->SetDirectory( f );
            files.push_back( f );
            trees.push_back( c );
        }
    }

    double mStep = ( mHigh - mLow ) / nMass;
    Long64_t n = t->GetEntries(), kept = 0, outMass = 0, outT = 0, tooFew = 0, otherPol = 0;
    double seenMlo = 1e30, seenMhi = -1e30, seenTlo = 1e30, seenThi = -1e30;
    int maxIdx = 0;
    for( int k : idx ) maxIdx = std::max( maxIdx, k );
    for( Long64_t e = 0; e < n; ++e ){
        t->GetEntry( e );
        if( usePol && pol != polValue ){ ++otherPol; continue; }
        if( nfs <= maxIdx || nfs > kMax ){ ++tooFew; continue; }
        TLorentzVector X;
        for( int k : idx ) if( k < nfs ) X += TLorentzVector( Px[k], Py[k], Pz[k], E[k] );
        double m = X.M();
        double mt = -( TLorentzVector( pxB, pyB, pzB, eB ) - X ).M2();
        seenMlo = std::min( seenMlo, m );   seenMhi = std::max( seenMhi, m );
        seenTlo = std::min( seenTlo, mt );  seenThi = std::max( seenThi, mt );
        if( m < mLow || m >= mHigh ){ ++outMass; continue; }
        int i = int( ( m - mLow ) / mStep );
        if( i >= nMass ) i = nMass - 1;

        int j = 0;
        if( useT ){
            if( mt < tLow || mt >= tHigh ){ ++outT; continue; }
            j = tLog ? int( std::log( mt / tLow ) / std::log( tHigh / tLow ) * nT )
                     : int( ( mt - tLow ) / ( tHigh - tLow ) * nT );
            if( j >= nT ) j = nT - 1;
        }
        trees[ i * nT + j ]->Fill();
        ++kept;
    }

    for( size_t b = 0; b < files.size(); ++b ){
        files[b]->cd();
        trees[b]->Write();
        files[b]->Close();
        delete files[b];
    }
    std::cout << "split_bins: " << inName << ": kept " << kept << " of " << n
              << " events in " << nMass << " x " << nT << " bins" << std::endl;
    if( usePol )
        std::cout << "  dropped: " << otherPol << " with " << polBranch << " != " << polValue << std::endl;
    std::cout << "  dropped: " << outMass << " outside the mass range, " << outT
              << " outside the -t range, " << tooFew << " with too few final-state particles"
              << " (NumFinalState <= " << maxIdx << ")" << std::endl;
    if( n > tooFew + otherPol )
        std::cout << "  seen in the file: mass " << seenMlo << " to " << seenMhi
                  << " GeV, -t " << seenTlo << " to " << seenThi << " GeV^2"
                  << " (mass of final-state particles " << massIdx << ")" << std::endl;
    if( kept == 0 && n > 0 )
        std::cerr << "split_bins: WARNING: no events kept from " << inName
                  << "; compare the ranges above with your mass and t binning" << std::endl;
}
