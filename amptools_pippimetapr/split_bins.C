// split_bins.C
// ROOT macro: split a ROOTDataReader-format tree into mass x (-t) bins in one pass.
// Used by amptools_auto.py in place of split_mass / split_t, because those
// compute the mass from every final-state particle after the recoil, which for
// "reaction EtaPi Beam Proton EtaPrime Pi- Pi+" would be M(eta' pi- pi+).
//
//   mass = M of the final-state particles listed in massIdx
//          (default "1,2" = EtaPrime, Pi-; the proton is final-state index 0)
//   -t   = -(beam - those particles)^2, i.e. the momentum transfer to the
//          recoiling p pi+ system (does not depend on the target mass)
//
// Writes <outBase>_m<i>_t<j>.root for every bin (empty bins give an empty tree),
// i = 0..nMass-1, j = 0..nT-1 (j = 0 only when nT <= 0, i.e. no t binning).
// Every branch is copied unchanged (Weight etc. included).  Bins are
// low <= x < high; events outside the ranges are dropped.
//
// Usage:
//   root -l -b -q 'split_bins.C("in.root","out/data",1.1,2.6,30, 0.1,0.3,1, false)'
//   root -l -b -q 'split_bins.C("in.root","out/data",1.1,2.6,30, 0,0,0)'      // mass only
//   arguments: in, outBase, massLow, massHigh, nMass, tLow, tHigh, nT,
//              tLog (logarithmic t bins), tree ("kin"), massIdx ("1,2")

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
                 const char* treeName = "kin", const char* massIdx = "1,2" )
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
    Long64_t n = t->GetEntries(), kept = 0, outMass = 0, outT = 0, tooFew = 0;
    double seenMlo = 1e30, seenMhi = -1e30, seenTlo = 1e30, seenThi = -1e30;
    int maxIdx = 0;
    for( int k : idx ) maxIdx = std::max( maxIdx, k );
    for( Long64_t e = 0; e < n; ++e ){
        t->GetEntry( e );
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
    std::cout << "  dropped: " << outMass << " outside the mass range, " << outT
              << " outside the -t range, " << tooFew << " with too few final-state particles"
              << " (NumFinalState <= " << maxIdx << ")" << std::endl;
    if( n > tooFew )
        std::cout << "  seen in the file: mass " << seenMlo << " to " << seenMhi
                  << " GeV, -t " << seenTlo << " to " << seenThi << " GeV^2"
                  << " (mass of final-state particles " << massIdx << ")" << std::endl;
    if( kept == 0 && n > 0 )
        std::cerr << "split_bins: WARNING: no events kept from " << inName
                  << "; compare the ranges above with your mass and t binning" << std::endl;
}
