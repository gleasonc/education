// extract_fit_results_etapi.cc
//
// Reads a list of AmpTools .fit files (one path per line) for the eta pi- fits
// defined in fit_etaprpim_SPD_TEMPLATE.cfg and writes one csv row per fit:
// per-wave intensities and production parameters, reflectivity totals, coherent
// sums over m for each L (Ssum+, Psum+, Dsum+, Ssum-, ..., and Ssum, Psum, Dsum
// over both reflectivities) and phase differences.
//
// Usage:
//   extract_fit_results_etapi <file_list> <out.csv> <accCorrected 0|1>
//                             [<phase pairs | auto>] [<cov.csv | ->] [<corr.csv | ->]
//
//   phase pairs: comma separated, e.g. "P1+:D1+,S0+:D1+"; "auto" (default) uses
//                every pair of non-zero waves within the same reflectivity,
//                decided from the first fit file
//   cov/corr:    long-format csv (file,par_i,par_j,value); "-" to skip
//
// Build (same way as plot_etapi_delta), e.g.:
//   g++ -O2 extract_fit_results_etapi.cc -I$AMPTOOLS -L$AMPTOOLS/lib
//       -lAmpTools $(root-config --cflags --libs) -o extract_fit_results_etapi

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <complex>
#include <cmath>
#include <limits>
#include <algorithm>
#include <utility>

#include "IUAmpTools/FitResults.h"

using namespace std;

// waves in fit_etaprpim_SPD_TEMPLATE.cfg
// naming: <L><m><reflectivity>, 'n' = negative m (Dn2+ is L=2, m=-2, eps=+)
static const char* kWaves[] = {
    "S0+", "Pn1+", "P0+", "P1+", "Dn2+", "Dn1+", "D0+", "D1+", "D2+",
    "S0-", "Pn1-", "P0-", "P1-", "Dn2-", "Dn1-", "D0-", "D1-", "D2-"
};
static const int kNumWaves = sizeof(kWaves) / sizeof(kWaves[0]);
static const double kNaN = numeric_limits<double>::quiet_NaN();

bool isPos( const string& w ){ return w[ w.size() - 1 ] == '+'; }

// both sums carrying a wave (Re first)
vector< string > waveAmps( const string& w ){
    string refl = isPos( w ) ? "Positive" : "Negative";
    vector< string > amps;
    amps.push_back( "EtaPi::" + refl + "Re::" + w );
    amps.push_back( "EtaPi::" + refl + "Im::" + w );
    return amps;
}

bool has( const vector< string >& list, const string& s ){
    return find( list.begin(), list.end(), s ) != list.end();
}

pair< double, double > safeIntensity( const FitResults& r, const vector< string >& fitAmps,
                                      const vector< string >& amps, bool acc ){
    vector< string > present;
    for( size_t a = 0; a < amps.size(); ++a )
        if( has( fitAmps, amps[a] ) ) present.push_back( amps[a] );
    if( present.empty() ) return make_pair( 0., 0. );
    return r.intensity( present, acc );
}

vector< string > split( const string& s, char d ){
    vector< string > out; string tok; istringstream ss( s );
    while( getline( ss, tok, d ) ) if( !tok.empty() ) out.push_back( tok );
    return out;
}

int main( int argc, char* argv[] ){

    if( argc < 4 ){
        cerr << "Usage: " << argv[0] << " <file_list> <out.csv> <accCorrected 0|1>"
             << " [phase pairs|auto] [cov.csv|-] [corr.csv|-]" << endl;
        return 1;
    }
    string listName = argv[1], outName = argv[2];
    bool acc = atoi( argv[3] ) != 0;
    string phaseArg = argc > 4 ? argv[4] : "auto";
    string covName  = argc > 5 ? argv[5] : "-";
    string corrName = argc > 6 ? argv[6] : "-";

    vector< string > files;
    { ifstream in( listName.c_str() ); string l;
      while( getline( in, l ) ) if( !l.empty() ) files.push_back( l ); }
    if( files.empty() ){ cerr << "No input files in " << listName << endl; return 1; }

    vector< string > waves( kWaves, kWaves + kNumWaves );

    vector< string > posAmps, negAmps, allAmps;
    for( int w = 0; w < kNumWaves; ++w ){
        vector< string > a = waveAmps( waves[w] );
        vector< string >& r = isPos( waves[w] ) ? posAmps : negAmps;
        r.insert( r.end(), a.begin(), a.end() );
        allAmps.insert( allAmps.end(), a.begin(), a.end() );
    }

    // coherent sums over m for each L: per reflectivity ("S+", ...) and both ("S", ...)
    const char* kL[] = { "S", "P", "D" };
    vector< string > lSumNames;
    vector< vector< string > > lSumAmps;
    for( int r = 0; r < 3; ++r ){            // r = 0: +, 1: -, 2: both
        for( int l = 0; l < 3; ++l ){
            vector< string > amps;
            for( int w = 0; w < kNumWaves; ++w ){
                if( waves[w][0] != kL[l][0] ) continue;
                if( r == 0 && !isPos( waves[w] ) ) continue;
                if( r == 1 &&  isPos( waves[w] ) ) continue;
                vector< string > a = waveAmps( waves[w] );
                amps.insert( amps.end(), a.begin(), a.end() );
            }
            lSumNames.push_back( string( kL[l] ) + "sum" + ( r == 0 ? "+" : r == 1 ? "-" : "" ) );
            lSumAmps.push_back( amps );
        }
    }

    // ---- decide phase-difference pairs ----
    vector< pair< string, string > > phases;
    if( phaseArg == "auto" ){
        FitResults first( files[0] );
        vector< string > fa = first.ampList();
        vector< string > live;
        for( int w = 0; w < kNumWaves; ++w ){
            string re = waveAmps( waves[w] )[0];
            if( has( fa, re ) && abs( first.productionParameter( re ) ) > 0 ) live.push_back( waves[w] );
        }
        for( size_t i = 0; i < live.size(); ++i )
            for( size_t j = i + 1; j < live.size(); ++j )
                if( isPos( live[i] ) == isPos( live[j] ) ) phases.push_back( make_pair( live[i], live[j] ) );
    }
    else if( phaseArg != "none" ){
        vector< string > pairs = split( phaseArg, ',' );
        for( size_t p = 0; p < pairs.size(); ++p ){
            vector< string > ab = split( pairs[p], ':' );
            if( ab.size() != 2 || !has( waves, ab[0] ) || !has( waves, ab[1] ) ){
                cerr << "Bad phase pair '" << pairs[p] << "' (use e.g. P1+:D1+)" << endl; return 1;
            }
            if( isPos( ab[0] ) != isPos( ab[1] ) )
                cerr << "Warning: " << pairs[p] << " mixes reflectivities; the phase is not physical" << endl;
            phases.push_back( make_pair( ab[0], ab[1] ) );
        }
    }

    // ---- header ----
    ofstream out( outName.c_str() );
    out << "file,likelihood,minuit_status,ematrix_status";
    for( int w = 0; w < kNumWaves; ++w ){
        const string& n = waves[w];
        out << "," << n << "," << n << "_err," << n << "_re," << n << "_im";
    }
    out << ",pos_refl,pos_refl_err,neg_refl,neg_refl_err,all,all_err";
    for( size_t k = 0; k < lSumNames.size(); ++k )
        out << "," << lSumNames[k] << "," << lSumNames[k] << "_err";
    for( size_t p = 0; p < phases.size(); ++p ){
        string c = phases[p].first + "_" + phases[p].second;
        out << ",phase_" << c << ",phase_" << c << "_err";
    }
    out << "\n";
    out.precision( 10 );

    ofstream cov, corr;
    if( covName  != "-" ){ cov.open(  covName.c_str()  ); cov  << "file,par_i,par_j,value\n"; cov.precision( 10 ); }
    if( corrName != "-" ){ corr.open( corrName.c_str() ); corr << "file,par_i,par_j,value\n"; corr.precision( 10 ); }

    // ---- loop over fits ----
    for( size_t f = 0; f < files.size(); ++f ){
        cout << "Processing " << files[f] << endl;
        FitResults r( files[f] );
        if( !r.valid() ){
            cerr << "Warning: invalid fit results in " << files[f] << ", skipping" << endl;
            continue;
        }
        vector< string > fa = r.ampList();

        out << files[f] << "," << r.likelihood() << ","
            << r.lastMinuitCommandStatus() << "," << r.eMatrixStatus();

        for( int w = 0; w < kNumWaves; ++w ){
            vector< string > a = waveAmps( waves[w] );
            pair< double, double > I = safeIntensity( r, fa, a, acc );
            complex< double > c = has( fa, a[0] ) ? r.productionParameter( a[0] )
                                                  : complex< double >( kNaN, kNaN );
            out << "," << I.first << "," << I.second << "," << c.real() << "," << c.imag();
        }

        pair< double, double > Ip = safeIntensity( r, fa, posAmps, acc );
        pair< double, double > In = safeIntensity( r, fa, negAmps, acc );
        pair< double, double > Ia = safeIntensity( r, fa, allAmps, acc );
        out << "," << Ip.first << "," << Ip.second << "," << In.first << "," << In.second
            << "," << Ia.first << "," << Ia.second;

        for( size_t k = 0; k < lSumAmps.size(); ++k ){
            pair< double, double > Il = safeIntensity( r, fa, lSumAmps[k], acc );
            out << "," << Il.first << "," << Il.second;
        }

        for( size_t p = 0; p < phases.size(); ++p ){
            string a1 = waveAmps( phases[p].first )[0], a2 = waveAmps( phases[p].second )[0];
            if( has( fa, a1 ) && has( fa, a2 ) ){
                pair< double, double > ph = r.phaseDiff( a1, a2 );
                out << "," << ph.first << "," << ph.second;
            }
            else out << "," << kNaN << "," << kNaN;
        }
        out << "\n";

        if( cov.is_open() || corr.is_open() ){
            const vector< string >& pn = r.parNameList();
            const vector< vector< double > >& m = r.errorMatrix();
            for( size_t i = 0; i < pn.size() && i < m.size(); ++i )
                for( size_t j = 0; j < pn.size() && j < m[i].size(); ++j ){
                    if( cov.is_open() )
                        cov << files[f] << "," << pn[i] << "," << pn[j] << "," << m[i][j] << "\n";
                    if( corr.is_open() ){
                        double d = sqrt( m[i][i] * m[j][j] );
                        corr << files[f] << "," << pn[i] << "," << pn[j] << ","
                             << ( d > 0 ? m[i][j] / d : 0. ) << "\n";
                    }
                }
        }
    }

    cout << "Wrote " << outName << endl;
    return 0;
}
