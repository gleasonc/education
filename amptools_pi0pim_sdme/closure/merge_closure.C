// Merge gen_amp output files into one AmpTools tree, adding the branches the
// SDME fit chain expects: pol/I (orientation) and Weight/F (= 1).
//
//   root -l -b -q 'merge_closure.C("data.root","data_files.txt","pol")'
//
// listFile: one "file.root pol" per line (pol = -1 for MC without orientation).
// gen_amp writes the tree "kin" with Beam + Proton Pi- Pi0 (NumFinalState = 3).

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "TFile.h"
#include "TTree.h"

void merge_closure(const char* outName, const char* listFile, const char* polBranch = "pol")
{
  const int maxN = 10;
  Int_t nFS;
  Float_t e[maxN], px[maxN], py[maxN], pz[maxN];
  Float_t eB, pxB, pyB, pzB, wIn, w;
  Int_t pol;

  TFile out(outName, "RECREATE");
  TTree* tOut = new TTree("kin", "kin");
  tOut->Branch("NumFinalState", &nFS, "NumFinalState/I");
  tOut->Branch("E_FinalState", e, "E_FinalState[NumFinalState]/F");
  tOut->Branch("Px_FinalState", px, "Px_FinalState[NumFinalState]/F");
  tOut->Branch("Py_FinalState", py, "Py_FinalState[NumFinalState]/F");
  tOut->Branch("Pz_FinalState", pz, "Pz_FinalState[NumFinalState]/F");
  tOut->Branch("E_Beam", &eB, "E_Beam/F");
  tOut->Branch("Px_Beam", &pxB, "Px_Beam/F");
  tOut->Branch("Py_Beam", &pyB, "Py_Beam/F");
  tOut->Branch("Pz_Beam", &pzB, "Pz_Beam/F");
  tOut->Branch("Weight", &w, "Weight/F");
  tOut->Branch(polBranch, &pol, (std::string(polBranch) + "/I").c_str());

  std::ifstream list(listFile);
  std::string line;
  Long64_t total = 0;
  while (std::getline(list, line)) {
    std::istringstream ss(line);
    std::string file;
    int p;
    if (!(ss >> file >> p)) continue;
    TFile in(file.c_str());
    TTree* t = in.IsZombie() ? nullptr : (TTree*)in.Get("kin");
    if (!t) { std::cerr << "merge_closure: no tree kin in " << file << std::endl; exit(1); }
    t->SetBranchAddress("NumFinalState", &nFS);
    t->SetBranchAddress("E_FinalState", e);
    t->SetBranchAddress("Px_FinalState", px);
    t->SetBranchAddress("Py_FinalState", py);
    t->SetBranchAddress("Pz_FinalState", pz);
    t->SetBranchAddress("E_Beam", &eB);
    t->SetBranchAddress("Px_Beam", &pxB);
    t->SetBranchAddress("Py_Beam", &pyB);
    t->SetBranchAddress("Pz_Beam", &pzB);
    bool hasW = t->GetBranch("Weight") != nullptr;
    if (hasW) t->SetBranchAddress("Weight", &wIn);
    pol = p;
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (nFS > maxN) { std::cerr << "merge_closure: too many particles" << std::endl; exit(1); }
      w = hasW ? wIn : 1.f;
      out.cd();
      tOut->Fill();
    }
    total += t->GetEntries();
    std::cout << file << ": " << t->GetEntries() << " events, " << polBranch << " = " << p << std::endl;
  }
  out.cd();
  tOut->Write();
  out.Close();
  std::cout << "Wrote " << outName << " (" << total << " events)" << std::endl;
}
