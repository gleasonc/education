// =====================================================================
// Exercise 1: What is inside the file?
//
// What this macro does:
//   1. Opens the ROOT file.
//   2. Lists what is stored in it (one TTree called "pipneta__B4").
//   3. Prints the names and types of all the "branches" (columns).
//   4. Prints the first few rows so you can see real numbers.
//   5. Counts how many rows pass a few simple conditions.
//
// How to run it (in a terminal):
//   cd tutorial/cpp
//   root -l -b -q ex1_explore.C
//
// If your data file is somewhere else, give its path in quotes:
//   root -l -b -q 'ex1_explore.C("/path/to/flat_pipneta__B4_v6.root")'
//
// What you should see: a list of 23 branches, 6,091,343 entries,
// and a small table of numbers. No plots are made in this exercise.
// =====================================================================

// "#include" lines tell ROOT which tools (classes) we are going to use.
// TFile = a ROOT file, TTree = a table of events stored in the file.
#include "TFile.h"
#include "TTree.h"
#include <iostream>   // lets us print text with std::cout

// A ROOT macro is a "function". ROOT runs the function whose name is the
// same as the file name (ex1_explore.C -> ex1_explore).
//
// The text in the brackets is the function's "argument": a value you can
// pass in when you run it. Here it is the path to the data file. The part
// after "=" is the default, used when you do not pass anything.
// "const char*" is C++'s way of writing "a piece of text".
void ex1_explore(const char* filename = "../../flat_pipneta__B4_v6.root")
{
   // ---- Step 1: open the file -------------------------------------
   // TFile::Open opens the file and gives us back a POINTER to it.
   // A pointer (the "*" in "TFile*") is like an address: it tells the
   // program where the file object lives in memory. We call it "f".
   TFile* f = TFile::Open(filename);

   // If the file could not be opened (wrong path, typo), f is "nullptr",
   // which means "points to nothing". Check for that and stop early.
   if (f == nullptr) {
      std::cout << "Could not open " << filename << std::endl;
      return;   // "return" leaves the function right here
   }

   // ---- Step 2: list the contents ---------------------------------
   // To use something a pointer points to, write "->" after it.
   // f->ls() means "call the ls (list) function of the file f".
   std::cout << "\n=== Contents of the file ===" << std::endl;
   f->ls();

   // ---- Step 3: get the tree --------------------------------------
   // f->Get("name") finds an object by its name. Get() does not know
   // what kind of object it is, so we tell it with (TTree*) in front.
   TTree* tree = (TTree*)f->Get("pipneta__B4");

   // A tree is a table. Each ROW is one "entry" (one candidate event,
   // also called a combo). Each COLUMN is a "branch" (one number we
   // know about that event, e.g. the beam energy).
   std::cout << "\n=== The tree ===" << std::endl;
   std::cout << "Tree name:         " << tree->GetName() << std::endl;
   std::cout << "Number of entries: " << tree->GetEntries() << std::endl;

   // ---- Step 4: list the branches ---------------------------------
   // tree->GetListOfBranches() gives us the list of columns.
   // GetEntries() on that list tells us how many there are.
   int nBranches = tree->GetListOfBranches()->GetEntries();
   std::cout << "Number of branches: " << nBranches << "\n" << std::endl;

   // A "for loop" repeats the code inside { } several times.
   //   int i = 0        start a counter called i at 0
   //   i < nBranches    keep going while i is less than nBranches
   //   i++              add 1 to i after each pass
   // So the code inside runs for i = 0, 1, 2, ..., nBranches-1.
   for (int i = 0; i < nBranches; i++) {
      // Take branch number i from the list.
      TBranch* br = (TBranch*)tree->GetListOfBranches()->At(i);
      // Each branch's title ends in "/D" (a decimal number, "Double")
      // or "/I" (a whole number, "Integer").
      std::cout << "  " << br->GetTitle() << std::endl;
   }

   // ---- Step 5: look at the first rows ----------------------------
   // tree->Scan(columns, cut, options, how many rows, first row) prints
   // a table. "" for the cut means "no cut, show every row".
   std::cout << "\n=== First 5 entries ===" << std::endl;
   tree->Scan("Beam_E:Chi2_NDF:Mandelstam_t:MissingMassSquared:accidweight:pol",
              "", "colsize=12", 5, 0);

   // ---- Step 6: count entries that pass a condition ---------------
   // tree->GetEntries("condition") counts the rows where the condition
   // is true. The condition is written like C++ (==, <, >, &&, ||).
   std::cout << "\n=== Counting with conditions ===" << std::endl;
   std::cout << "accidweight == 1 (in-time beam photon): "
             << tree->GetEntries("accidweight == 1") << std::endl;
   std::cout << "accidweight <  0 (out-of-time photon):  "
             << tree->GetEntries("accidweight < 0") << std::endl;
   std::cout << "Chi2_NDF < 10:                          "
             << tree->GetEntries("Chi2_NDF < 10") << std::endl;

   // "pol" is the orientation of the beam's polarization in degrees.
   // -1 means amorphous radiator (no polarization).
   std::cout << "\nEntries for each value of pol:" << std::endl;
   int polValues[5] = {-1, 0, 45, 90, 135};   // an "array": a list of 5 numbers
   for (int i = 0; i < 5; i++) {
      // Form() builds a piece of text, like "pol == 45", from a number.
      std::cout << "  pol = " << polValues[i] << ": "
                << tree->GetEntries(Form("pol == %d", polValues[i])) << std::endl;
   }

   // Always close the file when you are done with it.
   f->Close();
}
