"""Exercise 1: What is inside the file?

Opens the ROOT file with uproot, lists what is in it, prints every branch
(column) of the tree with its type, shows the first 5 entries and counts
entries that pass a few conditions. No plots in this exercise.

Run from the tutorial/ folder:
    python python/ex1_explore.py
    python python/ex1_explore.py /path/to/flat_pipneta__B4_v6.root
"""
import sys
import numpy as np
import uproot

# Use the path given on the command line, or the default one.
filename = sys.argv[1] if len(sys.argv) > 1 else "../flat_pipneta__B4_v6.root"

f = uproot.open(filename)

print("=== Contents of the file ===")
print(f.classnames())           # {name: type} for every object in the file

tree = f["pipneta__B4"]
print("\n=== The tree ===")
print("Tree name:         ", tree.name)
print("Number of entries: ", tree.num_entries)
print("Number of branches:", len(tree.keys()), "\n")

# Each branch with its type: double = decimal number, int32_t = whole number
for name, branch in tree.items():
    print(f"  {name:20s} {branch.typename}")

print("\n=== First 5 entries ===")
cols = ["Beam_E", "Chi2_NDF", "Mandelstam_t", "MissingMassSquared", "accidweight", "pol"]
first = tree.arrays(cols, entry_stop=5, library="np")   # dict of numpy arrays
print("Row " + "".join(f"{c[:12]:>14s}" for c in cols))
for i in range(5):
    print(f"{i:3d} " + "".join(f"{first[c][i]:14.6g}" for c in cols))

print("\n=== Counting with conditions ===")
# Read only the branches we need (reading all 23 would use a lot of memory).
d = tree.arrays(["accidweight", "Chi2_NDF", "pol"], library="np")
# (array == 1) is an array of True/False; summing it counts the Trues.
print("accidweight == 1 (in-time beam photon):", np.sum(d["accidweight"] == 1))
print("accidweight <  0 (out-of-time photon): ", np.sum(d["accidweight"] < 0))
print("Chi2_NDF < 10:                         ", np.sum(d["Chi2_NDF"] < 10))

print("\nEntries for each value of pol:")
values, counts = np.unique(d["pol"], return_counts=True)
for v, n in zip(values, counts):
    print(f"  pol = {v}: {n}")
