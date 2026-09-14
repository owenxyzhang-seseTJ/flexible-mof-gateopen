#!/usr/bin/env python3
"""Run GFN2-xTB in its own process to avoid mixed torch/tblite OpenMP runtimes."""
import argparse
from pathlib import Path
from prepare_guest_rdkit import refine_tblite
if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir",required=True,type=Path)
    args=p.parse_args()
    print(refine_tblite(args.output_dir))
