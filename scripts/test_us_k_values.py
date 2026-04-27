#!/usr/bin/env python3
"""Test different force constants (k) for umbrella sampling.

Runs short 5ns simulations for a single window with different k values
and compares CV distributions.

Usage:
    python scripts/test_us_k_values.py \
        --prmtop data/md_runs/Hsap_WT/Hsap_WT.prmtop \
        --rst7 data/analysis/final/us_start_Hsap_WT.rst7 \
        --center 10.0 \
        --k-values 200 400 1000 \
        --gpu 0 \
        --outdir data/analysis/us_k_test
"""
import argparse
import os
import json
import numpy as np
from pathlib import Path

# Import run_us from run_us_simple.py
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_us_simple import run_us


def analyze_cv_distribution(cv_file):
    """Analyze CV distribution from a .cv.dat file."""
    data = np.loadtxt(cv_file, comments='#', skiprows=1)
    if len(data) == 0:
        return None
    # Skip first 20% as warmup
    warmup = int(len(data) * 0.2)
    cv_vals = data[warmup:, 2]  # CV column
    return {
        'mean': float(np.mean(cv_vals)),
        'std': float(np.std(cv_vals)),
        'min': float(np.min(cv_vals)),
        'max': float(np.max(cv_vals)),
        'n_frames': len(cv_vals),
    }


def main():
    parser = argparse.ArgumentParser(description='Test US force constants')
    parser.add_argument('--prmtop', required=True)
    parser.add_argument('--rst7', required=True)
    parser.add_argument('--center', type=float, required=True)
    parser.add_argument('--k-values', nargs='+', type=float, default=[200, 400, 1000])
    parser.add_argument('--gpu', type=int, default=0)
    parser.add_argument('--prod-ns', type=float, default=5.0)
    parser.add_argument('--outdir', default='data/analysis/us_k_test')
    args = parser.parse_args()
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    results = {}
    
    for k in args.k_values:
        name = f"w{args.center:.0f}A_k{int(k)}"
        print(f"\n{'='*60}")
        print(f"Testing k={k} for center={args.center}Å")
        print(f"{'='*60}")
        
        run_us(
            prmtop_path=args.prmtop,
            rst7_path=args.rst7,
            center_A=args.center,
            k=k,
            name=name,
            outdir=str(outdir),
            gpu_id=args.gpu,
            em_steps=500,
            prod_ns=args.prod_ns,
        )
        
        cv_file = outdir / f"{name}_cv.dat"
        stats = analyze_cv_distribution(str(cv_file))
        if stats:
            results[f'k{int(k)}'] = stats
            print(f"\n  CV distribution: μ={stats['mean']:.2f}Å, σ={stats['std']:.2f}Å, range=[{stats['min']:.2f}, {stats['max']:.2f}]")
    
    # Summary comparison
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"{'k (kJ/mol/nm²)':>16} {'mean(Å)':>10} {'std(Å)':>10} {'range(Å)':>15}")
    for k_label, stats in results.items():
        k_val = int(k_label[1:])
        r = stats['max'] - stats['min']
        print(f"{k_val:>16} {stats['mean']:>10.2f} {stats['std']:>10.2f} {r:>15.2f}")
    
    # Recommendation
    k_list = [int(k[1:]) for k in results.keys()]
    stds = [results[f'k{k}']['std'] for k in k_list]
    
    print(f"\nRecommendation:")
    if all(s < 0.5 for s in stds):
        print("  All k values produce very narrow distributions. Consider reducing k to allow more exploration.")
    elif stds[0] > stds[-1] * 1.5:
        print("  Lower k values allow significantly wider exploration. k=200-400 may be preferable.")
    else:
        print("  CV distributions are similar across k values. Current k=1000 is acceptable.")
    
    # Save results
    summary_path = outdir / 'k_test_summary.json'
    with open(summary_path, 'w') as f:
        json.dump({
            'center': args.center,
            'prod_ns': args.prod_ns,
            'results': results,
        }, f, indent=2)
    print(f"\nSaved summary: {summary_path}")


if __name__ == '__main__':
    main()
