#!/usr/bin/env python3
"""
verify_results.py - Quick verification that results can be reproduced

This script runs a sanity check on 2 sessions to verify the analysis
produces consistent results without running the full 10-hour pipeline.
"""

import json
import sys
from pathlib import Path

def main():
    results_file = Path('results/all_results.json')
    
    if not results_file.exists():
        print("ERROR: results/all_results.json not found")
        print("Run the full analysis first: python analysis/run_analysis.py")
        sys.exit(1)
    
    with open(results_file) as f:
        results = json.load(f)
    
    # Expected values from the paper
    EXPECTED = {
        'n_sessions': 49,
        'n_survivors': 29,
        'survival_rate': 0.592,
        # Two example sessions with known results
        'session_794812542': {
            'survives': True,
            'combined_wpli_p': 2.82e-09,  # Surrogate floor
        },
        'session_829720705': {
            'survives': True,
            'combined_wpli_p': 2.82e-09,
        },
    }
    
    print("=" * 60)
    print("VERIFICATION RESULTS")
    print("=" * 60)
    
    all_pass = True
    
    # Check counts
    n_sessions = len(results['sessions'])
    n_survivors = len(results['survivors'])
    
    print(f"\nSessions analyzed: {n_sessions} (expected: {EXPECTED['n_sessions']})")
    if n_sessions != EXPECTED['n_sessions']:
        print("  ❌ MISMATCH")
        all_pass = False
    else:
        print("  ✓ Match")
    
    print(f"\nSurvivors: {n_survivors} (expected: {EXPECTED['n_survivors']})")
    if n_survivors != EXPECTED['n_survivors']:
        print("  ❌ MISMATCH")
        all_pass = False
    else:
        print("  ✓ Match")
    
    survival_rate = n_survivors / n_sessions
    print(f"\nSurvival rate: {100*survival_rate:.1f}% (expected: {100*EXPECTED['survival_rate']:.1f}%)")
    if abs(survival_rate - EXPECTED['survival_rate']) > 0.01:
        print("  ❌ MISMATCH")
        all_pass = False
    else:
        print("  ✓ Match")
    
    # Check specific sessions
    for session_key in ['session_794812542', 'session_829720705']:
        session_id = session_key.split('_')[1]
        print(f"\n{session_key}:")
        
        if session_id not in results['sessions']:
            print(f"  ❌ Session not found")
            all_pass = False
            continue
        
        session = results['sessions'][session_id]
        expected = EXPECTED[session_key]
        
        if session['survives'] != expected['survives']:
            print(f"  ❌ survives: {session['survives']} (expected: {expected['survives']})")
            all_pass = False
        else:
            print(f"  ✓ survives: {session['survives']}")
        
        # p-value should be at floor (allow some tolerance due to float precision)
        if session['combined_wpli_p'] > 1e-7:
            print(f"  ⚠ combined_wpli_p: {session['combined_wpli_p']:.2e} (expected floor)")
        else:
            print(f"  ✓ combined_wpli_p at floor")
    
    print("\n" + "=" * 60)
    if all_pass:
        print("✅ ALL CHECKS PASSED")
        print("Results match published values")
    else:
        print("❌ SOME CHECKS FAILED")
        print("Results may differ from published values")
    print("=" * 60)
    
    sys.exit(0 if all_pass else 1)

if __name__ == '__main__':
    main()
