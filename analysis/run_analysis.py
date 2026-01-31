#!/usr/bin/env python3
"""
run_analysis.py - Main analysis pipeline for hippocampal-cortical gamma synchrony

This script implements the falsification-forward analysis described in:
"Robust Hippocampal-Cortical Gamma Synchrony in Mouse Visual Processing"

Usage:
    python run_analysis.py [--output-dir OUTPUT_DIR] [--n-surrogates N] [--seed SEED]

The analysis:
1. Loads Allen Visual Coding Neuropixels sessions with CA1 + VISp coverage
2. Extracts 180s of LFP, split into 3x60s segments
3. Computes PLV and wPLI for gamma band (30-80 Hz)
4. Tests against 5000 phase-randomized surrogates
5. Combines evidence via Fisher's method
6. Classifies sessions as surviving or failing falsification criteria

Output: all_results.json with complete session-level data
"""

import argparse
import json
import logging
from pathlib import Path
from datetime import datetime

import numpy as np
from scipy import signal, stats
from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# =============================================================================
# ANALYSIS PARAMETERS
# =============================================================================

PARAMS = {
    'freq_low': 30,           # Hz - gamma band lower bound
    'freq_high': 80,          # Hz - gamma band upper bound
    'segment_duration': 60,   # seconds per segment
    'n_segments': 3,          # number of segments for temporal replication
    'n_surrogates': 5000,     # surrogates for permutation testing
    'alpha': 0.05,            # significance threshold
    'min_segments_pass': 2,   # minimum segments passing for survival
    'filter_order': 4,        # Butterworth filter order
    'random_seed': 7829,      # base random seed for reproducibility
}


# =============================================================================
# SIGNAL PROCESSING
# =============================================================================

def bandpass_filter(data, fs, low, high, order=4):
    """Apply zero-phase Butterworth bandpass filter."""
    nyq = fs / 2
    sos = signal.butter(order, [low/nyq, high/nyq], btype='band', output='sos')
    return signal.sosfiltfilt(sos, data)


def compute_plv(phase1, phase2):
    """Compute Phase Locking Value between two phase time series."""
    return np.abs(np.mean(np.exp(1j * (phase1 - phase2))))


def compute_wpli(sig1, sig2, fs):
    """
    Compute weighted Phase Lag Index.
    
    wPLI down-weights zero-lag correlations, making it robust to volume conduction.
    """
    # Cross-spectrum via Hilbert transform
    analytic1 = signal.hilbert(sig1)
    analytic2 = signal.hilbert(sig2)
    cross_spectrum = analytic1 * np.conj(analytic2)
    
    # wPLI computation
    imag_cross = np.imag(cross_spectrum)
    numerator = np.abs(np.mean(np.abs(imag_cross) * np.sign(imag_cross)))
    denominator = np.mean(np.abs(imag_cross))
    
    if denominator == 0:
        return 0.0
    return numerator / denominator


def generate_surrogate(data, rng):
    """
    Generate phase-randomized surrogate via FFT method.
    
    Preserves amplitude spectrum while randomizing phase.
    """
    n = len(data)
    fft_data = np.fft.fft(data)
    amplitudes = np.abs(fft_data)
    
    # Random phases (preserve conjugate symmetry for real output)
    random_phases = rng.uniform(0, 2*np.pi, n)
    if n % 2 == 0:
        random_phases[n//2] = 0  # Nyquist
    random_phases[0] = 0  # DC
    # Conjugate symmetry
    random_phases[n//2+1:] = -random_phases[1:n//2][::-1]
    
    surrogate_fft = amplitudes * np.exp(1j * random_phases)
    return np.real(np.fft.ifft(surrogate_fft))


# =============================================================================
# STATISTICAL ANALYSIS
# =============================================================================

def permutation_test(observed, surrogates):
    """
    Compute p-value from permutation test.
    
    Uses (1 + k) / (1 + N) formulation to avoid p=0.
    """
    k = np.sum(surrogates >= observed)
    n = len(surrogates)
    return (1 + k) / (1 + n)


def fisher_combine(p_values):
    """
    Combine p-values using Fisher's method.
    
    χ² = -2 Σ ln(pᵢ), distributed as χ²(2k) under null.
    """
    p_values = np.array(p_values)
    # Clip to avoid log(0)
    p_values = np.clip(p_values, 1e-300, 1.0)
    chi2 = -2 * np.sum(np.log(p_values))
    df = 2 * len(p_values)
    return stats.chi2.sf(chi2, df)


# =============================================================================
# DATA LOADING
# =============================================================================

def get_qualifying_sessions(cache):
    """Find sessions with both CA1 and VISp LFP coverage."""
    sessions = cache.get_session_table()
    
    qualifying = []
    for session_id in sessions.index:
        try:
            channels = cache.get_channels()
            session_channels = channels[channels['ecephys_session_id'] == session_id]
            structures = session_channels['ecephys_structure_acronym'].unique()
            
            if 'CA1' in structures and 'VISp' in structures:
                qualifying.append(session_id)
        except Exception as e:
            logger.warning(f"Error checking session {session_id}: {e}")
    
    return qualifying


def load_session_data(cache, session_id, duration_s=180):
    """
    Load LFP data for a session.
    
    Returns:
        ca1_lfp: LFP from CA1 channel
        visp_lfp: LFP from VISp channel
        fs: sampling frequency
        channel_info: dict with channel IDs
    """
    session = cache.get_session_data(session_id)
    
    # Get channel info
    channels = session.channels
    ca1_channels = channels[channels['ecephys_structure_acronym'] == 'CA1']
    visp_channels = channels[channels['ecephys_structure_acronym'] == 'VISp']
    
    if len(ca1_channels) == 0 or len(visp_channels) == 0:
        raise ValueError(f"Session {session_id} missing CA1 or VISp channels")
    
    # Select middle channel from each region (deterministic)
    ca1_idx = len(ca1_channels) // 2
    visp_idx = len(visp_channels) // 2
    
    ca1_channel_id = ca1_channels.index[ca1_idx]
    visp_channel_id = visp_channels.index[visp_idx]
    
    # Get probe IDs
    ca1_probe_id = ca1_channels.loc[ca1_channel_id, 'ecephys_probe_id']
    visp_probe_id = visp_channels.loc[visp_channel_id, 'ecephys_probe_id']
    
    # Load LFP
    ca1_lfp_obj = session.get_lfp(ca1_probe_id)
    visp_lfp_obj = session.get_lfp(visp_probe_id)
    
    fs = 1250  # Allen LFP sampling rate
    n_samples = int(duration_s * fs)
    
    ca1_lfp = ca1_lfp_obj.sel(channel=ca1_channel_id).values[:n_samples]
    visp_lfp = visp_lfp_obj.sel(channel=visp_channel_id).values[:n_samples]
    
    channel_info = {
        'ca1_channel_id': int(ca1_channel_id),
        'visp_channel_id': int(visp_channel_id),
        'ca1_probe_id': int(ca1_probe_id),
        'visp_probe_id': int(visp_probe_id),
    }
    
    return ca1_lfp, visp_lfp, fs, channel_info


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def analyze_segment(ca1_seg, visp_seg, fs, n_surrogates, rng):
    """
    Analyze a single segment.
    
    Returns dict with PLV, wPLI, and p-values.
    """
    # Bandpass filter
    ca1_filt = bandpass_filter(ca1_seg, fs, PARAMS['freq_low'], PARAMS['freq_high'])
    visp_filt = bandpass_filter(visp_seg, fs, PARAMS['freq_low'], PARAMS['freq_high'])
    
    # Extract phase
    ca1_phase = np.angle(signal.hilbert(ca1_filt))
    visp_phase = np.angle(signal.hilbert(visp_filt))
    
    # Observed metrics
    observed_plv = compute_plv(ca1_phase, visp_phase)
    observed_wpli = compute_wpli(ca1_filt, visp_filt, fs)
    
    # Surrogate testing
    surrogate_plv = np.zeros(n_surrogates)
    surrogate_wpli = np.zeros(n_surrogates)
    
    for i in range(n_surrogates):
        surr_ca1 = generate_surrogate(ca1_filt, rng)
        surr_phase = np.angle(signal.hilbert(surr_ca1))
        
        surrogate_plv[i] = compute_plv(surr_phase, visp_phase)
        surrogate_wpli[i] = compute_wpli(surr_ca1, visp_filt, fs)
    
    # P-values
    plv_p = permutation_test(observed_plv, surrogate_plv)
    wpli_p = permutation_test(observed_wpli, surrogate_wpli)
    
    return {
        'plv': float(observed_plv),
        'wpli': float(observed_wpli),
        'plv_p': float(plv_p),
        'wpli_p': float(wpli_p),
        'pass_plv': plv_p < PARAMS['alpha'],
        'pass_wpli': wpli_p < PARAMS['alpha'],
        'pass_both': plv_p < PARAMS['alpha'] and wpli_p < PARAMS['alpha'],
    }


def analyze_session(cache, session_id, base_seed):
    """
    Complete analysis for one session.
    
    Returns dict with all results.
    """
    logger.info(f"Analyzing session {session_id}")
    
    # Session-specific seed for reproducibility
    session_seed = base_seed + session_id % 10000
    rng = np.random.default_rng(session_seed)
    
    # Load data
    try:
        ca1_lfp, visp_lfp, fs, channel_info = load_session_data(cache, session_id)
    except Exception as e:
        logger.error(f"Failed to load session {session_id}: {e}")
        return None
    
    # Segment analysis
    segment_samples = int(PARAMS['segment_duration'] * fs)
    segment_results = []
    
    for seg_idx in range(PARAMS['n_segments']):
        start = seg_idx * segment_samples
        end = start + segment_samples
        
        ca1_seg = ca1_lfp[start:end]
        visp_seg = visp_lfp[start:end]
        
        result = analyze_segment(ca1_seg, visp_seg, fs, PARAMS['n_surrogates'], rng)
        result['segment'] = seg_idx
        segment_results.append(result)
        
        logger.info(f"  Segment {seg_idx}: PLV={result['plv']:.3f} (p={result['plv_p']:.4f}), "
                   f"wPLI={result['wpli']:.3f} (p={result['wpli_p']:.4f})")
    
    # Combine evidence
    plv_p_values = [s['plv_p'] for s in segment_results]
    wpli_p_values = [s['wpli_p'] for s in segment_results]
    
    combined_plv_p = fisher_combine(plv_p_values)
    combined_wpli_p = fisher_combine(wpli_p_values)
    
    # Count segments passing both criteria
    segments_pass_both = sum(s['pass_both'] for s in segment_results)
    
    # Survival determination
    survives = (
        combined_plv_p < PARAMS['alpha'] and
        combined_wpli_p < PARAMS['alpha'] and
        segments_pass_both >= PARAMS['min_segments_pass']
    )
    
    logger.info(f"  Combined: PLV p={combined_plv_p:.2e}, wPLI p={combined_wpli_p:.2e}, "
               f"Segments={segments_pass_both}/3 -> {'SURVIVES' if survives else 'FAILS'}")
    
    return {
        'session_id': int(session_id),
        'segments': segment_results,
        'combined_plv_p': float(combined_plv_p),
        'combined_wpli_p': float(combined_wpli_p),
        'segments_pass_both': int(segments_pass_both),
        'survives': bool(survives),
        **channel_info,
    }


def main():
    parser = argparse.ArgumentParser(description='Hippocampal-cortical gamma synchrony analysis')
    parser.add_argument('--output-dir', type=Path, default=Path('results'),
                       help='Output directory')
    parser.add_argument('--n-surrogates', type=int, default=5000,
                       help='Number of surrogates for permutation testing')
    parser.add_argument('--seed', type=int, default=7829,
                       help='Random seed for reproducibility')
    parser.add_argument('--manifest', type=Path, default=None,
                       help='Path to Allen SDK manifest (optional)')
    args = parser.parse_args()
    
    # Update params
    PARAMS['n_surrogates'] = args.n_surrogates
    PARAMS['random_seed'] = args.seed
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize cache
    logger.info("Initializing Allen SDK cache...")
    if args.manifest:
        cache = EcephysProjectCache.from_warehouse(manifest=str(args.manifest))
    else:
        cache = EcephysProjectCache.from_warehouse()
    
    # Find qualifying sessions
    logger.info("Finding sessions with CA1 + VISp coverage...")
    session_ids = get_qualifying_sessions(cache)
    logger.info(f"Found {len(session_ids)} qualifying sessions")
    
    # Analyze all sessions
    results = {
        'metadata': {
            'analysis_date': datetime.now().isoformat(),
            'parameters': PARAMS,
            'n_sessions': len(session_ids),
        },
        'sessions': {},
        'survivors': [],
    }
    
    for session_id in session_ids:
        session_result = analyze_session(cache, session_id, PARAMS['random_seed'])
        
        if session_result is not None:
            results['sessions'][str(session_id)] = session_result
            if session_result['survives']:
                results['survivors'].append(session_id)
    
    # Summary
    n_analyzed = len(results['sessions'])
    n_survivors = len(results['survivors'])
    logger.info(f"\n{'='*60}")
    logger.info(f"ANALYSIS COMPLETE")
    logger.info(f"Sessions analyzed: {n_analyzed}")
    logger.info(f"Survivors: {n_survivors} ({100*n_survivors/n_analyzed:.1f}%)")
    logger.info(f"{'='*60}")
    
    # Save results
    output_file = args.output_dir / 'all_results.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    logger.info(f"Results saved to {output_file}")
    
    # Also save channel IDs as CSV for easy reference
    import pandas as pd
    channel_data = []
    for sid, r in results['sessions'].items():
        channel_data.append({
            'session_id': sid,
            'ca1_channel_id': r['ca1_channel_id'],
            'visp_channel_id': r['visp_channel_id'],
            'ca1_probe_id': r['ca1_probe_id'],
            'visp_probe_id': r['visp_probe_id'],
            'survives': r['survives'],
        })
    df = pd.DataFrame(channel_data)
    df.to_csv(args.output_dir / 'channel_ids.csv', index=False)
    logger.info(f"Channel IDs saved to {args.output_dir / 'channel_ids.csv'}")


if __name__ == '__main__':
    main()
