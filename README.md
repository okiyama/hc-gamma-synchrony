# Hippocampal-Cortical Gamma Synchrony Analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

Falsification-forward analysis pipeline for detecting inter-regional gamma synchrony, applied to the Allen Institute Visual Coding Neuropixels dataset.

**Paper:** [bioRxiv preprint link]

## Key Finding

59% of sessions (29/49) with simultaneous hippocampal CA1 and primary visual cortex (VISp) recordings exhibit robust gamma-band (30-80 Hz) phase synchrony that survives:
- Volume conduction control (wPLI)
- 5000-surrogate permutation testing
- Temporal replication across 3 independent segments
- Fisher's method for combining evidence

## Repository Structure

```
schizago/
├── analysis/
│   ├── run_analysis.py          # Main analysis pipeline
│   ├── plv_wpli.py              # PLV and wPLI computation
│   ├── surrogate.py             # Phase randomization surrogates
│   └── fisher_combine.py        # Fisher's method implementation
├── figures/
│   ├── generate_figures.py      # Reproduce all paper figures
│   └── *.png, *.pdf             # Generated figures
├── results/
│   ├── all_results.json         # Complete session-level results
│   └── channel_ids.csv          # All channel selections
├── paper/
│   ├── paper.tex                # LaTeX source
│   ├── paper.pdf                # Compiled paper
│   └── references.bib           # Bibliography
├── requirements.txt             # Python dependencies
├── environment.yml              # Conda environment (alternative)
├── reproduce.sh                 # One-command reproduction script
└── README.md                    # This file
```

## Quick Start

### Option 1: Conda (Recommended)

```bash
# Clone repository
git clone https://github.com/[username]/schizago.git
cd schizago

# Create environment
conda env create -f environment.yml
conda activate schizago

# Run analysis (downloads data automatically, ~10-12 hours)
python analysis/run_analysis.py

# Generate figures
python figures/generate_figures.py
```

### Option 2: pip

```bash
git clone https://github.com/[username]/schizago.git
cd schizago

python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

pip install -r requirements.txt
python analysis/run_analysis.py
```

## Requirements

- Python 3.8+
- ~50 GB disk space (Allen Institute data cache)
- ~16 GB RAM recommended
- 10-12 hours compute time (5000 surrogates × 49 sessions × 3 segments)

## Data

Data is automatically downloaded from the Allen Institute Brain Observatory via `allensdk`. On first run, LFP data for qualifying sessions will be cached to `~/.allen/`.

**Manual download (optional):**
```python
from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache
cache = EcephysProjectCache.from_warehouse(manifest='/path/to/manifest.json')
```

## Reproducing Results

### Full Pipeline (10-12 hours)

```bash
./reproduce.sh
```

This will:
1. Install dependencies
2. Download required Allen Institute data
3. Run the complete analysis pipeline
4. Generate all figures
5. Produce `results/all_results.json`

### Verify Pre-computed Results

To verify our results without re-running the full analysis:

```bash
python analysis/verify_results.py
```

This compares the provided `results/all_results.json` against a quick sanity check on 2 sessions.

### Generate Figures Only

```bash
python figures/generate_figures.py --input results/all_results.json
```

## Analysis Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Frequency band | 30-80 Hz | Gamma range |
| Segment duration | 60 s | Balance temporal resolution and statistical power |
| Number of segments | 3 | Temporal replication |
| Surrogates | 5000 | Statistical rigor (floor p = 0.0002) |
| Filter | 4th-order Butterworth | Standard, zero-phase |
| Significance threshold | p < 0.05 | Per metric, combined via Fisher's method |
| Survival criteria | PLV p<0.05 AND wPLI p<0.05 AND ≥2/3 segments pass | Multi-criteria falsification |

## Results Summary

| Metric | Value |
|--------|-------|
| Sessions analyzed | 49 |
| Survivors | 29 (59.2%) |
| 95% CI | [45.2%, 71.8%] |
| Volume conduction caught | 4 sessions |
| No synchrony | 15 sessions |

See `results/all_results.json` for complete session-level data including channel IDs.

## Citation

```bibtex
@article{jocque2026hippocampal,
  title={Robust Hippocampal-Cortical Gamma Synchrony in Mouse Visual Processing: 
         A Falsification-Forward Analysis of the Allen Visual Coding Dataset},
  author={Jocque, Julian},
  journal={bioRxiv},
  year={2026},
  doi={10.1101/2026.XX.XX.XXXXXX}
}
```

## License

MIT License. See [LICENSE](LICENSE).

## Acknowledgments

- Allen Institute for Brain Science for the Visual Coding Neuropixels dataset
- This work was conducted independently without external funding

## Contact

Julian Jocque - julian@julianjocque.com
