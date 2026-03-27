# RTS-96 Cyber-Physical Power System — Reliability & Global Sensitivity Analysis

---

## Overview

This repository implements a full **reliability and global sensitivity analysis pipeline** for the IEEE RTS-96 24-bus power system modelled as a **cyber-physical system (CPPS)**. Every generator and transmission line has both a physical failure pathway and a cyber attack pathway, captured by a 6-state continuous-time Markov chain (CTMC). System-level reliability metrics — Loss of Load Probability (LOLP) and Expected Energy Not Served (EENS) — are computed via truncated N-1 contingency enumeration combined with an inner Monte Carlo simulation over wind and load uncertainty, solved at each step by a DC Optimal Power Flow (DC-OPF).

The sensitivity pipeline follows a two-step design to manage the extreme dimensionality:

1. **Morris elementary effects screening** (`morris_screening.py`) — cheaply identifies the ~50–80 parameters that actually influence reliability out of 563 total, discarding the rest
2. **Reduced-dimension Sobol analysis** (`run_reduced_sobol.py`) — runs statistically valid global sensitivity indices only on the parameters that survived Morris screening

---

## Table of Contents

- [System Description](#system-description)
- [Mathematical Background](#mathematical-background)
- [File Structure](#file-structure)
- [Dependencies](#dependencies)
- [Quick Start](#quick-start)
- [Script 1: full_component_sobol.py](#script-1-full_component_sobolpy)
- [Script 2: morris_screening.py](#script-2-morris_screeningpy)
- [Script 3: run_reduced_sobol.py](#script-3-run_reduced_sobolpy)
- [Complete Two-Step Workflow](#complete-two-step-workflow)
- [Output Files Reference](#output-files-reference)
- [ACCRE HPC Guide](#accre-hpc-guide)
- [Interpreting Results](#interpreting-results)
- [Parameter Reference](#parameter-reference)
- [Known Limitations](#known-limitations)

---

## System Description

The IEEE RTS-96 test system is a 24-bus, three-area reliability benchmark network widely used in power systems research.

| Property | Value |
|---|---|
| Buses | 24 |
| Generators | 32 (across 10 generator buses) |
| Transmission lines | 38 |
| Wind farms | 4 × 80 MW (buses 4, 5, 8, 9) |
| Load buses | 17 |
| Peak system load | 2,850 MW |
| Total installed generation | 3,405 MW (3,725 MW with wind) |
| Reserve margin | 19.5% (conventional only) |

**Generator types by capacity:**

| Capacity (MW) | FOR | MTTR (h) | Count |
|---|---|---|---|
| 12 | 2% | 60 | 5 |
| 20 | 10% | 50 | 4 |
| 50 | 1% | 20 | 6 |
| 76 | 2% | 40 | 4 |
| 100 | 4% | 50 | 3 |
| 155 | 4% | 40 | 4 |
| 197 | 5% | 50 | 3 |
| 350 | 8% | 100 | 1 |
| 400 | 12% | 150 | 2 |

---

## Mathematical Background

### 1. Rate Conversion from FOR and MTTR

Physical failure and repair rates for generators are derived from industry-standard reliability parameters:

```
μ = 1 / MTTR                         (repair rate, events/hour)
λ = (FOR × μ) / (1 − FOR)            (failure rate, events/hour)
p_up = μ / (λ + μ)                   (steady-state availability)
```

For transmission lines, rates are derived from outage frequency (per year) and repair duration:

```
λ = outage_rate_per_year / 8760       (failure rate, events/hour)
μ = 1 / repair_duration_hours         (repair rate, events/hour)
```

### 2. Six-State Cyber-Physical Markov Model

Each of the 70 components (32 generators + 38 lines) is modelled by a 6-state continuous-time Markov chain. The six states are:

| State | Label | Description |
|---|---|---|
| 0 | Operational | Fully in-service at rated capacity |
| 1 | Cyber-Manipulated | Adversary has altered controls; physically running |
| 2 | Cyber-Reconnaissance | Attacker probing; still in service but vulnerable |
| 3 | Partially Degraded | Reduced capacity; physical damage beginning |
| 4 | Physical Outage | Direct physical failure; completely offline |
| 5 | Cyber-Triggered Outage | Cyber event caused physical trip; offline |

**Service states** (component contributing to supply): {0, 1, 3, 5}  
**Outage states** (component offline): {4} for physical, {5} for cyber-triggered

The 6×6 rate matrix Q is built from eight transition parameters per component. Steady-state availability is solved analytically:

```
π Q = 0     subject to  Σ πᵢ = 1
p_up = π₀ + π₁ + π₂ + π₃
```

This solve is fully deterministic — no randomness — given the eight rate parameters.

### 3. Eight Rate Multipliers Per Component (563 Sobol Parameters)

Each component has eight uncertain transition rates. In the Sobol analysis, multipliers are applied to the nominal baseline rates:

```
actual_rate = nominal_rate × multiplier
```

| Rate Key | Bounds | Physical Meaning |
|---|---|---|
| `s_lambda_A` | [0.5, 2.0] | Cyber attack arrival rate |
| `s_theta` | [0.5, 1.5] | Attack routing probability (tighter bound: probability-valued) |
| `s_mu_M` | [0.5, 2.0] | Cyber monitoring/mitigation recovery rate |
| `s_mu_R` | [0.5, 2.0] | Cyber restoration rate (S2 → S0) |
| `s_mu_p` | [0.5, 2.0] | Cyber-triggered physical repair rate (S5 → S0) |
| `s_lambda_f` | [0.5, 2.0] | Physical fault rate (S0 → S4) |
| `s_lambda_p` | [0.5, 2.0] | Permanent failure rate (S3 → S4) |
| `s_mu` | [0.5, 2.0] | Physical repair rate (S4 → S0) |

**Total dimensions:**
- 32 generators × 8 = 256
- 38 lines × 8 = 304
- wind_k + wind_mean + load_sigma = 3
- **Total D = 563**

### 4. Contingency Enumeration (Layer 3)

All N-1 single-component outage states are enumerated analytically. For a given set of component availabilities {p_up}, the exact probability of each contingency k is:

```
P(contingency k) = ∏ (1 − p_up_i) for i ∈ outage_set
                 × ∏ p_up_j        for j ∉ outage_set
```

The enumerated contingencies (32 generator outages + 6–8 line outages + base case = 39–41 total) cover **>99.99%** of total system probability mass.

### 5. Inner Monte Carlo — Wind and Load Uncertainty (Layer 4)

For each contingency k, 50 independent scenarios of operational conditions are drawn:

**Wind speed** per farm: `v ~ Weibull(shape=wind_k, scale=λ_scale)` where `λ_scale = wind_mean / Γ(1 + 1/k)`

**Wind power curve:**
```
CF = 0                              if v < 3.0 m/s (cut-in) or v ≥ 25 m/s (cut-out)
CF = 1                              if v ≥ 12.0 m/s (rated)
CF = ((v − 3.0) / (12.0 − 3.0))³   if 3.0 ≤ v < 12.0 m/s
```

**Load scaling** (log-normal, mean = 1):
```
z ~ Normal(0, 1)
load_scale = exp(load_sigma × z − 0.5 × load_sigma²)
```

**DC-OPF** (HiGHS solver via `scipy.optimize.linprog`): minimize total load shed subject to power balance at every bus, generator capacity limits, and transmission line thermal limits.

**Aggregation over 50 scenarios:**
```
P_shed_k = count(shed > 0) / 50
E_shed_k = mean(shed_MW) over 50 scenarios
```

### 6. Reliability Metrics

```
LOLP = Σₖ  P(contingency k) × P_shed_k
EENS = Σₖ  P(contingency k) × E_shed_k × 8,760 hr/yr    [MWh/year]
```

### 7. Sobol Global Sensitivity Analysis

Sobol indices decompose output variance into contributions from individual inputs and their interactions.

**First-order index S1:**  
`S1_i = Var[E(Y | X_i)] / Var(Y)` — fraction of variance explained by X_i alone.

**Total-order index ST:**  
`ST_i = 1 − Var[E(Y | X_{~i})] / Var(Y)` — fraction including all interactions with other parameters.

**Saltelli estimator** (2010) requires `N × (D + 2)` model evaluations. Confidence interval scales as `CI ≈ 2 / √N`. Valid range for both S1 and ST is [0, 1].

### 8. Morris Elementary Effects

For each parameter i, the Elementary Effect at sample point x is:

```
EEᵢ = [ f(x₁, …, xᵢ + Δ, …, x_D) − f(x₁, …, x_D) ] / Δ
```

Repeat across r trajectories. Summarise as:

```
μ*ᵢ = (1/r) Σ |EEᵢ|    (mean absolute effect — how influential?)
σᵢ  = std(EEᵢ)          (standard deviation — does it interact?)
```

Cost: `r × (D + 1)` evaluations. At r=20, D=563: **11,280 evaluations** vs 36,160 for Sobol at N=64.

---

## File Structure

```
project/
│
├── full_component_sobol.py     # Core model + full 563D Sobol driver
├── morris_screening.py         # Step 1: Morris screening to find active parameters
├── run_reduced_sobol.py        # Step 2: Sobol on reduced active parameter set
│
├── README.md                   # This file
│
├── morris_output/              # Created by morris_screening.py
│   ├── morris_results.csv
│   ├── morris_active_params.json
│   ├── morris_mu_star_lolp.png
│   ├── morris_mu_star_eens.png
│   └── morris_top_bar.png
│
├── reduced_sobol_output/       # Created by run_reduced_sobol.py
│   └── reduced_sobol_indices.csv
│
└── sobol_results_full_component.csv     # Created by full_component_sobol.py
    component_availability_full_component.csv
    sobol_indices_summary_full_component.csv
    nominal_component_rates_full_component.csv
```

---

## Dependencies

```bash
pip install numpy pandas scipy matplotlib SALib
```

| Package | Version | Purpose |
|---|---|---|
| numpy | ≥ 1.24 | Numerical arrays |
| pandas | ≥ 2.0 | Data storage and CSV output |
| scipy | ≥ 1.11 | linprog (HiGHS), Weibull distribution, special functions |
| matplotlib | ≥ 3.7 | Morris diagnostic plots |
| SALib | ≥ 1.5 | Saltelli Sobol sampler/analyser, Morris sampler/analyser |

**Python version:** 3.9 or later recommended.

**On ACCRE (Vanderbilt HPC):**
```bash
module load Python/3.10
pip install --user numpy pandas scipy matplotlib SALib
```

---

## Quick Start

### Fastest possible test (verify the pipeline works — ~5 minutes):
```bash
# Full model, 2 Sobol samples, 6 lines, 20 inner trials
python full_component_sobol.py --n-samples 2 --inner-trials 20 --enum-top-lines 6

# Morris quick test, 5 trajectories
python morris_screening.py --trajectories 5 --inner-trials 20 --enum-top-lines 6

# Reduced Sobol quick test
python run_reduced_sobol.py \
    --active-params morris_output/morris_active_params.json \
    --n-samples 4 --inner-trials 20
```

### Recommended research run (on ACCRE, 20 cores):
```bash
# Step 1: Morris (~2–3 hours)
python morris_screening.py \
    --trajectories 20 --inner-trials 50 \
    --enum-top-lines 8 --n-jobs 20

# Step 2: Reduced Sobol (~6–8 hours)
python run_reduced_sobol.py \
    --active-params morris_output/morris_active_params.json \
    --n-samples 512 --inner-trials 50 \
    --enum-top-lines 8 --n-jobs 20
```

---

## Script 1: `full_component_sobol.py`

### Purpose

This is the **core model file**. It defines the entire RTS-96 cyber-physical system and runs the full 563-dimensional Sobol sensitivity analysis from scratch. It is both a standalone runner and the importable module that `morris_screening.py` and `run_reduced_sobol.py` depend on.

### What It Contains

| Section | Description |
|---|---|
| `rates_from_FOR_MTTR()` | Converts FOR + MTTR to λ/μ rates |
| `rates_from_lambdaP_duration()` | Converts line outage rate + duration to λ/μ |
| `weibull_scale_from_mean()` | Computes Weibull scale parameter from mean wind speed and shape k |
| `wind_power_cf()` | Cubic ramp power curve: cut-in 3 m/s, rated 12 m/s, cut-out 25 m/s |
| `lognormal_multiplier_from_z()` | Log-normal load scaling with mean = 1 |
| `stable_seed()` | Deterministic seed generation using SHA-256 (reproducible across platforms) |
| `build_system()` | Constructs all generators, lines, wind farms, and loads from hard-coded RTS-96 data |
| `build_q_matrix()` | Builds 6×6 transition rate matrix for one component |
| `component_steady_state()` | Solves πQ=0, returns p_up |
| `apply_component_multipliers()` | Scales nominal rates by Sobol-sampled multipliers |
| `DCOPF` | DC-OPF linear program class using HiGHS via scipy.optimize.linprog |
| `enumerate_truncated_contingencies()` | Exact N-1 contingency probability enumeration |
| `inner_monte_carlo_conditional()` | 50-scenario wind/load Monte Carlo per contingency |
| `compute_reliability()` | Full nested model: Markov → contingency → MC → OPF → LOLP/EENS |
| `build_sobol_problem()` | Builds the 563-dimensional SALib problem dict |
| `run_sobol()` | Saltelli sampling + parallel evaluation + Sobol index computation |

### Command-Line Arguments

```
python full_component_sobol.py [OPTIONS]
```

| Argument | Default | Description |
|---|---|---|
| `--n-samples` | 2 | Base Sobol sample size N. Total evaluations = N × (D+2) = N × 565 |
| `--inner-trials` | 50 | Monte Carlo scenarios per contingency in Layer 4 |
| `--max-outages` | 1 | Maximum contingency order (1 = N-1 only, 2 = N-2) |
| `--enum-top-lines` | 8 | Number of highest-risk transmission lines explicitly enumerated |
| `--seed` | 123 | Base random seed for reproducibility |

### Evaluation Cost Table

| N | Total Evaluations | CI(ST) | Status |
|---|---|---|---|
| 2 | 1,130 | ±1.41 | Code test only — indices unreliable |
| 64 | 36,160 | ±0.25 | Minimum viable for rough rankings |
| 128 | 72,320 | ±0.18 | Directional rankings |
| 512 | 289,280 | ±0.09 | Reliable — suitable for publication |
| 1,024 | 578,560 | ±0.06 | Publication grade |

### Example Commands

```bash
# Minimum viable run (rough rankings, ~5 hrs on 20 ACCRE cores)
python full_component_sobol.py --n-samples 64 --inner-trials 50 --enum-top-lines 8

# Production quality
python full_component_sobol.py --n-samples 512 --inner-trials 50 --enum-top-lines 8

# N-2 contingency analysis (expensive — adds line×line pairs)
python full_component_sobol.py --n-samples 64 --max-outages 2 --enum-top-lines 6
```

### Output Files

**`sobol_results_full_component.csv`**  
One row per Sobol sample evaluation (N × 565 rows total). Columns include all 563 input multipliers, derived `wind_scale`, and output columns: `LOLP`, `EENS_MWh_per_year`, `enum_mass`, `tail_mass`, `n_contingencies`, `mean_p_up_gen`, `mean_p_up_line`.

**`component_availability_full_component.csv`**  
One row per Sobol sample. Columns: `sample_id`, then `p_up_G1` ... `p_up_G32`, `p_up_L1` ... `p_up_L38`. Shows how each component's availability varies across the Sobol samples.

**`sobol_indices_summary_full_component.csv`**  
One row per parameter (563 rows). Columns: `parameter`, `S1_LOLP`, `S1_conf_LOLP`, `ST_LOLP`, `ST_conf_LOLP`, `S1_EENS`, `S1_conf_EENS`, `ST_EENS`, `ST_conf_EENS`.

**`nominal_component_rates_full_component.csv`**  
One row per component (70 rows). Columns: `component_id`, `component_type`, and all eight nominal transition rates before any multiplier is applied.

### How `compute_reliability()` Works — Step by Step

```
Input: 563-element params dict (all multipliers + wind_k + wind_mean + load_sigma)

1. build_system() — load RTS-96 data (fixed, same every call)

2. For each of 32 generators:
   a. apply_component_multipliers() — scale nominal rates by this sample's multipliers
   b. build_q_matrix() — 6×6 transition rate matrix
   c. component_steady_state() — solve πQ=0, get p_up_gen

3. For each of 38 lines:
   a. Same three steps → p_up_line

4. select_lines_for_enumeration() — pick top-N lines by risk score λ/limit

5. enumerate_truncated_contingencies() — N-1 exact probabilities
   → 39–41 (outage_set, probability) pairs covering 99.99% of probability space

6. For each contingency:
   a. Remove outaged generators from capacity
   b. Remove outaged lines from network
   c. inner_monte_carlo_conditional() — 50 wind/load scenarios × DC-OPF
   → E_shed_k, P_shed_k

7. LOLP = Σ P(k) × P_shed_k
   EENS = Σ P(k) × E_shed_k × 8760

Output: LOLP, EENS, enum_mass, tail_mass, n_contingencies,
        mean_p_up_gen, mean_p_up_line, p_up_gen_dict, p_up_line_dict
```

---

## Script 2: `morris_screening.py`

### Purpose

Morris screening is the **first step** in the efficient two-step sensitivity workflow. It identifies which of the 563 parameters actually influence LOLP and EENS, so the subsequent Sobol analysis can run on a much smaller parameter space.

### Why Morris Before Sobol?

Running Sobol directly on 563 parameters at N=512 requires **289,280 model evaluations** (~4–5 days on 20 ACCRE cores). In practice, most parameters do nothing — the failure rates of generators on lightly loaded buses have essentially zero effect on system-level reliability. Morris finds these unimportant parameters cheaply and eliminates them.

Morris cost at r=20: **r × (D+1) = 20 × 564 = 11,280 evaluations** (~2–3 hours). This is the screening investment that enables everything downstream.

### The Elementary Effect — What It Measures

For each parameter i, Morris draws a random starting point x, then nudges x_i by one grid step Δ and re-evaluates the model:

```
EEᵢ = [ f(x₁, …, xᵢ + Δ, …, x_D) − f(x₁, …, xᵢ, …, x_D) ] / Δ
```

This is done from r=20 independent random starting points (trajectories). The summary statistics are:

- **μ\* = mean|EEᵢ|** — overall influence. High μ\* means the parameter reliably changes the output.
- **σ = std(EEᵢ)** — non-linearity and interactions. High σ with high μ\* means the parameter interacts strongly with others.

### Interpreting the μ\* vs σ Plot

```
         High σ
            │
            │   ●  Interactive &        ← keep these (both μ* and σ high)
            │      important
            │
            │           ●  Also keep   ← important but linear
            │
  ──────────┼──────────────────────── μ*
            │
  Low σ     │   ● ● ● ●  Negligible    ← fix these at nominal (1.0)
            │
```

Parameters clustered near zero on the μ\* axis can be safely fixed. Parameters with any significant μ\* — even if σ is low — should be kept for Sobol.

### Active Parameter Threshold

After computing μ\* for all 563 parameters, a parameter is classified as **active** (kept for Sobol) if its μ\* exceeds the threshold percentile in either LOLP or EENS:

```
active if:  mu_star_LOLP ≥ percentile(mu_star_LOLP, threshold_pct)
         OR mu_star_EENS ≥ percentile(mu_star_EENS, threshold_pct)
```

The three wind/load parameters (`wind_k`, `wind_mean`, `load_sigma`) are always kept regardless.

**Tuning the threshold:**
- `--threshold-percentile 75` (default) — keeps top 25% most influential parameters
- `--threshold-percentile 80` — stricter, smaller active set, faster Sobol
- `--threshold-percentile 60` — conservative, larger active set, slower Sobol

After running, check `morris_results.csv` and plot the μ\* distribution. If you see a clear elbow (a cluster near zero then a sharp jump), tighten the threshold. If it is gradual, keep 75.

### Command-Line Arguments

```
python morris_screening.py [OPTIONS]
```

| Argument | Default | Description |
|---|---|---|
| `--trajectories` | 20 | Number of Morris trajectories r. Cost = r × (D+1) = r × 564 |
| `--inner-trials` | 50 | Inner MC scenarios per contingency |
| `--max-outages` | 1 | Maximum contingency order |
| `--enum-top-lines` | 8 | Lines explicitly enumerated in contingency analysis |
| `--n-jobs` | 1 | Parallel workers (set to CPU count on ACCRE) |
| `--seed` | 42 | Base random seed |
| `--num-levels` | 4 | Number of levels in Morris grid (4 is standard) |
| `--threshold-percentile` | 75.0 | μ\* percentile cutoff for active/inactive classification |
| `--output-dir` | morris_output | Directory for all output files |

### Example Commands

```bash
# Quick test — verify it runs (5 minutes)
python morris_screening.py --trajectories 5 --inner-trials 20 --enum-top-lines 6

# Standard screening run (recommended)
python morris_screening.py --trajectories 20 --inner-trials 50 --enum-top-lines 8

# On ACCRE with 20 cores
python morris_screening.py \
    --trajectories 20 --inner-trials 50 \
    --enum-top-lines 8 --n-jobs 20 \
    --output-dir morris_output

# Conservative threshold — keep more parameters
python morris_screening.py \
    --trajectories 20 --inner-trials 50 \
    --threshold-percentile 60 \
    --n-jobs 20
```

### Trajectory Count Guide

| r | Evaluations | Wall time (20 cores) | Screening quality |
|---|---|---|---|
| 5 | 2,820 | ~30 min | Exploratory only |
| 10 | 5,640 | ~60 min | Rough screening |
| **20** | **11,280** | **~2–3 hrs** | **Recommended** |
| 50 | 28,200 | ~5–6 hrs | High confidence screening |

### Output Files

**`morris_output/morris_results.csv`**  
One row per parameter (563 rows), sorted by μ\* LOLP descending. Columns: `parameter`, `mu_star_LOLP`, `sigma_LOLP`, `rank_LOLP`, `mu_star_EENS`, `sigma_EENS`, `rank_EENS`.

**`morris_output/morris_active_params.json`**  
JSON file consumed by `run_reduced_sobol.py`. Contains:
```json
{
  "n_active": 58,
  "n_inactive": 505,
  "threshold_pct_lolp": 0.000123,
  "threshold_pct_eens": 0.412,
  "active_params": ["s_lambda_f_G13", "s_lambda_f_L37", ...],
  "inactive_params": ["s_lambda_A_G5", "s_theta_G8", ...]
}
```

**`morris_output/morris_mu_star_lolp.png`**  
Scatter plot of μ\* vs σ for all 563 parameters, coloured by type (generator physical, generator cyber, line physical, line cyber, wind/load). Top-15 parameters are labelled. Threshold line shown.

**`morris_output/morris_mu_star_eens.png`**  
Same as above but for EENS as the output metric.

**`morris_output/morris_top_bar.png`**  
Horizontal bar chart showing top-40 parameters side-by-side for LOLP (left) and EENS (right), coloured by parameter type.

### What the Colour Coding Means in the Plots

| Colour | Parameter Type |
|---|---|
| Red | Generator physical rates (λ_f, λ_p, μ) |
| Orange | Generator cyber rates (λ_A, θ, μ_M, μ_R, μ_p) |
| Blue | Line physical rates |
| Teal | Line cyber rates |
| Purple | Wind/load (wind_k, wind_mean, load_sigma) |

If the red/orange cluster dominates the upper-right corner of the μ\*–σ plot, generators drive reliability risk. If blue/teal dominates, transmission corridor vulnerability is the primary concern.

### Importable Functions

`morris_screening.py` also exports two utility functions used by `run_reduced_sobol.py`:

```python
from morris_screening import build_reduced_sobol_problem, expand_params_with_inactive

# Build SALib problem dict with only active parameters
problem = build_reduced_sobol_problem("morris_output/morris_active_params.json")

# Fill inactive params back at 1.0 so compute_reliability() receives all 563
full_params = expand_params_with_inactive(reduced_row_dict, "morris_output/morris_active_params.json")
```

---

## Script 3: `run_reduced_sobol.py`

### Purpose

This is the **second step** of the two-step workflow. It reads the active parameter list produced by Morris screening and runs the full Sobol global sensitivity analysis on **only those parameters**, holding all inactive parameters fixed at their nominal multiplier (1.0). Because D is now ~50–80 instead of 563, Sobol at N=512 requires ~26,000–43,000 evaluations instead of 289,280 — roughly **7–11× fewer**.

### How Inactive Parameters Are Handled

Every inactive parameter is set to its nominal multiplier of **1.0** — meaning its corresponding transition rate is used exactly as calibrated from the FOR/MTTR data, with no uncertainty. This is the correct choice because Morris confirmed these parameters do not influence the output; their uncertainty is irrelevant to reliability.

The `expand_params_with_inactive()` function in `morris_screening.py` handles this automatically:

```python
# In the evaluation loop:
reduced_params = {name: value for name, value in zip(active_names, row)}
full_params = expand_params_with_inactive(reduced_params, active_params_json)
# full_params now has all 563 keys — active ones from the sample, inactive ones at 1.0
compute_reliability(full_params, ...)
```

### Command-Line Arguments

```
python run_reduced_sobol.py [OPTIONS]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `--active-params` | Yes | — | Path to `morris_active_params.json` from Step 1 |
| `--n-samples` | No | 512 | Sobol base sample count N |
| `--inner-trials` | No | 50 | Inner MC scenarios per contingency |
| `--max-outages` | No | 1 | Max contingency order |
| `--enum-top-lines` | No | 8 | Lines to enumerate |
| `--n-jobs` | No | 1 | Parallel workers |
| `--seed` | No | 123 | Base random seed |
| `--output-dir` | No | reduced_sobol_output | Output directory |

### Evaluation Cost vs Quality

Assumes D_reduced ≈ 60 after Morris screening at 75th percentile threshold.

| N | Evaluations (D=60) | CI(ST) | Wall time (20 cores) |
|---|---|---|---|
| 32 | 1,984 | ±0.35 | ~20 min |
| 128 | 7,936 | ±0.18 | ~1.5 hrs |
| **512** | **31,744** | **±0.09** | **~6–8 hrs** |
| 1,024 | 63,488 | ±0.06 | ~12–14 hrs |

Compare to full 563D Sobol at same CI target:
- CI ±0.09 requires N=512 → 289,280 evaluations → **4–5 days**
- Same CI on reduced D=60 requires N=512 → 31,744 evaluations → **6–8 hours**

### Example Commands

```bash
# Minimum viable after Morris (rough rankings)
python run_reduced_sobol.py \
    --active-params morris_output/morris_active_params.json \
    --n-samples 128 --inner-trials 50

# Recommended — CI < 0.10
python run_reduced_sobol.py \
    --active-params morris_output/morris_active_params.json \
    --n-samples 512 --inner-trials 50 \
    --enum-top-lines 8 --n-jobs 20

# Publication quality
python run_reduced_sobol.py \
    --active-params morris_output/morris_active_params.json \
    --n-samples 1024 --inner-trials 50 \
    --enum-top-lines 8 --n-jobs 20
```

### Output Files

**`reduced_sobol_output/reduced_sobol_indices.csv`**  
One row per active parameter, sorted by ST_LOLP descending. Columns: `parameter`, `S1_LOLP`, `S1_conf_LOLP`, `ST_LOLP`, `ST_conf_LOLP`, `S1_EENS`, `S1_conf_EENS`, `ST_EENS`, `ST_conf_EENS`.

The script also prints mean CI for both metrics and indicates whether CI < 0.10 has been achieved.

---

## Complete Two-Step Workflow

### Step-by-Step on ACCRE

**1. Set up the environment:**
```bash
module load Python/3.10
pip install --user numpy pandas scipy matplotlib SALib
```

**2. Upload your files to ACCRE:**
```
full_component_sobol.py
morris_screening.py
run_reduced_sobol.py
```

**3. Create SLURM script for Morris screening (`morris_job.slurm`):**
```bash
#!/bin/bash
#SBATCH --job-name=morris_screen
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=morris_%j.log

module load Python/3.10
cd /path/to/your/scripts

python morris_screening.py \
    --trajectories 20 \
    --inner-trials 50 \
    --enum-top-lines 8 \
    --n-jobs 20 \
    --seed 42 \
    --output-dir morris_output
```

**4. Submit Morris job:**
```bash
sbatch morris_job.slurm
# Check progress:
tail -f morris_<jobid>.log
```

**5. After Morris completes, check how many parameters survived:**
```bash
python -c "
import json
with open('morris_output/morris_active_params.json') as f:
    d = json.load(f)
print(f'Active: {d[\"n_active\"]}  Inactive: {d[\"n_inactive\"]}')
"
```

**6. Create SLURM script for reduced Sobol (`reduced_sobol_job.slurm`):**
```bash
#!/bin/bash
#SBATCH --job-name=reduced_sobol
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=reduced_sobol_%j.log

module load Python/3.10
cd /path/to/your/scripts

python run_reduced_sobol.py \
    --active-params morris_output/morris_active_params.json \
    --n-samples 512 \
    --inner-trials 50 \
    --enum-top-lines 8 \
    --n-jobs 20 \
    --seed 123 \
    --output-dir reduced_sobol_output
```

**7. Submit reduced Sobol job:**
```bash
sbatch reduced_sobol_job.slurm
```

**Total wall time: ~10–12 hours across both jobs for publication-quality results.**

---

## Output Files Reference

| File | Script | Description |
|---|---|---|
| `sobol_results_full_component.csv` | full_component_sobol.py | All inputs + outputs for every Sobol sample |
| `component_availability_full_component.csv` | full_component_sobol.py | Per-sample p_up for all 70 components |
| `sobol_indices_summary_full_component.csv` | full_component_sobol.py | S1/ST/CI for all 563 parameters |
| `nominal_component_rates_full_component.csv` | full_component_sobol.py | Baseline transition rates for all components |
| `morris_output/morris_results.csv` | morris_screening.py | μ*, σ, rank for all 563 parameters |
| `morris_output/morris_active_params.json` | morris_screening.py | Active/inactive parameter lists |
| `morris_output/morris_mu_star_lolp.png` | morris_screening.py | μ* vs σ scatter plot (LOLP) |
| `morris_output/morris_mu_star_eens.png` | morris_screening.py | μ* vs σ scatter plot (EENS) |
| `morris_output/morris_top_bar.png` | morris_screening.py | Top-40 parameters bar chart |
| `reduced_sobol_output/reduced_sobol_indices.csv` | run_reduced_sobol.py | S1/ST/CI on active parameters only |

---

## Interpreting Results

### Reading Sobol Indices

**A reliable result** (N ≥ 64 on reduced problem):
- ST values between 0 and 1 for most parameters
- CI(ST) < 0.25 for rough rankings; < 0.10 for publication
- Σ S1 ≤ 1.0 (if > 1.0, significant under-sampling)

**Typical findings in this system:**

| Finding | Interpretation |
|---|---|
| `s_lambda_f_G13` ST is highest for EENS | G13's physical failure rate drives energy not served — invest in physical hardening of Bus 13 units |
| `s_lambda_f_L37` top for both LOLP and EENS | L37 (Bus 20→23) is the most critical transmission corridor — structural redundancy needed |
| `s_mu_R_L14` top for LOLP | L14's cyber restoration speed drives outage duration — priority incident response for this line |
| Generator cyber rates rank higher than physical rates | Cyber attack pathways are as important as physical failure modes for this system |

### CI and Sample Size Guidance

```
CI > 0.50  →  N too small. Rankings are noise. Do not interpret.
CI 0.25–0.50 →  Very rough — only top-3 to top-5 can be trusted.
CI 0.10–0.25 →  Rough rankings. Top-10 can be trusted directionally.
CI < 0.10  →  Reliable. Suitable for publication with caveats.
CI < 0.06  →  High confidence. All ranked parameters are trustworthy.
```

### The μ\* vs σ Plot (Morris)

Look for a **natural gap** in the μ\* distribution. Parameters clearly to the right of the gap are active; everything in the cluster near zero is negligible. If no gap is visible, the threshold percentile setting will determine the cut.

High σ relative to μ\* indicates strong non-linear interactions. These parameters are particularly important for Sobol — the ST >> S1 gap will be large for them.

---

## Parameter Reference

### Wind and Load Parameters

| Parameter | Range | Description |
|---|---|---|
| `wind_k` | [1.5, 3.0] | Weibull shape. Low k = erratic wind; high k = steady near mean |
| `wind_mean` | [6.0, 10.0] m/s | Mean wind speed across all four farms |
| `load_sigma` | [0.03, 0.12] | Log-normal std dev for load fluctuation |

### Generator Rate Keys (per generator G1–G32)

| Key | Meaning | Note |
|---|---|---|
| `s_lambda_A_Gx` | Cyber attack rate | High = frequent attack attempts |
| `s_theta_Gx` | Attack routing probability | Fraction taking destructive path |
| `s_mu_M_Gx` | Monitoring recovery rate | Speed of cyber manipulation recovery |
| `s_mu_R_Gx` | Restoration rate | Speed of cyber reconnaissance recovery |
| `s_mu_p_Gx` | Cyber-triggered repair rate | Speed of recovery from cyber-induced trip |
| `s_lambda_f_Gx` | Physical fault rate | Direct physical failure frequency |
| `s_lambda_p_Gx` | Permanent failure rate | Failure from degraded state |
| `s_mu_Gx` | Physical repair rate | Speed of physical outage repair |

### Line Rate Keys (per line L1–L38)

Same eight keys as generators, with `Lx` suffix. Key lines to watch:

| Line | Route | Capacity | Critical Parameter | Why |
|---|---|---|---|---|
| L14 | Bus 9→11 | 400 MW | `s_mu_R_L14` | MTTR = 768 hrs (32 days) makes cyber recovery speed critical |
| L37 | Bus 20→23 | 500 MW | `s_lambda_f_L37` | Only export corridor from Bus 23 generation pocket (505 MW) |
| L28 | Bus 16→17 | 500 MW | `s_lambda_p_L28` | Backbone corridor for northern load zone |

---

## Known Limitations

| Limitation | Impact | Future Fix |
|---|---|---|
| N=2 pilot Sobol has CI = ±1.41 — indices are invalid as absolute numbers | Rankings are directional only; specific ST values unreliable | Run N ≥ 64 using the two-step Morris → reduced Sobol workflow |
| N-1 contingencies only (max_outages=1) | N-2 compound failures (e.g. G13 + L37 simultaneous) not captured | Set `--max-outages 2` — significant cost increase |
| DC-OPF ignores reactive power and voltage constraints | Load shed slightly underestimated vs full AC-OPF | Replace `linprog` with AC-OPF (PyPSA or pandapower) |
| 50 inner Monte Carlo scenarios per contingency gives ±14% EENS error | Inner integration is the largest accuracy bottleneck | Replace with Smolyak sparse grid (Level-2, d=3 → 65 deterministic points, ~5× lower error) |
| All wind farms share the same wind_k and wind_mean | Spatial diversity in wind resources not captured | Extend to per-farm wind parameters (+6 dimensions) |
| Cyber parameters are assumed independent across components | In reality, coordinated attacks affect multiple components simultaneously | Introduce correlated attack scenarios |

---

## Citation and Acknowledgements

This work extends the cyber-physical reliability framework from:

> Rostami et al. (2023). "Cyber-Physical Security Assessment for Electric Power Systems." *IEEE Transactions on Smart Grid*.

Sensitivity analysis methodology follows:

> Saltelli, A. et al. (2010). "Variance based sensitivity analysis of model output. Design and estimator for the total sensitivity index." *Computer Physics Communications*, 181(2), 259–270.

> Morris, M.D. (1991). "Factorial sampling plans for preliminary computational experiments." *Technometrics*, 33(2), 161–174.

RTS-96 test system data from:

> IEEE Task Force on Test Systems for Economic Analysis (1996). "The IEEE Reliability Test System-1996."

## Reference for Sobol and Morris:

> https://salib.readthedocs.io/en/latest/

> https://gsa-module.readthedocs.io/en/stable/implementation/morris_screening_method.html

*Vanderbilt University · Power Systems Cybersecurity Research · March 2026*
