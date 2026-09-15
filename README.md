# Longitudinal subclonal concordance of paired cfDNA and gDNA

A Bayesian hierarchical Dirichlet process (HDP) mixture model for tracking subclonal
structure in paired cell-free DNA (cfDNA) and genomic DNA (gDNA) 5-hydroxymethylcytosine
(5-hmC) profiles collected longitudinally from multiple myeloma patients on a KRd
(carfilzomib/lenalidomide/dexamethasone) trial.

## The question

Bone marrow biopsy (gDNA) is the reference for assessing myeloma subclonal composition, but
it is invasive and cannot be repeated often. A peripheral blood draw (cfDNA) can be. The
question this model addresses is whether cfDNA **tracks** the marrow: if the subclonal
composition inferred from blood moves with the composition inferred from marrow in the same
patient over time, cfDNA becomes a candidate surrogate for serial monitoring.

## The model

Each observation is one (patient, timepoint, source) triple; features are gene-body 5-hmC
counts, VST-transformed jointly across sources.

Let `i` index patients, `j` genes, `t` timepoints, and `k` subclonal components.

**gDNA** is a finite mixture over subclonal signatures:

```
Y^g_{ijt}  ~  sum_{k=1..K} omega^g_{ikt} * N( theta_{ijk}, sigma_k^2 )
```

**cfDNA** is the same mixture, diluted by a non-tumour background:

```
Y^cf_{ijt} ~ (1 - omega_{0it}) * sum_{k=1..K} omega^cf_{ikt} * N( theta_{ijk}, sigma_k^2 )
           +      omega_{0it}  *                              N( theta_{0ij}, sigma_0^2 )
```

Key structural choices:

- **Signatures `theta_{ijk}` are shared across sources** and subject-indexed, with
  hierarchical shrinkage `theta_{ijk} ~ N(thetabar_{jk}, varsigma_k^2)` and a horseshoe
  prior on the population signature `thetabar_{jk}`. Concordance therefore lives entirely
  in the mixture weights, not in a cross-source signature correlation.
- **Contamination `omega_{0it}` is time-varying** (per subject-timepoint, `Beta(a_0, b_0)`),
  so it is a trajectory rather than one number per patient.
- **The cfDNA background has an estimated mean** `theta_{0ij}`, shrunk toward a cohort
  profile. This matters: VST-transformed data are not mean-centred, so a background pinned
  at zero can only fit noise.

Two models are compared: `M1` (tracking — a shared global composition `pi` anchors both
sources) against `M0` (non-tracking — independent source-specific compositions).

Inference is by collapsed/Rao-Blackwellised Gibbs sampling with Pólya-urn predictive
allocation, optional simulated annealing and parallel tempering, with the allocation and
signature updates in Rcpp.

## ⚠️ The real data is NOT part of this repository

**No patient data is included here, and none will be.** The KRd trial 5-hmC counts and
sample keys are protected human-subjects data and are not redistributable. Specifically,
these files are absent and must be supplied by an authorised user:

```
KRd trial/5hmC data/kRd-cfDNA_genebody_count.RDS
KRd trial/5hmC data/kRd-gDNA_genebody_count.RDS
KRd trial/5hmC data/kRd-cfDNA_sample_key.csv
KRd trial/5hmC data/kRd-gDNA_sample_key.csv
```

Everything in `results/` is likewise excluded, since fitted objects are derived from patient
data.

**What this means in practice:** `R/02_preprocess.R` and every script downstream of it will
fail without those inputs. The parts of the repository that run standalone are the
**simulation study** and the **sampler validation suite** — see "Running without the real
data" below. Those exercise the full model and sampler on synthetic data generated from the
model's own DGP, which is sufficient to reproduce the methodological claims, verify the
sampler, and check the C++/R implementations against each other. Reproducing the *applied*
results requires data access through the study investigators.

## Repository layout

```
R/
  00_bootstrap.R          package loading + Rcpp compilation (source this first)
  01_lib_core.R           sampler, R fallbacks, relabelling, diagnostics, WAIC/LOO
  02_preprocess.R         raw counts -> DE gene selection -> results/model_data.rds
  04_*..31_*              downstream analysis, figures, tables, manuscript outputs
  cpp/sampler_core.cpp    Rcpp core: allocations, signature update, log-likelihood
  simulation/             data-generating process, scenario configs, validation drivers
  archive/                superseded scripts, kept for provenance
slurm/                    UChicago RCC (Midway3) batch scripts + their R entry points
```

`R/02_preprocess.R` is self-contained: it performs the differentially-expressed gene
selection itself (pooled limma-voom, `~ timepoint + source` with `duplicateCorrelation` on
patient, adj. *P* < 0.05 and |logFC| > 0.5, giving ~310 genes) and writes both the gene list
and the model input.

## Requirements

R (≥ 4.2; developed against 4.4.2) with:

- **Core sampler:** `Rcpp`, `coda`, `posterior`, `loo`
- **Preprocessing:** `DESeq2`, `limma` (Bioconductor), `dplyr`, `tidyr`, `readxl`
- **Figures/tables:** `ggplot2`, `ggrepel`, `patchwork`, `gridExtra`, `grid`, `scales`,
  `xtable`, `pdftools`
- **Downstream analyses:** `lme4`, `survival`

A C++ compiler is required. The Rcpp path is the only supported one for production runs; the
pure-R fallbacks exist for verification and are far too slow for real fits.

## Running without the real data

This is the recommended entry point for anyone without data access.

```bash
# 1. Confirm the toolchain and that Rcpp compiles
Rscript -e 'source("R/00_bootstrap.R"); cat("rcpp:", .rcpp_available, "\n")'

# 2. Sampler validation suite (12 tests, synthetic data, a few minutes)
Rscript R/simulation/run_local_validation.R

# 3. A quick end-to-end smoke fit
Rscript R/simulation/run_smoke_test.R
```

The validation suite covers parameter recovery, contamination identifiability, parallel
tempering, and C++/R agreement. On the last point, note what equivalence means here:
`update_theta` and `compute_loglik` are **stream-identical** between the two
implementations and are checked numerically (T11), whereas `compute_allocations` **cannot**
be — R's `sample.int` consumes a different number of uniforms than the C++ Fisher-Yates
shuffle, and R uses Walker alias sampling where the C++ uses an inverse-CDF search. It is
therefore checked **distributionally** (T12). A draw-for-draw mismatch in the allocations is
expected and is not a bug.

## Running the full pipeline (requires data access)

From the project root, with the `KRd trial/5hmC data/` directory in place:

```bash
Rscript R/02_preprocess.R          # -> results/model_data.rds, candidate_genes_union.{rds,csv}
```

Then on the cluster:

```bash
mkdir -p logs/krd
sbatch slurm/krd_M1.sbatch         # fit M1 (tracking), 4 chains
sbatch slurm/krd_M0.sbatch         # fit M0 (non-tracking), 4 chains
sbatch slurm/krd_combine.sbatch    # after BOTH complete: diagnostics + model comparison
```

Cross-validation and sensitivity analyses:

```bash
sbatch slurm/07_kfold_cv.sbatch  && sbatch slurm/07b_kfold_aggregate.sbatch
sbatch slurm/09_lopo_cv.sbatch   && sbatch slurm/09b_lopo_aggregate.sbatch
sbatch slurm/06_kappa_sensitivity.sbatch && sbatch slurm/06b_kappa_aggregate.sbatch
```

Simulation study (array jobs; `slurm/submit_all.sh` drives the whole set):

```bash
sbatch slurm/01_sim_scenario1.sbatch   # ... through 05_sim_scenario5.sbatch
sbatch slurm/05_sim_aggregate_results.sbatch
```

Finally, to regenerate manuscript figures and tables locally:

```bash
Rscript R/30_generate_all_manuscript_outputs.R
```

## HPC notes (UChicago RCC / Midway3)

The batch scripts target Midway3 and submit to Yuan Ji's allocation:

```
#SBATCH --account=pi-jiyuan
#SBATCH --partition=caslake
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0
```

Adjust `--account` if you are submitting under a different allocation, and
`--mail-user` to your own address.

Two things to watch:

- **Walltime.** Seven scripts request `--time=36:00:00`, which sits exactly at the caslake
  maximum. If SLURM rejects a job as over-limit, drop to `35:59:00` or split the run. Each
  such script carries an inline note.
- **Array throttling.** The simulation arrays are throttled with `%N` to stay under the
  per-user queued-job cap. Raise it only if you know your QOS allows it.

Install the R dependencies into your user library once before the first batch submission;
the scripts do not install packages themselves.

## Reproducibility status

Two caveats worth stating plainly, both documented at greater length in the manuscript:

1. **The applied contamination and model-comparison numbers are being recomputed.** The
   model specification changed (shared subject-indexed signatures, time-varying
   contamination, estimated background mean), and a defect was found and fixed in the
   collapsed cfDNA allocation step — the background had been compared against
   *un-normalised* Pólya-urn masses, which handed every tumour component a spurious
   `+log(p)` bonus and drove the contamination fraction toward zero. On synthetic data with
   a known background, the corrected sampler recovers `omega_0 = 0.232` against a truth of
   `0.251` (correlation 0.96); before the fix it returned `0.023`. Earlier reports of
   "negligible contamination" are retracted and should not be cited.
2. **The leave-one-patient-out `Delta elpd` is not yet reportable.** Its prior-marginalised
   estimator has degenerate importance weights at feasible Monte Carlo sample sizes
   (effective sample size in the low single digits out of `M = 400`), so its standard errors
   are lower bounds. The *sign* of the comparison is stable; no standard-error multiple
   derived from it should be interpreted.

## Citation and contact

Manuscripts describing the method and the clinical application are in preparation. For data
access enquiries, contact the study investigators; for questions about the code, open an
issue.

Saurabh Bhandari — `sbhandari52@uchicago.edu`
