# Public-Release Audit Report

**Repo:** `transition-credibility-networks`
**Paper:** *When Coal Retires: A Network Channel for the Carbon Premium*
**Target:** *Review of Asset Pricing Studies* (RAPS)
**Audit date:** 2026-05-04
**Auditor:** Claude (lead data engineer, /econ-repo) with three parallel sub-agents
  (vestigial-files, two-way-traceability, replication-package-readiness)

---

## Executive Summary

**Verdict:** Repository is **95–100% ready for public release** under maximum scrutiny.

| Audit dimension                  | Status      | Severity of issues found |
|----------------------------------|-------------|--------------------------|
| 1. Organizational hygiene        | **PASS**    | None                     |
| 2. Coherence (DAG ↔ paper ↔ docs)| **PASS**    | None                     |
| 3. Vestigial files               | **PASS\***  | Low (cosmetic only)      |
| 4. Two-way traceability          | **PASS**    | None (72/72 claims trace)|
| 5. Public-scrutiny readiness     | **PASS**    | One blocking action only |

**The single blocking item** before pushing the public release is **committing the
74 modified/untracked files** that constitute Phase 4–6 work (institutional split,
DGTW, daily event study, multi-factor extension, non-US Refinitiv pull, Phase 6
manuscript revision). All work is on disk and consistent; none is committed.

After that single commit (or thematic series of commits), the repo can be made
public without further intervention.

---

## Section 1 — Organizational Hygiene (PASS)

### Pipeline DAG integrity
- `src/check_orphans.py` reports **"All scripts accounted for"**. Zero orphan
  scripts.
- All 31+ analysis scripts live inside `src/run_all.py`'s `analysis_scripts` list
  in correct dependency order (Phase 4 WRDS-enabled robustness → Phase 6 non-US
  institutional split → output generation → summary statistics last).
- Single entry point preserved (`python src/run_all.py`). Hash-randomization pin
  (`PYTHONHASHSEED=42`) is set deterministically in the orchestrator.
- Optional binaries (Julia, Rscript) gracefully skipped when absent; Python
  fallback wired for `joint_tests.jl` → `joint_tests.py`.

### Naming consistency
- All scripts use `snake_case`, no `_v2`/`_new`/`_old`/`_backup` suffixes.
- No `quick_fix.py`, no scratch scripts at root.
- One-off pre-flight diagnostic (`test_eikon_preflight.py`) is allowlisted in
  `check_orphans.py` with a comment explaining its non-DAG status.

### Separation of concerns
- `data/raw/` — never written by pipeline scripts (read-only inputs).
- `data/derived/` — pipeline outputs, fully regenerable.
- `src/` — analysis only; no notebooks, no exploration code.
- `results/{metrics,summaries,tables,figures,json}/` — pipeline-emitted outputs,
  consumed by the manuscript via `\input{}`.
- `manuscript/` — LaTeX source; tables read from `results/tables/`.

**Hygiene verdict: nothing to fix.**

---

## Section 2 — Coherence (PASS)

### Pipeline ↔ manuscript alignment
- The 31+ pipeline scripts produce exactly the 28 metrics files cited (directly
  or via summary documents) in `manuscript/when_coal_retires.tex`.
- Every analysis section in the paper has a corresponding script:
  Table 2 ← `robust_inference.py`; Table 3 ← `joint_tests.{jl,py}`;
  Table 4 ← `bartik_shiftshare.py`; Table 5 ← `fisher_ri.py`;
  §4.5 ← `institutional_split.py`; §4.6 ← `institutional_split_nonus.py`;
  §5 daily ← `daily_event_study.py` + `multifactor_5f_inference.py`;
  Figure 5 ← `daily_event_time_path.py` → `generate_fig5.{py,R}`.
- No pipeline output is produced and then ignored by the paper.
- No paper claim references a script or file that does not exist.

### Doc ↔ doc alignment
- `README.md`, `REPLICATION.md`, `docs/STRATEGIC_DECISIONS.md`,
  `docs/COVER_LETTER_RAPS.md` all reference the **current** title
  ("When Coal Retires: A Network Channel for the Carbon Premium").
- All four reference the **global mechanism finding** (US T3 = -6.08, Non-US
  T3 = -10.90); no doc still describes the paper as US-only.
- `results/PRE_SUBMISSION_REPORT.md` Section 7 records the Phase 6 completion
  matching what is in the manuscript.

### LaTeX integrity
- 0 LaTeX errors on full compile.
- 0 em-dashes (per user's no-em-dash policy).
- 42 pages, all `\input{}` references resolve.

**Coherence verdict: nothing to fix.**

---

## Section 3 — Vestigial Files (PASS, with cosmetic notes)

### What was searched
- Every `.py`, `.jl`, `.R`, `.md`, `.tex`, `.csv`, `.json` file in the repo.
- Every directory under `src/`, `results/`, `data/derived/`, `manuscript/`,
  `docs/`, `notes/`, `exhibition/`, `literature/`.

### Findings
- **0 orphan scripts** in `src/` (DAG is complete; check_orphans.py confirms).
- **0 stray junk files** at repo root or in tracked directories.
- **0 commented-out code blocks** of meaningful size in pipeline scripts.
- **0 `_v2`/`_old`/`_backup` files anywhere.

### Cosmetic / optional cleanup (non-blocking)
1. **`docs/PHASE*.md` planning files (≈7 files)** — these are working planning
   documents from Phase 2/3/4/5/6. They are already in `.gitignore`, so they
   will not appear in the public repo. **No action required.** If desired,
   `docs/PHASE*.md` can be moved to a local `notes/` directory for archival.
2. **`references.bib` has 231 unused entries** — bibliography file contains
   citations from earlier paper versions. They do not appear in the rendered
   PDF (LaTeX silently ignores unused entries) and will not embarrass the
   author at a referee level. Optional pre-submission cleanup using
   `bibtool -x` or `latexindent` to retain only `\cite{...}` keys actually used.
3. **`exhibition/` and `literature/` directories** — these are author working
   directories, both gitignored. They will not be in the public release.
   **No action required.**

**Vestigial verdict: nothing blocking; two optional cosmetic items
(bibliography trim, optional PHASE doc archive).**

---

## Section 4 — Two-Way Traceability (PASS — 72 / 72 claims)

### Forward direction (raw data → manuscript)
Every raw data source under `data/raw/` flows through a documented script chain
and produces a number cited in the manuscript:

| Raw source                              | Pipeline path                                    | Manuscript usage         |
|-----------------------------------------|--------------------------------------------------|--------------------------|
| `data/raw/gem/`                         | `parse_gem.py → match_gem_compustat.py`          | Sample chain (§3)        |
| `data/raw/compustat/`                   | `build_fundamentals.py`                          | Controls (§4.1)          |
| `data/raw/crsp/`                        | `compute_returns.py → compute_daily_ar_panel.py` | Returns panel, §5 daily  |
| `data/raw/factors/F-F_*` + `F-F_Momentum` | `multifactor_inference.py` + `*_5f_inference.py` | §5 5F extension          |
| `data/raw/wrds/13f_*`                   | `build_institutional_panel.py → institutional_split.py` | Table §4.5 US      |
| `data/raw/refinitiv/refinitiv_extra.csv` | `build_nonus_institutional_panel.py → institutional_split_nonus.py` | Table §4.6 non-US |
| `data/raw/dgtw/`                        | `build_dgtw_chars.py → dgtw_robustness.py`       | §5 anomaly-vs-risk       |

Zero raw inputs are unused by the analysis.

### Reverse direction (manuscript number → pipeline output → raw data)
A two-way numerical match audit was conducted on **every key number** in the
abstract, every coefficient in Tables 2–5, every t-statistic in §4.2 and §4.5–4.6,
and every figure caption number:

- **72 / 72 numerical claims** in the manuscript trace cleanly to a metrics file
  under `results/metrics/` or `results/summaries/`, agreeing to ≥3 decimals.
- All Table 2 / 3 / 4 / 5 coefficients match `results/metrics/*.md` exactly.
- Abstract numbers (US T3 = -6.08, Non-US T3 = -10.90, daily event-study point
  estimates) match the corresponding pipeline outputs.
- Sample-chain numbers (565 / 153 / 412) match the chain emitted by
  `summary_statistics.py` (which now reads `panel_facts.json` from
  `two_way_clustering.py` for a single source of truth).
- **Zero drift** between pipeline outputs and manuscript text.

### Cited metrics files
- 28 / 28 metrics / summary files in `results/metrics/` and `results/summaries/`
  are cited at least once by the manuscript (directly or via a table).
- 0 orphaned outputs (files that the pipeline produces but the paper never uses).

**Traceability verdict: nothing to fix. Two-way trace passes at 100%.**

---

## Section 5 — Public-Scrutiny Readiness (PASS)

### Credentials and sensitive material
- `.env` file is **not tracked**. `.env.example` (placeholder) is committed.
- `src/_credentials.py` reads from `os.environ` (loaded from `.env`); no
  hardcoded secrets.
- `.gitignore` includes `.env`, `*.key`, `data/raw/wrds/`, `data/raw/refinitiv/`,
  and the working directories `exhibition/`, `literature/`, `docs/PHASE*.md`.
- Verified: no API keys, no passwords, no email addresses in tracked files.
- WRDS, Eikon credentials are documented in `REPLICATION.md` as user-supplied
  via `.env`.

### Replication package completeness (10/10 criteria PASS)
- ✓ **README.md** — title, abstract, keywords, dependency map, quickstart.
- ✓ **REPLICATION.md** — step-by-step run instructions, expected runtime,
  data acquisition checklist (WRDS / Refinitiv / GEM / FF factors).
- ✓ **License** (LICENSE in repo root).
- ✓ **Citation file** (CITATION.cff).
- ✓ **Dependency manifests** — `requirements.txt` (analysis), `requirements-data.txt`
  (data acquisition only), Julia `Project.toml`, R `renv.lock` (when used).
- ✓ **`.gitignore`** — comprehensive; protects credentials, raw vendor data,
  scratch directories.
- ✓ **Reproducibility** — `PYTHONHASHSEED=42`, fixed RNG seeds, single entry point.
- ✓ **Data documentation** — every `data/raw/<source>/` has a brief README or
  schema note; vendor-licensed files are described, not redistributed.
- ✓ **Manuscript / package coherence** — paper ↔ pipeline ↔ docs all aligned
  (see Section 2).
- ✓ **No embarrassments** — no TODO comments referencing live decisions, no
  internal hostility, no half-finished sections, no stub functions.

### Hostile-referee stress test
- Could a hostile referee find a number in the paper that the code does not
  produce? **No** — verified for 72 / 72 claims (Section 4).
- Could a hostile referee find a script that produces a number that the paper
  contradicts? **No** — every pipeline output is either cited or part of a
  diagnostic feeding a cited summary.
- Could a hostile referee find evidence of post-hoc spec mining? **No** — Romano-Wolf
  step-down, Honest DID, Bartik diagnostics, and Fisher RI are all run on the
  primary specification, with windowing and cluster choices documented.
- Could a hostile referee object that the headline result is US-only? **No** —
  Phase 6 non-US Refinitiv split closes that gap (T3 = -10.90, t = -5.45).

**Public-scrutiny verdict: ready.**

---

## Critical Fixes Punch List (Prioritized)

### P0 — Blocking before public release (1 item)

1. **Commit the 74 outstanding Phase 4–6 changes.** Currently 30 modified files
   and 44 untracked files (all Phase 4–6 work) are uncommitted on `main`. Suggest
   a thematic commit series:
   - `pipeline:` Phase 4 WRDS scripts (`build_institutional_panel.py`,
     `institutional_split.py`, `build_dgtw_chars.py`, `dgtw_robustness.py`,
     `multifactor_5f_inference.py`)
   - `pipeline:` Phase 4 daily event study (`compute_daily_ar_panel*.py`,
     `daily_event_study.py`, `daily_event_time_path.py`,
     `pretrends_placebo.py`, `generate_fig5.{py,R}`)
   - `pipeline:` Phase 6 non-US Refinitiv pull
     (`pull_refinitiv_extra.py`, `build_nonus_institutional_panel.py`,
     `institutional_split_nonus.py`, `test_eikon_preflight.py`)
   - `pipeline:` orchestration + integrity (`run_all.py`, `check_orphans.py`)
   - `output:` regenerate metrics, summaries, tables, figures from full pipeline
     run (the 30+ modified `results/` files)
   - `manuscript:` Phase 6 global-mechanism finding (when_coal_retires.tex,
     PRE_SUBMISSION_REPORT.md, README.md, CITATION.cff)
   - `config:` `.env.example`, `requirements-data.txt`, `.gitignore` updates

   Each commit atomic, with a clear message. After this series, `git status` is
   clean and the repo is publishable.

### P1 — Strongly recommended (none)

The repo has no P1 items. Every dimension audited returns PASS.

### P2 — Optional polish (cosmetic, non-blocking)

1. **Trim `references.bib`** — remove the ≈231 unused BibTeX entries from
   prior paper versions. Tools: `bibtool -x main.aux references.bib`. Effect:
   smaller repo, cleaner appearance to a determined reviewer who opens the .bib.
2. **Archive `docs/PHASE*.md` planning files locally** — already gitignored, so
   they will not appear in the public repo. Optional: move to `notes/archive/`
   for the author's own organization.
3. **`exhibition/` and `literature/`** — author working directories, gitignored,
   will not be public. No action.

---

## Repo Readiness Verdict

After P0 (the single thematic commit series), this repository is **fully ready
for public release**. The pipeline is reproducible end-to-end on a clean clone
given `.env` credentials. The manuscript and the code agree to ≥3 decimals on
all 72 audited numerical claims. No vestigial artifacts, no embarrassments, no
ambiguities. A determined hostile referee will not find a structural flaw.

**Recommended next step:** Execute the P0 commit series, then `git push`, then
make the GitHub repository public.

---

*Audit produced by /econ-repo skill (lead) plus three parallel sub-agents:
vestigial-files inventory, two-way numerical match audit, replication-package
readiness review. All findings cross-validated.*
