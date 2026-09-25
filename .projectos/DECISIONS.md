# Decisions - grnbench
Append-only ADR-lite. Supersede instead of editing history.

    ### D-001 | YYYY-MM-DD | agent | ACTIVE | short title
    - decision: what was chosen
    - why: reasons / evidence
    - alternatives: what was rejected and why
    - supersedes: D-000 (optional)

STATUS: PROPOSED | ACTIVE | SUPERSEDED by D-xxx | REVERTED


### D-001 | 2026-09-05 | backfill | ACTIVE | Counts or nothing: raw integer counts only
- decision: Ingest rejects anything that does not look like raw integer counts; normalized/log inputs raise IngestError loudly
- why: Preprocessing expects raw integer counts and handles depth normalization explicitly. Feeding normalized data breaks QC gates (mito/ribo %), TF discovery, and strategy reversibility.
- alternatives: Could auto-detect and auto-normalize: loses transparency and couples steps unnecessarily. Could accept but warn: silently wrong answers worse than loud failure.

### D-002 | 2026-09-05 | backfill | ACTIVE | No imputation: strategies cannot smooth or impute data
- decision: The strategy registry is explicitly prohibited from containing smoothing/imputation methods; a test enforces this.
- why: Imputation hides data loss and changes the inference problem. Counts carry more information than smoothed/imputed values. Making the choice visible in preprocessing options is better than burying it.
- alternatives: Could allow imputation as an optional strategy: too easy to use without noticing its assumptions. Could document it as a separate tool: loses the force of design constraints.

### D-003 | 2026-09-05 | backfill | ACTIVE | HVG never drops TFs: kept features = HVGs ∪ TFs detected in data
- decision: HVG filtering is followed by union with TF list; no TF is removed even if not highly variable
- why: A TF can be present in the data but show low variance (housekeeping, stochastic dropout in scRNA-seq). Removing it violates the model contract and silently breaks inference on that regulator. TF list defines the biological task scope.
- alternatives: Could require users to filter TF list: places burden on user and error is silent. Could autoselect TFs by variance: biologically wrong; causal regulators can be lowly variable.

### D-004 | 2026-09-05 | backfill | ACTIVE | QC gates fail loudly: extreme thresholds abort execution
- decision: Thresholds on cell/gene counts, mitochondrial %, ribosomal %, doublet rate and batch imbalance fail with exit code 2 unless --force is given and recorded.
- why: Silent filtering produces garbage outputs with correct-looking manifests. Loud failures force users to choose: accept risk explicitly or tune thresholds. The force flag is recorded in the manifest.
- alternatives: Could auto-tune thresholds: loses transparency and shifts blame to the tool. Could warn silently: the worst outcome (wrong result, confident user).

### D-005 | 2026-09-05 | backfill | ACTIVE | Full provenance: every step logged with parameters, versions, thresholds, input fingerprints
- decision: Manifests record parameter values, effective gene/TF lists, SHA256 input hashes, dependency versions, candidate universe and seed.
- why: Reproducibility is the baseline for any benchmark claim. Provenance must include enough to recreate the same measurement even if no code is available. Hidden choices (default thresholds, autosort) are the enemy of reproducibility.
- alternatives: Could log only major parameters: misses the edge cases that usually break reproduction. Could assume code version suffices: code can change its defaults between runs.

### D-006 | 2026-09-05 | backfill | ACTIVE | Unavailability is data, not a crash: missing deps stay listed with reason
- decision: Methods with missing external dependencies (arboreto, java/jar, Docker image, broken plugin import) are skipped loudly and stay in the available list with their reason (e.g. "arboreto not installed").
- why: Discovering which methods can run is a runtime property that varies per environment. Silently omitting methods loses information. Explicitly stating "method X unavailable: reason Y" is honest and helps debugging.
- alternatives: Could crash if any requested method is missing: forces the user to install everything. Could silently drop unavailable methods: lies about what was computed.

### D-007 | 2026-09-15 | backfill | ACTIVE | Null-calibrated scores: every metric includes (hits+1)/(n+1) p-value vs fixed-universe permutation
- decision: AUPRC, AUROC, EPR and precision@k are reported alongside p-values from score permutation on the *same* fixed candidate universe. A method can lead the table and still be indistinguishable from its own null.
- why: Leaderboards assert a total order over methods. That order is untested. Null calibration reveals whether observed differences exceed what random ranking produces. Conservative (hits+1)/(n+1) p-values avoid overstating significance.
- alternatives: Could report only observed metrics: no indication of whether differences matter. Could use resampling of the reference instead of the scores: conflates ranking uncertainty with reference incompleteness.

### D-008 | 2026-09-15 | backfill | ACTIVE | Resolved comparisons only: paired t-test with Holm correction over independent experimental units (datasets)
- decision: grnbench compare runs paired Student t on declared independent units, Holm-corrects both directions of every pair, Bonferroni simultaneous intervals, and refuses to call an unresolved pair "equivalent".
- why: BEELINE-style benchmarks have independent units (datasets/organisms), not edges (millions) or cells (billions). Confusing these inflates false positives from 0 resolved to 80%. Holm is valid under dependence if individual p-values are valid.
- alternatives: Could use Wilcoxon signed-rank: different estimand (median vs mean difference). Could use false-discovery rate: less conservative but requires more resolved pairs to publish. Could skip testing: leaderboard asserting a total order anyway.

### D-009 | 2026-09-05 | backfill | ACTIVE | Plugin-first architecture: custom models discoverable via entry point, CLI spec file, or Python file reference
- decision: User methods are registered via three pluggable routes: [project.entry-points."grnbench.methods"] in pyproject.toml, --plugin path/to/model.py:ClassName for scripts, --plugin spec.json for Docker/CLI with env variable placeholders.
- why: No fixed method list means users can add their own without forking the code. All three paths go through the same edge-schema registry, manifests and evaluation, so plugins are first-class citizens. Lowering barriers to benchmarking increases adoption.
- alternatives: Could require packaging as a PyPI package: too much friction. Could only support one plugin path: loses flexibility. Could run plugins in different evaluation branches: violates consistency.

### D-010 | 2026-09-05 | backfill | ACTIVE | Not a leaderboard: geneRNIB handles rankings; grnbench audits ranking validity
- decision: grnbench is positioned as an audit tool, not a ranking system. geneRNIB is the live benchmark with a leaderboard. grnbench asks: "Is this ranking real?" not "Who wins?"
- why: Leaderboards assert a total order. That order is untested and unstable across incomplete gold standards. Benchmarking infrastructure should be pluggable (plugins, strategies, null models); benchmarking philosophy should be honest (admit when design cannot rank).
- alternatives: Could build a leaderboard anyway: then own the false positives. Could try to scale to geneRNIB: different goals and tooling (Python vs modular CI). Could deprecate the tool: leaves a gap for auditing harnesses.

### D-011 | 2026-09-13 | backfill | SUPERSEDED by D-012 | Direction: metascience study on when GRN benchmarks support superiority claims
- decision: Pursue the audit-recommended conditional GO: restricted methodological pilot E0-E4 with gates G0-G5; software/portfolio as fallback
- why: Audit 2026-09-13: current benchmark paper NO-GO (overlap + flawed statistics); metascience direction is the recommended conditional GO
- alternatives: Publish current plugin/null benchmark (NO-GO); software-only / geneRNIB component (engineering GO, lower scientific value)

### D-012 | 2026-09-23 | claude | ACTIVE | Direction: universal GRN data-cleaning tool with usable interface, pre-mounted models and side-by-side comparison
- decision: grnbench is a product: any scRNA-seq input -> intelligent, explained cleaning -> GRN-ready data -> run the built-in models -> compare them, all from a usable interface (not CLI only). Metascience paper is NOT the goal.
- why: User direction 2026-09-23 (explicit correction of the audit-driven D-011 proposal).
- alternatives: D-011 metascience study (rejected by user); CLI-only software (insufficient: interface required)
- supersedes: D-011
- update 2026-09-23: user: 'deve diventare un tool da pubblicare' -> publication target = software tool paper + released package. Interface default chosen by agent (user dismissed the question): local web app `grnbench ui` (Streamlit), revisable.

### D-013 | 2026-09-23 | claude | ACTIVE | Cleaning engine = knowledge base of stages/decisions/options + dataset profile + rule selector + tested candidate pip...
- decision: Method set by user: (1) enumerate cleaning stages, (2) critical/decision points, (3) all solutions per problem, implemented (missing deps listed with reason), (4) dataset profile with visual data analysis + context, (5) literature guidelines per option distilled into a machine-readable knowledge base (best-practice combos from top publications), (6) rule-based selector assembles best pipeline + alternatives, runs and scores them, shows pros/cons, user can edit. Must be practical and integrable.
- why: User direction 2026-09-23

### D-014 | 2026-09-23 | claude | ACTIVE | Interface = local Streamlit app (grnbench ui), thin layer over the library
- decision: 7-step Streamlit app; every screen reproducible via plan.json + exact grnbench prep commands; logic in grnbench.ui.logic (tested), app smoke-tested with AppTest
- why: User dismissed the interface question; agent default. Publishable tools commonly ship Streamlit; pip-installable; reuses pipeline directly.
- alternatives: FastAPI+React web app; desktop app; terminal wizard

### D-015 | 2026-09-24 | claude | ACTIVE | Validate the engine on real data before trusting synthetic gains
- decision: Headline claims about the automatic cleaning must come from real datasets with reference networks (CELLxGENE + CollecTRI halves) and be paired with the synthetic bench; synthetic-only gains are reported as upper bounds of the simulator's assumptions.
- why: Interim v2: synthetic +42% vs real ~ties; rules alone hurt on real data
