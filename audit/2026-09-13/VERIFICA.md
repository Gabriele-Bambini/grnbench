# Verifica dell'audit — 13 settembre 2026

**Suite esistente: 199 passed, 3 skipped, 431 warnings; exit 0; 83.52 secondi.**

Runtime: `~\.codex\tmp\grnbench-audit-20260913\Scripts\python.exe`, CPython 3.13.5. Ambiente creato appositamente fuori dal progetto con `uv venv --system-site-packages`, installazione editable di `.[dev]`. Non è un ambiente ermetico; versioni installate in `environment-observed.txt`. Nessuna installazione nel Python globale.

```powershell
rtk proxy uv venv ~\.codex\tmp\grnbench-audit-20260913 --system-site-packages
rtk proxy uv pip install --python ~\.codex\tmp\grnbench-audit-20260913\Scripts\python.exe --link-mode copy -e .[dev]
# Test eseguiti con OMP_NUM_THREADS=2 e OPENBLAS_NUM_THREADS=2.
rtk proxy ~\.codex\tmp\grnbench-audit-20260913\Scripts\python.exe -m pytest -q --junitxml=audit/2026-09-13/pytest.xml
```

Log completo in `pytest.log`, JUnit in `pytest.xml`.

Tre test saltati perché arboreto non installato:

- `tests.test_classical::test_random_forest_matches_genie3_where_both_are_available`
- `tests.test_infer.TestTreeMethods::test_genie3_schema_and_signal`
- `tests.test_infer.TestTreeMethods::test_grnboost2_schema`

GENIE3 e GRNBoost2 quindi non verificati nel runtime dell'audit. ARACNe3 è indicato non disponibile per assenza di ARACNE3_JAR. La CLI `methods` eseguita conferma disponibili Spearman, Pearson, MI, CLR, partial correlation, random forest e lasso; la disponibilità non convalida la correttezza scientifica.

400 warning sono deprecazioni sklearn di `n_alphas`; gli altri riguardano conversioni di indici AnnData e un gruppo pseudobulk piccolo intenzionale nei test. Nessuna riparazione effettuata.

## Controesempi verificati separatamente

```powershell
rtk proxy ~\.codex\tmp\grnbench-audit-20260913\Scripts\python.exe audit/2026-09-13/reproduce_eval_findings.py
rtk proxy ~\.codex\tmp\grnbench-audit-20260913\Scripts\python.exe audit/2026-09-13/reproduce_mi_findings.py
```

- `eval_findings.json`: **6/6 difetti riprodotti**. Exit 0 significa che la fixture conferma il difetto, non che il difetto è stato corretto.
- `reproduce_mi_findings.json`: MI di geni costanti ≈ln(8), CLR=1; punteggi diversi dopo permutazione comune delle cellule. Nessun flip di ranking nella specifica fixture Bernoulli.
- Fixture lette e rieseguite dal coordinatore anche nell'ambiente dei test, dopo l'esecuzione indipendente degli agenti.
- `source_sha256.json`: hash completi dello snapshot di sorgenti, test, benchmark, esempi, script HPC e principali file di configurazione/documentazione. Serve a identificare il codice auditato, non dimostra l'origine dei risultati README.

Questi controesempi dimostrano che il successo della suite corrente non implica correttezza del disegno statistico. Non è stata eseguita una replica BEELINE reale, una calibrazione Monte Carlo del framework riparato o una validazione biologica.

## File e contesto

- Applicazione, test esistenti e README originali non modificati.
- Aggiunti dossier, fixture e risultati in `audit/2026-09-13/`.
- Bootstrap Graphify eseguito secondo istruzioni globali: aggiunti i riferimenti locali in `.graphify-context.json`, `AGENTS.md`, `CLAUDE.md`, `.agent/rules/graphify.md` e registrazione al contesto condiviso. Nessun hook Git installato, perché questo snapshot non contiene `.git`.
- Refresh AST Graphify riuscito: **808 nodi, 1670 archi, 44 comunità**, output `graphify-out/`. Non ha richiesto estrazione semantica tramite LLM.
- Nessun commit, push, messaggio a terzi o modifica alla memoria globale.

Il dossier separa bug riprodotti, limiti metodologici, risultati soltanto dichiarati e proposta futura. La ricerca bibliografica è mirata e datata, non certifica l'assenza assoluta di prior art.
