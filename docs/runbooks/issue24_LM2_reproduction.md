# Runbook issue #24: riproduzione di EXP-20260925-LM-2 (seconda persona, GB) sull'HPC Bocconi

> **Eseguito il 2026-10-08** (job Slurm 690538, `COMPLETED 0:0`, 35 min con 3 worker). Esito: `compare_runs.py`
> `DIFFERENT` sia esatto sia `--rtol 1e-9`, quindi C4 resta `reproduced-computational`. Dettagli nella sezione 4b del record
> e in `benchmarks/results/EXP-20260925-LM-2-repro-GB/`. I passi 1-24 sono stati eseguiti cosi' come scritti, senza STOP.
> Ogni "Atteso" del cluster marcato `UNVERIFIED` e' uscito come previsto: listing live = 20 id registrati,
> `ok=28 ok-eol=20 MISMATCH=0 UNLISTED=0`, `ids+truth identical: True`, `code 1d07221...` nel log, 20/20 righe, 0 errori.
> L'unica differenza e' la durata del job (35 min, stima 1-2 h).
> La sezione 4b del record (PR #82) cita questo file alla commit `0b5e5d0`: dopo uno squash-merge quella commit si raggiunge
> solo con `git fetch origin pull/83/head`; il percorso `docs/runbooks/issue24_LM2_reproduction.md` resta su `main`.
> Causa della differenza, trovata da LM (commento su #82 del 2026-10-10): arrotondamenti float32 diversi fra piattaforme
> nell'X normalizzata, amplificati dai pesi quasi uguali delle coppie reciproche A->B / B->A di Pearson; fix in PR #86.
> Cose imparate, utili alla prossima esecuzione:
> - `hpc.sh run` fa `cd ~/<dir>` prima del comando, quindi una cartella nuova va creata passando da un albero che esiste gia'.
> - Se le cartelle in home sono symlink (`~/gb_runs -> projects/...`), `find ~/X` non le attraversa: usare `find -L ~/X/`.
> - Da una nuova sessione ssh `ps`/`pgrep` possono non vedere i download in `nohup`: per il passo 13 meglio guardare i log
>   (`grep -c "^ok "`, `done: 6 files`).

Preparato il 2026-10-07 leggendo: issue #24 (corpo + tutti e 3 i commenti), `experiments/EXP-20260925-LM-2.md` e `.datasets.tsv`,
`experiments/README.md` regola 8, `experiments/datasets.tsv`, `hpc/hpc.sh`, `hpc/engine_bench.sbatch`, `hpc/fetch_cistarget.sh`,
`hpc/fetch_engine_datasets.sh`, `benchmarks/{fetch_cellxgene,data_manifest,real_bench}.py`, `benchmarks/compare_runs.py` (solo sulla PR #43),
`docs/HPC.md`, la skill `grnbench-hpc`, i record EXP-20260926-LM-1 e EXP-20261002-LM-1 (come lavora Luca sull'HPC).
Letto: `origin/main` = `cfeee83`, PR #43 head = `2cceed0b5820d74f28c63c217d166dac5a7d7dc5` (aperta, non unita).

**Quando e' stato scritto (2026-10-07), nulla di questo runbook era stato eseguito sul cluster** (l'esecuzione del 2026-10-08 e' nella nota in testa). Ogni "Atteso" porta una etichetta:
`(codice)` = ricavato leggendo codice o registro, `(provato)` = eseguito in un sandbox Linux senza cluster (Python 3.11, venv dal lock, `ssh` e `sbatch` finti che
eseguono in locale; vedi sotto), `UNVERIFIED` = non verificabile da qui, lo verifichi tu al passo indicato.
Provato nel sandbox (revisione scettica del 2026-10-07): le quote e le espansioni di `~`/`$HOME` dei comandi `hpc.sh run/submit/fetch` dei passi 6-8, 11, 12 (solo il confronto `diff`), 15, 17, 20, 22;
lo sbatch di `245ac12` che passa a `real_bench.py` di `1d07221` esattamente i flag che questo accetta (con `MODELS` `a+b+c` che diventa `a,b,c`, `JOBS=3`, `CODE_COMMIT` copiato);
`real_bench.py run --help` e il fallimento `unrecognized arguments: --max-attempts` (F1); `plan` e `report --no-final-report` di `1d07221` (stringhe attese); `compare_runs.py` della PR #43 (`2cceed0`)
sui risultati committati, esatto e `--rtol 1e-9` (output del passo 23); `data_manifest.py --check` (exit 0 anche con cartella vuota: serve il conteggio del passo 14).
Non provato: tutto cio' che tocca la rete (CELLxGENE, aertslab), Slurm, conda del cluster, `ssh` vero.

## Come si usa

- Un solo terminale bash (Git Bash su Windows) sulla TUA macchina, VPN Bocconi su, chiave gia' caricata da te in ssh-agent (`ssh-add`: lo digiti tu, la passphrase la scrivi solo tu;
  un agente non esegue mai `ssh-add`, non scrive e non chiede passphrase o password: AGENTS.md sezione 10). Se chiudi il
  terminale rifai il blocco variabili del passo 2 e, se il job e' gia' partito, `JOB=<numero scritto nel commento di lancio su #24>` (o letto con `bash hpc/hpc.sh status` da `$T`).
  Dopo il passo 24 fai `unset PYTHONPATH`.
- Ogni passo ha: **Comando**, **Atteso**, **STOP se**. Se scatta uno STOP: ti fermi, scrivi su #24 cosa hai visto (senza utente, matricola,
  nome di host, percorsi home: usa `~/...`) e aspetti. Non improvvisare.
- Un agente puo' digitare i comandi, ma i passi li decidi e li controlli tu (regola 8 di `experiments/README.md`, decisione di Luca su #24 del 2026-10-03).

## Regole che valgono per tutto il runbook

1. `Permission denied` o `Not allowed at this time`: STOP. `Not allowed` = lockout, aspetta 10 minuti, correggi la causa, prova UNA volta. Mai un ciclo di retry sull'autenticazione.
2. Nessun calcolo sul nodo di login (solo download, sha256, `sbatch`, `squeue`, file piccoli). Il calcolo e' nel job. Le uniche cose leggere che questo runbook lancia sul login sono quelle che la skill `grnbench-hpc` e `docs/HPC.md` gia' prevedono: `plan` (passo 15), `report` su 20 righe (passo 21), `run --help` (passo 10) e, solo se l'env manca, il self-test di `setup_env.sh` (passo 7).
3. Un solo job per esperimento (QOS `stud`: 1 in esecuzione + 1 in coda, 8 CPU). Non lanciare altro (per esempio l'esperimento di #52) mentre questo e' in coda o gira.
4. Mai rilanciare se esiste una issue `running` per questa riproduzione: si commenta l'issue (passo 3 e 16).
5. Nessun ID/matricola, utente, passphrase, percorso assoluto della home, nome di host in repo, record, PR, commenti: scrivi `~/...` e `<alias>`.
   Attenzione: l'output di `hpc.sh check/sync`, la riga `node ...` del log e `tasks.json` contengono utente/host/percorsi. Restano nel tuo terminale.
6. I numeri vengono dagli script. Niente numeri scritti a mano, nessun file di risultato modificato.
7. Mai scrivere nella sezione 1 del record, ne' in `benchmarks/results/EXP-20260925-LM-2/` (cartella di Luca).
8. I riferimenti CollecTRI NON si riscaricano da OmniPath/decoupler (D-020: i byte rigenerati sono diversi, il controllo sha256 fallirebbe per costruzione).
   Il punto 2 del corpo di #24 ("references step di `fetch_engine_datasets.sh`") e' superato dal commento del 2026-10-02 e dalla risposta di Luca del 2026-10-03. Non usare `fetch_engine_datasets.sh`.

## Tre cose trovate leggendo, che cambiano il piano scritto in #24

**F1 (bloccante, risolta qui). Lo sbatch di `main` di oggi non puo' pilotare il codice di `1d07221`.**
Il commento del 2026-10-02 dice che `real_bench.py` a `1d07221` accetta tutti i flag passati dallo sbatch di `main`. Era vero quel giorno (`main` = `245ac12`).
Dopo, #29 (`252d28d`, 2026-10-02) ha aggiunto allo sbatch `--max-attempts "$MAX_ATTEMPTS"`, e `real_bench.py run` a `1d07221` non ha quel flag (argparse: `unrecognized arguments`).
Con `cfeee83` il job fallirebbe in pochi secondi (provato lanciando `real_bench.py run ... --max-attempts 2` del codice di `1d07221`: `error: unrecognized arguments: --max-attempts 2`, exit 2). Soluzione di questo runbook: la "copia di strumenti" sul cluster e' `245ac12` (l'ultimo `main` il cui sbatch conosce `MODELS`,
`KEEP_ARTIFACTS`, `REPO_DIR`, `PYTHONPATH`, `CODE_COMMIT` ma non `--max-attempts`) in una cartella dedicata `~/grnbench-lm2-tools`. Non tocca `~/grnbench`.
Verificato (codice): fra `245ac12` e `origin/main`, `hpc.sh`, `fetch_cistarget.sh`, `setup_env.sh`, `fetch_cellxgene.py` sono identici; lo sbatch differisce solo per `MAX_ATTEMPTS`/`--max-attempts`
e un commento; `data_manifest.py` solo nel parser degli argomenti (stessa semantica di `--check`) e in una riga di `SOURCES` (colonna `source`, irrilevante per `--check`); le 48 righe del registro che servono qui sono identiche
(anche `EXP-20260925-LM-2.datasets.tsv`, `.excluded_ids.txt` e i risultati committati di LM-2 sono identici). Il passo 1 lo ricontrolla sul tuo clone;
se un giorno lo sbatch di `main` tornasse compatibile (il conteggio del passo 1 diventa 0) puoi usare `origin/main` al posto di `245ac12`. In ogni caso va scritto come deviazione nella sezione 4b.

**F2. `benchmarks/compare_runs.py` NON e' su `main`**: esiste solo nella PR #43 (aperta, review in corso). Il passo 23 dice come averlo.

**F3 (solo da annotare). Il testo della sezione 1 non coincide con il manifest.** La sezione 1 dice "human 16 / mouse 4, 8 protocols"; `EXP-20260925-LM-2.datasets.tsv` ha 13 umani / 7 topo / 7 protocolli.
Il manifest e gli sha256 sono la registrazione; la sezione 1 non si modifica (regola). Scrivilo in 4b come osservazione e avvisa Luca su #24.

## Mappa dei nomi

| cosa | locale (tua macchina) | cluster |
|---|---|---|
| clone principale (board, record, compare) | `$MAIN` | - |
| strumenti: `245ac12` | worktree `$T` (`../grnbench-lm2-tools`) | `~/grnbench-lm2-tools` |
| codice da riprodurre: `1d07221` | worktree `$R` (`../grnbench-lm2-repro`) | `~/grnbench-lm2-repro` |
| `compare_runs.py` (PR #43) | worktree `$C` (`../grnbench-compare`) | - |
| dati (cartella dedicata, solo questo esperimento) | - | `~/grnbench_data_lm2repro` |
| output del job | `$RUN` = `$HOME/grnbench_repro/EXP-20260925-LM-2-repro-GB` (fuori dal repo) | `~/gb_runs/EXP-20260925-LM-2-repro-GB` |
| scratch dei worker | - | `~/gb_tmp` |

Tempi indicativi: preparazione locale e cluster 20-30 min; download 10-40 min (UNVERIFIED, ~2 GB); job 40-70 min sul portatile di Luca con 3 worker, sul cluster 35 min con 3 worker (misurato il 2026-10-08, job 690538; EXP-20260926-LM-1: 1h10 con 8 worker per 4 modelli pesanti); confronto 10 min.

---

# FASE A - Prima di toccare il cluster (macchina locale)

## 1. Identita' dei commit, lock, compatibilita' sbatch

**Comando** (dal clone principale):
```bash
git fetch origin
git rev-parse '1d07221^{commit}' '245ac12^{commit}'
git rev-parse '1d07221^{tree}' '0bd3e99^{tree}'
git diff --quiet 1d07221 origin/main -- requirements/lock-py311.txt pyproject.toml && echo LOCK-IDENTICAL
for r in origin/main:hpc/engine_bench.sbatch 245ac12:hpc/engine_bench.sbatch 1d07221:benchmarks/real_bench.py; do printf '%s -> ' "$r"; git show "$r" | grep -c -- '--max-attempts'; done
```
**Atteso** (codice):
```
1d0722179f9109ddcc93eb5665b2deefd52d4106
245ac12729c162cb88023256243ce78abe5ee7c8
c972094fc4c14ce5c90db5a67cf8e3be3a454f35
c972094fc4c14ce5c90db5a67cf8e3be3a454f35
LOCK-IDENTICAL
origin/main:hpc/engine_bench.sbatch -> 1
245ac12:hpc/engine_bench.sbatch -> 0
1d07221:benchmarks/real_bench.py -> 0
```
Le due righe `tree` uguali dicono che `1d07221` (squash) ha lo stesso albero della testa della PR #7 (`0bd3e99`). Se `0bd3e99` non e' nel tuo clone salta quella riga (l'ho controllata io).
Nota: il record (sezione 4) cita `96c9f74` come "PR #7 head"; il suo albero e' diverso da `1d07221`, ma fra i due cambiano solo record, README, CLAIMS e CHANGELOG, mai `src/`, `benchmarks/`, `hpc/`, `requirements/` (verificato con `git diff 96c9f74 1d07221 --stat`).

**STOP se**: `LOCK-IDENTICAL` non compare (l'env costruito dal lock di `main` non sarebbe quello del codice sotto test); le due righe `tree` differiscono; il conteggio di `245ac12` o di `1d07221` non e' 0.
Se quello di `origin/main` e' 0 vedi F1 (puoi usare `origin/main` come strumenti, sostituendo `245ac12` ovunque e annotandolo).

## 2. Variabili, worktree, `.env.hpc`

**Comando** (dal clone principale; `.env.hpc` esiste gia' e punta al tuo alias, vedi `docs/HPC.md` sezione 1):
```bash
MAIN="$PWD"; T="$(dirname "$MAIN")/grnbench-lm2-tools"; R="$(dirname "$MAIN")/grnbench-lm2-repro"; C="$(dirname "$MAIN")/grnbench-compare"
RUN="$HOME/grnbench_repro/EXP-20260925-LM-2-repro-GB"
A='source /software/miniconda3/etc/profile.d/conda.sh && conda activate grnbench'
git worktree add "$T" 245ac12
git worktree add "$R" 1d07221
for d in "$T" "$R"; do grep -v '^GRNBENCH_HPC_DIR=' "$MAIN/.env.hpc" > "$d/.env.hpc"; echo "GRNBENCH_HPC_DIR=$(basename "$d")" >> "$d/.env.hpc"; done
git -C "$T" status --short; git -C "$R" status --short
git -C "$R" check-ignore .env.hpc; tail -n 1 "$T/.env.hpc" "$R/.env.hpc"
```
**Atteso**: i due `status --short` vuoti; `check-ignore` stampa `.env.hpc`; le ultime righe sono `GRNBENCH_HPC_DIR=grnbench-lm2-tools` e `GRNBENCH_HPC_DIR=grnbench-lm2-repro`.
**STOP se**: `status` non e' vuoto (per esempio conversioni di fine riga): non sincronizzare un albero sporco, il codice sotto test deve essere pulito. `.env.hpc` assente: `cp .env.hpc.example .env.hpc` e imposti TU l'alias (`GRNBENCH_HPC_HOST`).
Non committare mai `.env.hpc`.

## 3. Board, issue, PR #43

**Comando** (clone principale, sul branch `main`: AGENTS.md sezione 3 dice `git pull` su `main`; se hai lavoro su un altro branch, non cambiarlo ora: `git fetch origin` e leggi `git show origin/main:.projectos/STATE.md`):
```bash
git pull && python -m grnbench.board
gh issue view 24 --comments | tail -n 60
gh pr view 43 --json state,mergedAt,headRefOid
```
**Atteso**: la board mostra `nothing running` e #24 fra gli esperimenti reclamati; nessun commento nuovo su #24 oltre ai 3 noti (2026-10-02 x2, 2026-10-03 di Luca); la PR #43 e' `OPEN` con `headRefOid` `2cceed0b...` (al momento della lettura) oppure `MERGED`.
**STOP se**: esiste una issue `running` per questa riproduzione (o un tuo altro job): commenta l'issue, non rilanciare. `gh` non e' autenticato (`gh auth status`): fallo prima. Un commento nuovo cambia il piano: leggilo e riallinea il runbook.

## 4. Python locale pronto per il confronto

**Comando**: attiva il venv del repo (creato dal lock, AGENTS.md sezione 4) e:
```bash
python -c "import grnbench, anndata, numpy, pandas, scipy; print('deps ok')"
```
**Atteso** (provato): `deps ok`. `import grnbench` tira dentro anndata e il resto dello stack: un venv con solo numpy/pandas/scipy fallisce al passo 23. **STOP se**: `ModuleNotFoundError`: costruisci il venv con i comandi esatti di AGENTS.md sezione 4 (`pip install -r requirements/lock-py311.txt`, `pip install --no-deps -e .`) e riprova. Non e' un esperimento: non serve la suite di test.

---

# FASE B - Cluster: accesso, codice, riferimenti, ambiente

Da qui, salvo dove scritto, i comandi `bash hpc/hpc.sh ...` si lanciano da `cd "$T"` (il worktree degli strumenti).

## 5. Il cluster risponde

**Comando**:
```bash
cd "$T" && bash hpc/hpc.sh check
```
**Atteso**: `connected: ...@...` (utente@host: resta sul tuo terminale), la versione di Slurm, e `no ~/grnbench-lm2-tools yet: run sync`.
**STOP se**: `Could not resolve hostname` / timeout = VPN giu': esegui il comando che `hpc.sh` stampa come prima riga, poi riprova UNA volta. `Permission denied (publickey)`: chiave non nell'agent o `User` sbagliato; la carichi tu (`ssh-add`), un solo nuovo tentativo. `Not allowed at this time`: lockout, aspetta 10 minuti.

## 6. Sincronizza i due alberi

**Comando**:
```bash
cd "$T" && bash hpc/hpc.sh sync
(cd "$R" && bash hpc/hpc.sh sync)
cd "$T" && bash hpc/hpc.sh run 'cat GIT_COMMIT' && (cd "$R" && bash hpc/hpc.sh run 'cat GIT_COMMIT')
```
**Atteso** (codice): `synced 245ac12729c162cb88023256243ce78abe5ee7c8 to <host>:~/grnbench-lm2-tools`, poi `synced 1d0722179f9109ddcc93eb5665b2deefd52d4106 to <host>:~/grnbench-lm2-repro`, e le due stesse stringhe da `cat`. Nessun suffisso `-dirty`.
**STOP se**: compare `-dirty`; `tar`/ssh danno errore. (`sync` non cancella file remoti: le due cartelle sono nuove, quindi sono pulite.)
Nota: `hpc.sh run` fa `cd ~/<dir>` prima di tutto, quindi fallisce finche' la cartella non esiste: per questo il sync viene prima di ogni altro `run`.

## 7. Stato remoto: coda, spazio, cartelle libere, env, file di riferimento

**Comando** (da `$T`):
```bash
bash hpc/hpc.sh status
bash hpc/hpc.sh run 'df -h ~ | tail -n 1'
bash hpc/hpc.sh run 'ls -d ~/grnbench_data_lm2repro ~/gb_runs/EXP-20260925-LM-2-repro-GB 2>&1; true'
bash hpc/hpc.sh run 'source /software/miniconda3/etc/profile.d/conda.sh && conda env list | cut -d" " -f1 | grep -cx grnbench'
bash hpc/hpc.sh run 'hostname; ls -l /tmp/grnbench_registered_references_LM/'
```
**Atteso**: coda vuota (solo l'intestazione); spazio libero di ordine GB (nell'esperimento LM-1 gli artefatti sono stati 190 MB; qui servono ~2 GB di dati piu' qualche centinaio di MB, UNVERIFIED se esiste una quota per utente); le due cartelle `No such file or directory`; `1` per l'env;
l'`hostname` e' lo stesso nodo di login indicato da Luca nel suo commento del 2026-10-03 su #24 (confronto a occhio, non incollarlo), con `collectri_human.tsv` (521401 byte) e `collectri_mouse.tsv` (529071 byte) (byte dal registro).
**STOP se**:
- c'e' gia' un job in coda che non e' questo (non annullare quello degli altri, e non il tuo senza capire);
- una delle due cartelle esiste: non riusarla alla cieca (potrebbe contenere righe di un altro commit). Cambia i nomi ovunque (suffisso `-2`) e annotalo;
- il conteggio dell'env non e' `0` ne' `1`: guarda `conda env list` a mano e riporta;
- `0` per l'env: costruiscilo con `bash hpc/hpc.sh run 'bash hpc/setup_env.sh'` (primo uso, documentato in `docs/HPC.md` sezione 2; alla fine gira un self-test di circa 1 minuto, l'unica eccezione al "niente calcolo sul login"); attesa `SETUP_DONE`;
- il nodo e' diverso da quello di Luca o `/tmp/grnbench_registered_references_LM` manca: chiedi a Luca su #24 di mettere i file sul tuo nodo (lo ha offerto). NON generare i riferimenti da OmniPath.

## 8. Copia dei riferimenti registrati e controllo sha256 (fallo subito: `/tmp` puo' essere ripulita)

**Comando**:
```bash
bash hpc/hpc.sh run 'mkdir -p ~/grnbench_data_lm2repro/references && cp /tmp/grnbench_registered_references_LM/collectri_human.tsv /tmp/grnbench_registered_references_LM/collectri_mouse.tsv ~/grnbench_data_lm2repro/references/ && cd ~/grnbench_data_lm2repro/references && sha256sum collectri_human.tsv collectri_mouse.tsv && ls -l'
grep -E '^references/collectri_(human|mouse)\.tsv' "$T/experiments/datasets.tsv" | cut -f1,2,3
bash hpc/hpc.sh run 'cd ~/grnbench_data_lm2repro/references && printf "%s  %s\n" 4e5597b16a0e795804e1ee0d47bbb420e144e9c0878f8b2ea05cac34be7e303f collectri_human.tsv 514abd67c4c3b02c896e186cefe1c48070e77c3bafe71df871d0c9e86ba8dab4 collectri_mouse.tsv | sha256sum -c -'
```
**Atteso** (registro; la terza riga e' il controllo automatico, non a occhio: stampa `collectri_human.tsv: OK` e `collectri_mouse.tsv: OK`): le due righe di `sha256sum` coincidono con quelle del registro (seconda riga di comando) e con il commento di Luca:
```
4e5597b16a0e795804e1ee0d47bbb420e144e9c0878f8b2ea05cac34be7e303f  collectri_human.tsv   (521401 byte)
514abd67c4c3b02c896e186cefe1c48070e77c3bafe71df871d0c9e86ba8dab4  collectri_mouse.tsv   (529071 byte)
```
**STOP se**: `sha256sum -c` stampa `FAILED`, una dimensione differisce, o un file manca: STOP e scrivi a Luca su #24. Non rigenerare, non "aggiustare" il registro (D-020: un hash registrato non si sostituisce mai).
Lascia queste copie dove sono: servono anche a #52.

## 9. L'env del cluster coincide con il lock (best effort)

**Perche'**: se il tuo env e' precedente a #3 manca `scikit-image` e il doublet-caller ricade in silenzio sull'euristica: cambierebbero i numeri.
**Comando**:
```bash
bash hpc/hpc.sh run "$A && cd ~/grnbench-lm2-repro && python --version && python -m pip install --dry-run --no-deps -r requirements/lock-py311.txt 2>&1 | grep -v 'Requirement already satisfied' | head -n 20"
```
**Atteso**: `Python 3.11.x` e nessun'altra riga (al massimo avvisi di pip). (UNVERIFIED: richiede pip >= 22.2 nell'env; se `--dry-run` non e' riconosciuto, controlla a mano `python -m pip list --format=freeze | grep -i -E 'scikit-image|scrublet|scanpy|numpy|pandas|scipy'` contro il lock.)
**STOP se**: compare `Would install ...`: l'env non e' quello del lock. Aggiornalo PRIMA del job con `bash hpc/hpc.sh run 'bash hpc/setup_env.sh'` (idempotente), ripeti il comando. Mai rieseguire `setup_env.sh` mentre il job gira (AGENTS.md sezione 13).
Annota la versione di Python per la sezione 4b.

## 10. Il codice sotto test accetta i flag dello sbatch e viene importato dalla cartella giusta

**Comando**:
```bash
bash hpc/hpc.sh run "$A && cd ~/grnbench-lm2-repro && PYTHONPATH=\$HOME/grnbench-lm2-repro/src python benchmarks/real_bench.py run --help | grep -o -- '--[a-z][a-z-]*' | sort -u | tr '\n' ' '; echo"
bash hpc/hpc.sh run "$A && PYTHONPATH=\$HOME/grnbench-lm2-repro/src python -c 'import grnbench; print(grnbench.__file__)'"
bash hpc/hpc.sh run 'grep -c -- "--max-attempts" hpc/engine_bench.sbatch; true'
```
**Atteso** (codice): `--chunk --help --index --jobs --keep-artifacts --models --out --refs --tasks --workroot` (nessun `--max-attempts`); un percorso che finisce con `/grnbench-lm2-repro/src/grnbench/__init__.py` (non incollarlo: contiene il tuo home); `0` per lo sbatch degli strumenti.
**STOP se**: compare `--max-attempts` fra i flag (non e' il codice di `1d07221`), mancano flag, l'import punta a un'altra cartella (l'install editable dell'env prevarrebbe: lo sbatch usa `PYTHONPATH`, ma devi saperlo), o il conteggio sbatch non e' 0 (F1).

---

# FASE C - Dati (nodo di login, solo download e sha256)

## 11. Avvia i due download in background

**Comando** (da `$T`; `;` e non `&&` prima dell'ultimo comando, cosi' ssh torna subito):
```bash
bash hpc/hpc.sh run 'mkdir -p logs ~/grnbench_data_lm2repro; nohup bash -c "source /software/miniconda3/etc/profile.d/conda.sh && conda activate grnbench && python benchmarks/fetch_cellxgene.py --out ~/grnbench_data_lm2repro/cellxgene --from-manifest experiments/EXP-20260925-LM-2.datasets.tsv" > logs/lm2_fetch_cxg.log 2>&1 < /dev/null &'
bash hpc/hpc.sh run 'mkdir -p logs; nohup bash hpc/fetch_cistarget.sh ~/grnbench_data_lm2repro/cistarget > logs/lm2_fetch_cistarget.log 2>&1 < /dev/null &'
```
**Atteso**: i due comandi tornano subito senza output (provato: con un `ssh` finto la forma `nohup bash -c "...~..." > log 2>&1 < /dev/null &` torna in meno di 0,1 s ed espande `~` nel `bash -c` interno). (UNVERIFIED: il comportamento del `ssh` vero con processi in background; se non torna entro ~15 s, Ctrl-C e' sicuro, `nohup` lascia andare il download: verifica col passo 13.)
**STOP se**: errore ssh (regole del passo 5).

## 12. Selezione live di CELLxGENE contro i 20 id registrati (non sovrascrive nulla)

**Comando** (cartella `--out` DIVERSA dai dati: il dry-run scrive `manifest.tsv`):
```bash
bash hpc/hpc.sh run "$A && python benchmarks/fetch_cellxgene.py --out ~/gb_tmp/lm2_listing --max 20 --max-mb 60 --organisms 'Homo sapiens,Mus musculus' --exclude-ids experiments/EXP-20260925-LM-2.excluded_ids.txt --exclude-assay targeted --dry-run"
bash hpc/hpc.sh run 'cut -f2,10 ~/gb_tmp/lm2_listing/manifest.tsv | sort > ~/gb_tmp/lm2_listing/live.txt; cut -f2,10 experiments/EXP-20260925-LM-2.datasets.tsv | tr -d "\r" | sort > ~/gb_tmp/lm2_listing/registered.txt; diff ~/gb_tmp/lm2_listing/live.txt ~/gb_tmp/lm2_listing/registered.txt && echo SAME-20-IDS-AND-URLS'
```
**Atteso**: `20 datasets, 677 MB total` seguito da 20 righe di tabella, poi `SAME-20-IDS-AND-URLS` (come nella prima riproduzione, record 4 passo 1). (UNVERIFIED oggi: il listing pubblico puo' essere cambiato.)
**Se differisce** (non blocca il download, che usa il manifest registrato e controlla gli sha256): annota le righe `diff` per 4b e avvisa su #24; vuol dire che la procedura di selezione non e' piu' riproducibile dal listing, non che i dati registrati siano sbagliati (mia lettura, UNVERIFIED rispetto a Luca).
**STOP se**: errore HTTP (403/timeout) sull'API: riportalo; senza rete verso CELLxGENE falliranno anche i download del passo 11.

## 13. Attendi i download

**Comando** (ripeti ogni 5-10 minuti):
```bash
bash hpc/hpc.sh run 'grep -c "^ok " logs/lm2_fetch_cxg.log; tail -n 1 logs/lm2_fetch_cistarget.log; du -sh ~/grnbench_data_lm2repro/* ; ps -u "$(whoami)" -o etime,args | grep -E "[f]etch_cellxgene|[f]etch_cistarget|[c]url"; true'
```
**Atteso** (codice): alla fine `20`; l'ultima riga del log cisTarget `done: 6 files in ...`; `cellxgene` ~0.7 GB, `cistarget` ~1.3 GB. Il log CELLxGENE inizia con `20 datasets, 677 MB total`. Il segnale di fine sono i log (`20` e `done: 6 files`), non l'elenco dei processi: da una nuova sessione ssh `ps` puo' non vedere i download in `nohup` (visto il 2026-10-08), quindi "nessun processo elencato" non prova che abbiano finito.
**STOP se**: `Traceback`, `HTTP Error`, `curl: (22)` nei log (`tail -n 20 logs/lm2_fetch_*.log`). Il download riprende da dove era (salta i file completi): rilancia UNA volta il comando del passo 11 che si e' fermato (quello CELLxGENE o quello cisTarget; non entrambi alla cieca se uno e' finito); se fallisce ancora riporta su #24. Nessun dataset viene sostituito da un altro (regola di esclusione a): un dataset mancante = esito da riportare, non da aggirare.

## 14. Controllo dei dati contro i sha256 registrati

**Comando**:
```bash
bash hpc/hpc.sh run "$A && python benchmarks/data_manifest.py ~/grnbench_data_lm2repro --check experiments/datasets.tsv"; echo "exit=$?"
```
**Atteso** (codice + registro): 20 righe `ok-eol  cellxgene/.../truth.json` (tutti i `truth.json` sono registrati con fine riga CRLF, il cluster li scrive LF: accettato dallo script), poi
```
ok=28 ok-eol=20 MISMATCH=0 UNLISTED=0
exit=0
```
cioe' 20 `data.h5ad` + 6 cisTarget + 2 riferimenti = 28 `ok`, 20 `truth.json` `ok-eol`, totale 48 file (come EXP-20260926-LM-1: `ok=30 ok-eol=20` con 4 riferimenti).
**STOP se**: `MISMATCH` > 0 o `UNLISTED` > 0 o `exit=1` (regola di #24: ogni MISMATCH ferma la run; riporta il path e nient'altro, i file restano dov'e'); `ok + ok-eol` < 48 (file mancanti: il controllo senza `--require-all` non li segnala, torna al passo 13). Incolla l'ultima riga in 4b.

---

# FASE D - Piano e lancio

## 15. `plan` col codice sotto test, e controllo degli id

**Comando**:
```bash
bash hpc/hpc.sh run "$A && mkdir -p ~/gb_runs/EXP-20260925-LM-2-repro-GB ~/gb_tmp && cd ~/grnbench-lm2-repro && PYTHONPATH=\$HOME/grnbench-lm2-repro/src python benchmarks/real_bench.py plan --data-root ~/grnbench_data_lm2repro --out ~/gb_runs/EXP-20260925-LM-2-repro-GB/tasks.json"
cd "$T" && bash hpc/hpc.sh run 'cat ~/gb_runs/EXP-20260925-LM-2-repro-GB/tasks.json' | python -c '
import json, sys
new = json.load(sys.stdin)
old = json.load(open("benchmarks/results/EXP-20260925-LM-2/tasks.json", encoding="utf-8"))
key = lambda ts: [(t["id"], json.dumps(t["truth"], sort_keys=True)) for t in ts]
print("rerun", len(new), "committed", len(old), "ids+truth identical:", key(new) == key(old))
'
```
**Atteso**: `20 real-data tasks -> .../tasks.json`, poi `rerun 20 committed 20 ids+truth identical: True` (come record 4 passo 3).
**STOP se**: il numero non e' 20 (la cartella dati non e' dedicata o manca un dataset) o `False`. `tasks.json` contiene percorsi con il tuo home: non committarlo e non incollarlo.

## 16. Ultimo controllo prima del lancio

**Comando**:
```bash
cd "$MAIN" && python -m grnbench.board
cd "$T" && bash hpc/hpc.sh status
```
(Facoltativo, non crea job: `bash hpc/hpc.sh submit '--test-only <stessi argomenti del passo 17>'`; UNVERIFIED se il cluster lo permette; se risponde con un errore QOS, vale come STOP.)
**Atteso**: nessuna issue `running` (la tua arrivera' dopo), coda vuota.
**STOP se**: appare una issue `running` per questa riproduzione o un job in coda: NON rilanciare, commenta l'issue.

## 17. Lancio: un solo job, poi label `running` e commento

**Comando** (da `$T`; gli argomenti in apice singolo, cosi' `$HOME` lo espande il cluster e non la tua macchina):
```bash
cd "$T"
JOB=$(bash hpc/hpc.sh submit '--time=24:00:00 --export=ALL,REPO_DIR=$HOME/grnbench-lm2-repro,OUT=$HOME/gb_runs/EXP-20260925-LM-2-repro-GB,SCRIPT=benchmarks/real_bench.py,TASKS=$HOME/gb_runs/EXP-20260925-LM-2-repro-GB/tasks.json,REFS=$HOME/grnbench_data_lm2repro/references,CHUNK=100000,JOBS=3,MODELS=pearson_signed+pearson_motif+motif_prior,KEEP_ARTIFACTS=1,GRNBENCH_CISTARGET_DIR=$HOME/grnbench_data_lm2repro/cistarget,TMPDIR=$HOME/gb_tmp hpc/engine_bench.sbatch')
echo "job: [$JOB]"
```
**Atteso**: `job: [<numero>]` (solo cifre; se c'e' `;cluster` usa la parte prima del `;`). Le variabili sono quelle del commento di #24 piu' `OUT`/`TASKS`/`REFS`/`CHUNK=100000` (un solo job, indice 0) usate da Luca; `--time=24:00:00` e' lo standard del team.
**STOP se**: `JOB` vuoto o errori `QOSMax...` (hai gia' un job: non sottomettere altro). Se il job risulta `FAILED` dopo pochi secondi vedi la tabella "Se qualcosa va storto" in fondo: non ha scritto nulla, si puo' risottomettere lo stesso comando corretto.

Subito dopo (AGENTS.md sezione 6; l'issue esiste gia' ed e' tua, il record e' gia' su `main`):
```bash
gh issue edit 24 --add-label running
gh issue comment 24 --body "launched: Bocconi HPC (partition/QOS stud), Slurm job $JOB, 3 dataset workers (8 CPUs requested). started: <UTC>. expected finish: <UTC>, UNVERIFIED estimate 1-2 h (the laptop run took ~70 min). code under test: 1d0722179f9109ddcc93eb5665b2deefd52d4106 in ~/grnbench-lm2-repro; tooling: 245ac12729c162cb88023256243ce78abe5ee7c8 in ~/grnbench-lm2-tools (main at the time of the plan; current main's sbatch passes --max-attempts, which 1d07221 does not accept). data: ~/grnbench_data_lm2repro, checked with data_manifest.py --check (last line: <incolla>). outputs: ~/gb_runs/EXP-20260925-LM-2-repro-GB. check: bash hpc/hpc.sh status $JOB"
```
**STOP se**: `gh` non autenticato. Nessun nome di host, utente o percorso home nel commento.

---

# FASE E - Il job

## 18. Controllo iniziale (2-3 minuti dopo il lancio)

**Comando**:
```bash
bash hpc/hpc.sh status "$JOB"
bash hpc/hpc.sh log "$JOB" 40
```
**Atteso** (codice): stato `RUNNING` (o `PENDING` se c'e' coda; vedi sotto). Nel log, la riga 1 `node <nome>  cpus 3  task 0  chunk 100000` (il nome del nodo resta nel terminale) e la riga 2:
```
code 1d0722179f9109ddcc93eb5665b2deefd52d4106  python Python 3.11.x  models pearson_signed+pearson_motif+motif_prior
```
poi, via via che finiscono, righe `cxg__...: <secondi>s` (i piu' lenti ~5 min; la prima riga puo' tardare qualche minuto).
**STOP se**: `CODE_COMMIT` diverso da `1d07221...` o modelli diversi dai tre: `bash hpc/hpc.sh cancel "$JOB"` e riporta. `unrecognized arguments` = hai usato lo sbatch sbagliato (F1). `No space left on device` = `TMPDIR` non applicato (precedente LM-1 deviazione 1). `cannot activate conda env` = env mancante.
`PENDING` con `QOSMaxJobsPerUserLimit` / `QOSMaxSubmitJobPerUserLimit`: hai gia' un job; stop. `PENDING (Resources)`/`(Priority)` e' normale: commenta l'ETA (`bash hpc/hpc.sh run 'squeue --me --start'`) su #24 e chiudi la sessione; riprendi dal passo 19.

## 19. Attesa con polling (mai uno sleep cieco)

**Comando** (nel tuo terminale, oppure fai controllare all'agente ogni 10-15 minuti):
```bash
cd "$T"; d=60
while :; do
  out=$(bash hpc/hpc.sh status "$JOB") || { echo "ssh fallito: STOP e riporta"; break; }
  rows=$(bash hpc/hpc.sh run 'ls ~/gb_runs/EXP-20260925-LM-2-repro-GB/rows 2>/dev/null | wc -l') || { echo "ssh fallito: STOP"; break; }
  echo "$(date -u +%H:%M:%S) rows=$rows"; echo "$out" | tail -n 3
  echo "$out" | grep -Eq 'COMPLETED|FAILED|CANCELLED|TIMEOUT|OUT_OF_ME|NODE_FAIL' && break
  echo "$out" | grep -q 'QOSMax' && { echo "bloccato dal limite QOS: STOP e riporta"; break; }
  sleep $d; d=$(( d < 900 ? d * 2 : 900 ))
done
```
Ogni ora circa, o a ogni cambio di stato: `gh issue comment 24 --body "checked: RUNNING <hh:mm>, <n>/20 rows"` (AGENTS.md sezione 6).
**Atteso**: `rows` sale fino a `20`; alla fine una riga `sacct` con `COMPLETED` (`0:0`). Tempo: 40-70 min sul portatile; sul cluster UNVERIFIED.
**STOP se**: nessuna riga nuova da 60 minuti (a `1d07221` un worker ucciso per memoria puo' bloccare il pool: lo ha corretto #37 dopo): leggi `bash hpc/hpc.sh log "$JOB" 100`, riporta su #24, e solo dopo decidi `cancel` + risottomissione dello stesso comando (riprende dalle righe gia' scritte);
un errore ssh nel ciclo (il ciclo si ferma da solo: non rilanciarlo in loop se e' un errore di autenticazione).

## 20. Il job e' finito: lo stato di Slurm non basta

**Comando**:
```bash
bash hpc/hpc.sh status "$JOB"
bash hpc/hpc.sh run 'ls ~/gb_runs/EXP-20260925-LM-2-repro-GB/rows | wc -l; cat ~/gb_runs/EXP-20260925-LM-2-repro-GB/CODE_COMMIT; grep -L "\"error\": null" ~/gb_runs/EXP-20260925-LM-2-repro-GB/rows/*.json | wc -l; ls ~/gb_runs/EXP-20260925-LM-2-repro-GB/artifacts | wc -l'
bash hpc/hpc.sh log "$JOB" 5000 | grep -c -E '^cxg__.*: [0-9.]+s$'
bash hpc/hpc.sh log "$JOB" 5000 | grep -c Traceback
```
**Atteso**: `COMPLETED` / `0:0` e il tempo trascorso (`Elapsed`); poi, in ordine, `20`, `1d0722179f9109ddcc93eb5665b2deefd52d4106`, `0`, `20` (artefatti: implicito dal n = 20 dell'endpoint primario, UNVERIFIED); `20`; `0`.
`real_bench.py run` esce 0 anche se una riga ha un errore: per questo si guardano le righe.
A `1d07221` NON esistono `COMPLETE.json`, `RUN_CONFIG.json` ne' `--partial` (le ha introdotte #29 su `main`; la skill `grnbench-hpc` e `docs/HPC.md` di `main` li citano): la loro assenza qui NON e' uno STOP. Il controllo di completezza di questo runbook sono le righe sopra (20 righe, 0 errori) e il passo 14 per i dati.
**STOP se**: righe diverse da 20; un `"error"` non nullo (guarda quale, senza incollare percorsi: `grep -o '"error": "[^"]\{0,200\}' <file>`); `TIMEOUT`/`NODE_FAIL`/`OUT_OF_ME`/`FAILED`. Regole:
- `TIMEOUT` o `NODE_FAIL`: stesso comando del passo 17 una volta (il runner salta le righe esistenti), annotalo come deviazione operativa.
- Righe con errore per cause di infrastruttura (disco pieno, OOM): come LM-1 deviazione 1, sposta (non cancellare) le righe errate in `attempt1_<jobid>/` nella cartella di output, avvisa su #24, risottomissione UNA volta. Mai per errori del codice o dei dati: sono un risultato.
- Un errore che si ripete: e' un'osservazione da registrare, non da nascondere.

---

# FASE F - Risultati, confronto, decisione

## 21. Report sul cluster (con il codice sotto test)

**Comando**:
```bash
bash hpc/hpc.sh run "$A && cd ~/grnbench-lm2-repro && PYTHONPATH=\$HOME/grnbench-lm2-repro/src python benchmarks/real_bench.py report --out ~/gb_runs/EXP-20260925-LM-2-repro-GB --no-final-report > /dev/null && head -n 6 ~/gb_runs/EXP-20260925-LM-2-repro-GB/REPORT_REAL.md"
```
**Atteso** (dal report committato, codice): 
```
# Engine mass test on real datasets - 20 datasets

- ran end to end: 20
- refused as expected (not raw counts): 0
- errors: 0
- doublet calling fell back to the heuristic (scrublet failed) in 4/20 datasets (baseline or auto arm)
```
e nel file, nella tabella "Declared comparisons, AUPRC ratio, auto pipeline": `| pearson_motif | pearson_signed | 20 | +0.182 | [+0.138, +0.227] | 1.36e-07 | 1 | pearson_motif |`. Scrive `REPORT_REAL.md` e `all_rows_real.csv`.
Qui il report e' leggero (20 righe, pandas) e la skill `grnbench-hpc` sezione 7 lo lancia sul login. `--no-final-report` salta figure e LaTeX (richiederebbero un secondo job e non entrano nella decisione): in 4b scrivi "figures/LaTeX not regenerated; originals in `benchmarks/results/EXP-20260925-LM-2/viz` and `latex`" (AGENTS.md sezione 7, "Not produced: say why").
**STOP se**: il comando fallisce con un errore del codice (non una deviazione da "sistemare": riportalo); `errors:` non 0.

## 22. Checksum, recupero, verifica

**Comando**:
```bash
cd "$T" && bash hpc/hpc.sh run 'cd ~/gb_runs/EXP-20260925-LM-2-repro-GB && sha256sum rows/*.json REPORT_REAL.md all_rows_real.csv CODE_COMMIT > SHA256SUMS && wc -l < SHA256SUMS && du -sh rows artifacts && du -sh .'
mkdir -p "$RUN"
for f in rows REPORT_REAL.md all_rows_real.csv CODE_COMMIT SHA256SUMS; do bash hpc/hpc.sh fetch "~/gb_runs/EXP-20260925-LM-2-repro-GB/$f" "$RUN" || break; done
(cd "$RUN" && sha256sum -c SHA256SUMS | grep -vc ': OK$'); ls "$RUN/rows" | wc -l
bash hpc/hpc.sh log "$JOB" 5000 | grep '^code ' > "$RUN/job_code_line.txt"; cat "$RUN/job_code_line.txt"
```
**Atteso**: `23` righe di checksum; poi tre righe di `du`: `rows` pochi MB, `artifacts` ordine centinaia di MB (LM-1: 190 MB, UNVERIFIED qui) e il totale `.` (un solo `du -sh . rows artifacts` stampa solo `.`: GNU `du` non rielenca le sottocartelle gia' contate, provato); il `fetch` stampa `fetched ... -> ...` per ognuno dei 5; `0` (nessun checksum non-OK; su macOS usa `shasum -a 256 -c`); `20`; la riga `code 1d07...`.
**STOP se**: un checksum non e' OK o il `tar` si interrompe: rifai il `fetch` di quel pezzo UNA volta. Gli artefatti (grandi) restano sul cluster e non si scaricano: nella sezione 4b il percorso `~/gb_runs/EXP-20260925-LM-2-repro-GB`.
`$RUN` e' fuori dal repo: niente di questo va in git (tranne il testo che incollerai in 4b).

## 23. Procurarsi `compare_runs.py` (non e' su `main` finche' la PR #43 non e' unita)

**Comando** (clone principale):
```bash
cd "$MAIN" && git fetch origin
if git cat-file -e origin/main:benchmarks/compare_runs.py 2>/dev/null; then echo "ON-MAIN: usa un checkout di origin/main"; else echo "NON-SU-MAIN: PR #43"; fi
```
- `ON-MAIN` (la PR e' stata unita): `git worktree add "$C" origin/main` e usa `$C`.
- `NON-SU-MAIN`:
```bash
git fetch origin pull/43/head:pr-43-compare
git rev-parse pr-43-compare
git worktree add "$C" pr-43-compare
```
**Atteso**: `git rev-parse` stampa `2cceed0b5820d74f28c63c217d166dac5a7d7dc5` (testa letta il 2026-10-06). Se e' cambiata: `git diff 2cceed0 pr-43-compare -- benchmarks/compare_runs.py`, leggi cosa cambia (il significato di `primary: MATCH` e' quello su cui poggia la regola di #24) e annota lo sha esatto in 4b.
Nota (letto): fra `origin/main` e la PR #43 cambiano solo `compare_runs.py`, i suoi test, MAP, STATE, CHANGELOG, CONTRIBUTIONS: `src/`, il record e i risultati di LM sono identici, quindi `$C` ha lo stesso record e gli stessi risultati committati di `main`.

**Controllo del tuo ambiente** (positivo, sui file committati contro se stessi; non scrive nulla):
```bash
cd "$C" && export PYTHONPATH=src
python -c "import grnbench; print(grnbench.__file__)"
python benchmarks/compare_runs.py benchmarks/results/EXP-20260925-LM-2 benchmarks/results/EXP-20260925-LM-2 --record experiments/EXP-20260925-LM-2.md
```
**Atteso** (provato nel sandbox con la testa `2cceed0`, Python 3.11.15, venv dal lock, exit 0): il percorso di `grnbench` sta sotto `$C/src`; poi
```
tasks: same (20 committed, 20 rerun)
rows: identical - compared 1180 numeric, 40 bool and 400 other fields; 140 timing fields not compared; exact
report: identical (line endings and timing columns ignored)
primary: MATCH - pearson_motif > pearson_signed | auprc_ratio | auto
  committed: n=20 mean +0.182 [+0.138, +0.227] raw one-sided p(A>B) 3.39e-08
  rerun:     n=20 mean +0.182 [+0.138, +0.227] raw one-sided p(A>B) 3.39e-08
verdict: MATCH
```
Con `--rtol 1e-9` (provato) e' identico tranne la riga `rows:`, che termina con `; tolerance rtol=1e-09 atol=0` invece di `; exact`.
**STOP se**: `ModuleNotFoundError` (rivedi il passo 4), exit 2 o conteggi diversi (`1180 / 40 / 400 / 140`): lo strumento o i risultati committati sul tuo clone non sono quelli del piano; riporta su #24 prima di usarlo per decidere.

## 24. Il confronto vero: due volte (esatto e `--rtol 1e-9`), dal piano di #24

**Comando** (da `$C`, con `PYTHONPATH=src`):
```bash
cd "$C" && export PYTHONPATH=src
python benchmarks/compare_runs.py benchmarks/results/EXP-20260925-LM-2 "$RUN" --record experiments/EXP-20260925-LM-2.md --json "$RUN/compare_exact.json" | tee "$RUN/compare_exact.txt"; echo "exit=${PIPESTATUS[0]}"
python benchmarks/compare_runs.py benchmarks/results/EXP-20260925-LM-2 "$RUN" --record experiments/EXP-20260925-LM-2.md --rtol 1e-9 --json "$RUN/compare_rtol1e-9.json" | tee "$RUN/compare_rtol1e-9.txt"; echo "exit=${PIPESTATUS[0]}"
```
(Il piano diceva `--json <RUN>/compare.json`; con due esecuzioni servono due nomi: annotalo.)
**Atteso, se tutto combacia** (forma provata sui file committati contro se stessi al passo 23; per la tua rerun i valori sono quelli che trovi): per l'esatto, le righe del passo 23 con `rerun:` uguale a `committed:` e `verdict: MATCH`, `exit=0`. Per `--rtol 1e-9` la riga `rows:` termina con `; tolerance rtol=1e-09 atol=0`.
Differenze nel tempo (`seconds`, `*_seconds`) non contano; la riga `report:` e' solo informativa e non entra nel verdetto.
**Come si leggono i casi**:
- esatto `MATCH`: valori identici bit a bit anche su un altro OS/CPU del portatile Windows di Luca (esito forte).
- esatto `DIFFERENT` ma `--rtol 1e-9` `MATCH`: deriva nelle ultime cifre (altro OS/BLAS/CPU del portatile Windows di Luca): e' proprio il caso per cui il piano prevede la seconda esecuzione. Registra entrambi.
- `--rtol 1e-9` con `primary: DIFFERENT` o `NOT EVALUABLE` (exit 1): e' un RISULTATO, non un errore. Non ripetere il confronto cambiando opzioni e non rilanciare il job: riga 2 della tabella del passo 25 (C4 resta `reproduced-computational`), tutto in 4b.
- `exit=2` (cartella senza `rows/`, riga di record/primary non parsabile, JSON di riga illeggibile): e' un errore d'uso, non un verdetto. Leggi il messaggio su stderr, correggi il COMANDO (per esempio il percorso di `$RUN`) e ripeti il confronto; non e' una rerun dell'esperimento.
Non cambiare soglie dopo aver visto il risultato (niente `--atol`, niente altre tolleranze): se serve un'analisi in piu' va su #24, non nella regola.
**STOP se**: `compare_runs.py` non esiste o fa `ImportError` (torna al passo 23); `exit=2` che non si risolve correggendo il comando; i due output (esatto e `--rtol 1e-9`) non sono stati salvati in `$RUN/compare_*.txt` prima di decidere. Se l'esito e' `DIFFERENT` o `NOT EVALUABLE`: riporta su #24 gli output verbatim e poi registra tutto in 4b secondo la tabella del passo 25; non e' un fallimento da riparare.

## 25. Regola di decisione (fissata in #24 prima della run)

| esito del confronto `--rtol 1e-9` | cosa succede |
|---|---|
| riga `primary: MATCH - pearson_motif > pearson_signed \| auprc_ratio \| auto` | C4 diventa **`reproduced`** (decisione di Luca su #24, 2026-10-03, allineata da #35 in `paper/README.md` regola 3, `experiments/README.md` regola 8, AGENTS.md sezione 8: una rerun diretta e controllata da te, con un agente che digita i comandi, su dati recuperati e verificati sha256 da te, con sezione 4 firmata). Ogni `DIFFERENT` in `rows:` o `tasks:` va comunque riportato per intero in 4b e segnalato a Luca nella PR. |
| `primary: DIFFERENT`, `NOT EVALUABLE`, o `tasks:` diverso | C4 resta `reproduced-computational`. La differenza si registra in 4b, senza "ripararla" ne' rilanciare in silenzio (AGENTS.md sezione 7). |
| dati (passo 14) o job (passo 20) falliti, nessun confronto possibile | nessuna conclusione sul claim; l'issue resta aperta con la causa (oppure `failed` se non si puo' proseguire). |

Condizioni per la parola `reproduced` (AGENTS.md sezione 8): hai deciso e controllato tu ogni passo dal record e da #24; i dati li hai copiati/scaricati tu e controllati contro gli sha256 registrati (passi 8 e 14); firmi la sezione 4.

## 26. Cosa scrivere dopo

**Commento finale su #24** (nessun host/utente/percorso home), con: job `$JOB`, stato `sacct` ed `Elapsed`, `code 1d07221...` dal log, ultima riga di `data_manifest.py`, e le righe `rows:`/`primary:`/`verdict:` dei DUE confronti incollate com'e'.

**Record `experiments/EXP-20260925-LM-2.md`**: SOLO sotto la sezione 4 esistente, una sottosezione nuova "4b" (non toccare la sezione 1 ne' il testo gia' presente della 4, salvo le due righe di testa elencate sotto; il testo della prima riproduzione resta com'e'). Branch `gb/<argomento>` da `main`, mai su `main`.
Perche' le chiavi sono "strane" (provato in un clone di scratch: `grnbench report EXP-20260925-LM-2` sul record con la 4b a chiavi `environment`/`data`/`deviations` cambia il Methods e le Limitations del `report.tex`; con le chiavi `reproduction ...` il `report.tex` e `report_fallback.md` risultano identici): `grnbench report EXP-20260925-LM-2` legge il record come un dizionario piatto di righe `- chiave: valore` e, a chiavi doppie, vince l'ULTIMA. Una riga `- environment:`, `- data:`, `- deviations:` nella 4b sovrascriverebbe quelle della sezione 2 e finirebbe nel Methods / Configuration / Limitations del LaTeX di Luca alla prossima rigenerazione. Per questo ogni chiave della 4b ha il prefisso `reproduction`, e gli output di `compare_runs.py` vanno in un blocco recintato (``` ... ```, che il parser ignora). Non aggiungere altre righe `- status:`, `- machine:`, `- type:`, `- seeds:`, `- code:`, `- primary comparison:` nella 4b.
Schema:
```
### 4b. Second-person reproduction (GB, <YYYY-MM-DD>)

- directed and checked by: GB (@Gabriele-Bambini), signs this section. Every step decided and checked by GB from the record and issue #24. Commands typed by: a Claude Code agent session run for GB. Data copied/fetched by GB and checked on GB's own cluster account.
- reproduction job: issue #24. Slurm job <JOB>, <State> <ExitCode>, elapsed <Elapsed> (from sacct), 3 dataset workers, QOS stud, Bocconi HPC (Linux; no node name).
- code under test: 1d0722179f9109ddcc93eb5665b2deefd52d4106 (tree identical to the PR #7 head 0bd3e99; the original run used 67ca880). CODE_COMMIT of the rerun: <copiato>.
- reproduction tooling: 245ac12729c162cb88023256243ce78abe5ee7c8 (deviation: current main's sbatch passes --max-attempts, which 1d07221 does not accept; the rest of the tooling used here is identical to main), compare_runs.py from PR #43 head <sha> (not on main at the time).
- reproduction environment: conda env grnbench from requirements/lock-py311.txt (identical between 1d07221 and main), Python <versione>, pip --dry-run showed <nessuna differenza>.
- reproduction data: references = registered CollecTRI files copied from the login node (copy provided by LM, #24 comment of 2026-10-03), sha256 checked = 4e5597b1... (human), 514abd67... (mouse), equal to experiments/datasets.tsv; CELLxGENE: fetch_cellxgene.py --from-manifest experiments/EXP-20260925-LM-2.datasets.tsv; cisTarget: hpc/fetch_cistarget.sh; data_manifest.py --check experiments/datasets.tsv: <ultima riga>. Live-listing dry-run vs the 20 registered ids: <SAME / differenze>.
- reproduction plan: 20 tasks, ids and ground truth identical to the committed tasks.json.
- reproduction run: 20/20 rows, 0 errors, models pearson_signed+pearson_motif+motif_prior, KEEP_ARTIFACTS=1, TMPDIR on the home filesystem, report with --no-final-report (figures/LaTeX not regenerated; originals in benchmarks/results/EXP-20260925-LM-2/).
- reproduction comparison (compare_runs.py, outputs verbatim):
  ```
  <incolla compare_exact.txt>
  <incolla compare_rtol1e-9.txt>
  ```
- reproduction verdict: <MATCH / DIFFERENT> by the rule fixed in #24 (--rtol 1e-9 primary line). C4: <reproduced / reproduced-computational>.
- reproduction deviations: sbatch/tooling commit above; two JSON names instead of compare.json; <ogni risottomissione o spostamento di righe, con il job id>; the section 1 text (human 16 / mouse 4, 8 protocols) disagrees with the manifest (13 / 7, 7 protocols): noted, section 1 not edited.
- what this adds: another OS and CPU than the original Windows laptop (and a human deciding each step). What it does not test: <altro hardware, altro BLAS oltre a questo cluster>.
- reproduction large outputs (edges, processed data): ~/gb_runs/EXP-20260925-LM-2-repro-GB on the cluster, outside git.
```
Se la regola e' soddisfatta, nello stesso PR (e solo allora):
- le due righe di testa della sezione 4 che oggi dicono `reproduced-computational` e "human reproduction pending": `- status:` (conserva la storia: `reproduced`, data, "second person: GB, issue #24"; la rerun computazionale del 2026-09-26 resta citata) e `- reproduced by:` (e' la riga che `grnbench report` stampa come "Reproduced by": se la lasci com'e', il LaTeX continua a dire "human reproduction pending");
- la riga dell'indice in `experiments/README.md`, la riga C4 di `paper/CLAIMS.md` (status e "reproduced by"), `.projectos/STATE.md` (passo 1 completato), una riga in `## Unreleased` di `CHANGELOG.md`, `paper/CONTRIBUTIONS.md` (CRediT: Validation);
- NON rigenerare `benchmarks/results/EXP-20260925-LM-2/latex/` (e' la cartella di Luca, `experiments/README.md` regola 5): il suo `report.tex` resta con il vecchio stato finche' Luca non lo rigenera con `grnbench report EXP-20260925-LM-2`; scrivilo nel PR.
Se il confronto NON e' `MATCH`: della 4b basta la parte sopra (i fatti e gli output verbatim); `status` e `reproduced by` di testa non cambiano.
Il PR: compila `.github/pull_request_template.md` (what, why, experiment IDs, test con output reale, rischi, handoff; AGENTS.md sezione 5 e 12); la descrizione e' l'handoff. Il PR tocca `paper/`: richiede la review di LM (colonna "reviews" di AGENTS.md sezione 2). Lo stato `reproduced` lo da' la tua rerun diretta e firmata, non una review di agente (AGENTS.md sezione 5): la review di LM controlla il PR e i numeri.
Nota per il PR: `experiments/README.md` regola 5 dice "a rerun is a new experiment ID that links the old one", mentre #24, la regola 8 e `experiments/TEMPLATE.md` (sezione 4 "Reproduction") mettono la riproduzione di seconda persona nella sezione 4 del record esistente, come fatto per #8 e #15; qui si segue #24 (decisione di Luca) e si segnala la tensione fra le due regole a Luca.
I `compare_*.json`, i log e `tasks.json` non vanno in git. **STOP se**: stai per modificare la sezione 1, scrivere in `benchmarks/results/EXP-20260925-LM-2/` o impostare `reproduced` senza il `MATCH` del passo 24 o senza la review di LM.

## 27. Chiusura e pulizia

**Comando**:
```bash
gh issue edit 24 --remove-label running --add-label done
gh issue close 24 --comment "result: <MATCH|DIFFERENT> exact / <MATCH|DIFFERENT> rtol 1e-9, n=20, C4 -> <status>; details in record section 4b"
```
(`done` quando la run e il confronto sono stati fatti, qualunque sia il verdetto; `failed` se non si e' arrivati a un confronto: mia lettura di AGENTS.md sezione 6, UNVERIFIED rispetto a Luca.)
**STOP se**: non hai ancora i due output del passo 24 salvati e incollati in 4b (la chiusura riporta solo cio' che e' stato registrato); il job risulta ancora in coda o in esecuzione (`bash hpc/hpc.sh status` non vuoto: non chiudere, la label `running` resta); `gh` risponde con un errore di permessi (non cambiare account, chiedi a Luca).
Pulizia sul cluster: `bash hpc/hpc.sh run 'rm -rf ~/gb_tmp/lm2_listing'`. NON cancellare `~/grnbench_data_lm2repro` (le copie registrate dei riferimenti servono a #52 e a un'eventuale rerun) ne' `~/gb_runs/EXP-20260925-LM-2-repro-GB` (percorso citato nel record). Lo scratch dei worker lo rimuove lo sbatch.

---

# Se qualcosa va storto (tabella rapida)

| sintomo | causa probabile | cosa fare |
|---|---|---|
| `hpc.sh run`: `cd: ... No such file` | sync non fatto | passo 6 |
| job `FAILED` in pochi secondi, log con `unrecognized arguments: --max-attempts` | usato lo sbatch di `main` attuale | F1: strumenti a `245ac12`; ricontrolla passo 10; risottometti (non ha scritto nulla) |
| job `FAILED` in pochi secondi per virgolette | `$HOME` espanso in locale | rivedi gli apici singoli del passo 17 (precedente: EXP-20260926-LM-1, job fallito in 9 s) |
| `OSError: No space left on device` | `TMPDIR` non arrivato o home piena | verifica `TMPDIR=$HOME/gb_tmp` nell'`--export`; `df -h ~`; precedente LM-1 deviazione 1 |
| `QOSMaxJobsPerUserLimit` | hai gia' un job | non sottomettere altro |
| `Could not resolve hostname` | VPN | comando stampato da `hpc.sh`, una volta |
| `Permission denied` / `Not allowed at this time` | chiave / lockout | regola 1: STOP, 10 minuti, nessun ciclo |
| `ok-eol` per `truth.json` | solo fine riga | accettato (passo 14) |
| `MISMATCH` | file diverso dal registrato | STOP e riporta; mai sostituire righe del registro |

# Non verificato (UNVERIFIED), in un posto solo

- Qualunque cosa sul cluster: raggiungibilita' di CELLxGENE e aertslab dal nodo di login, presenza di `/tmp/grnbench_registered_references_LM` sul TUO nodo, stato dell'env `grnbench` (e se coincide col lock), spazio/quota in home, stato della coda e dei limiti QOS di oggi.
- (misurato il 2026-10-08) Il tempo del job sul cluster con 3 worker: 35 min (la stima era 1-2 h). Artefatti: 177 MB.
- Il comportamento di ssh con `nohup ... &` (passo 11), `pip install --dry-run` sulla tua versione di pip (passo 9), `sbatch --test-only` (passo 16).
- L'esito del listing live di CELLxGENE rispetto ai 20 id (passo 12) e la mia lettura che una differenza non blocchi.
- (era UNVERIFIED, ora provato nel sandbox Linux) I numeri attesi di `compare_runs.py` (1180/40/400/140, `3.39e-08`, `MATCH`): rieseguiti sui risultati committati con la testa `2cceed0` e il venv dal lock. Resta non provato che sul TUO cluster la rerun combaci: e' proprio il risultato dell'esperimento.
- Lo stato futuro della PR #43 (testa `2cceed0`, aperta al 2026-10-06 21:42Z).
- Che `done` vs `failed` valga come scritto al passo 27.
