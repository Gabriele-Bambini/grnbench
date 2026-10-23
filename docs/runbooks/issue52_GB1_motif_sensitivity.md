# Runbook: EXP-20261003-GB-1 (issue #52), scomposizione del lift del motif prior, stage 1 (Pearson)

> **Stato al 2026-10-10: lanciato** (job 692696, 15:25Z, da `main` `42c6fb2` dopo lo squash-merge di #55; passi 4-13
> eseguiti come scritti, senza STOP: hash della sezione 1 `7823bf8a...`, diff di identita' solo `autoprep.py | 1 +`,
> input e riferimenti OK, `ok=8 ok-eol=0 MISMATCH=0 UNLISTED=0 MISSING=0`, piano 40 unita', 0 escluse).
> La strada del codice e' la A del passo 4 (merge prima del lancio). Le note qui sotto sono quelle scritte prima del lancio.
> Risultati (aggiunto il 2026-10-23): PR #92 (C6 inconclusive, C4 mixed, combinato inconclusive).
> - Passi 7-8 fatti: artifact di LM copiati in `~/grnbench_artifacts/` con `INPUT_SHA256SUMS` costruito dalla sorgente e
>   verificato sulla copia (EXP-20261002-LM-1 424 file, EXP-20260925-LM-2 315, tutti OK).
> - Passi 9-10 fatti: riferimenti registrati (sha256 OK) e cisTarget in `~/grnbench_data_gb1`. I sei file cisTarget sono
>   stati copiati da `~/grnbench_data_lm2repro/cistarget`, gia' verificato per #24, invece di riscaricarli.
>   `data_manifest.py ~/grnbench_data_gb1 --check experiments/datasets.tsv`: `ok=8 ok-eol=0 MISMATCH=0 UNLISTED=0`.
>   Al lancio va comunque rifatto il passo 10 con lo snapshot e `--require-all`.
> - Passo 3, risolto il 2026-10-10: LM ha risposto (c) su #55 (non sara' il secondo umano prima di lunedi'), quindi vale la
>   strada di AGENTS.md sezione 5. La review di agente di `661514b` resta, con la nota 7 che LM ha chiesto. Le note 7 e 8 sono
>   nel record (sezione 2.2, testa `b3ec128`).
>
> **Gate 1 e piattaforma (aggiornato il 2026-10-10).** Prima LM aveva scritto (#82, review di #83) che il gate 1 sarebbe
> probabilmente fallito su Linux per gli artifact di C4, fatti su Windows. LM l'ha poi corretto su #55. La matrice processata
> di `prostate-gland__20de8570` fatta su Linux, valutata sul portatile Windows, da' la riga Linux fino all'ultimo bit in
> entrambi i bracci (differenza 0). Una matrice identica si valuta in modo identico sulle due piattaforme, e le matrici
> tenute sono byte fissi: per C4 il gate 1 dovrebbe passare sul cluster (UNVERIFIED finche' il job non gira). La
> differenza di #24/#82 nasce prima, nella normalizzazione rifatta su un'altra piattaforma, che il gate 1 non ripete.
> Se comunque alcuni dataset di C4 falliscono il gate 1, LM li rilancia sul portatile di LM-2 (nota 8, punto 5 della nota 2.1).
> Un gate 1 fallito e' un risultato da riportare, non da aggirare cambiando la tolleranza. #86 e #88 non devono essere nel
> commit del job (nota 8).
>
> **Fatto al lancio del 2026-10-10 (scritto prima del lancio, da una copia pulita del repo su `main` = `42c6fb2`):** passi 1, 2, 4 (strada A),
> 5, 6, 10, 11, 12, 13; poi 14-20. I passi 3, 7, 8 e 9 sono fatti (sopra). Tutti i controlli del passo 4 danno l'atteso su `42c6fb2`
> (il punto 5 e' ora `118 0`, vedi il passo 4). Per non dipendere da cio' che finisce su `main` prima del lancio, lancia da quel commit:
> `git switch --detach 42c6fb2` al posto di `git pull`. `hpc.sh sync` registra `git rev-parse HEAD` (e `-dirty` se l'albero non e pulito),
> quindi funziona anche da HEAD staccato.

Letto sullo stato di GitHub e del repo al 2026-10-07 (PR #55 testa `661514b3713b4937d1bef0aee599591dd2962d36`, `main` = `cfeee83`).
Autorita per i comandi: sezione 2 del record `experiments/EXP-20261003-GB-1.md` (sulla branch `gb/motif-sensitivity`).
Dove questo runbook aggiunge un controllo o cambia la forma di un comando lo dice; la sostanza del comando non cambia.
Riletto da un secondo agente (scettico) contro record, script `motif_sensitivity.py`, `hpc.sh`, `engine_bench.sbatch`, issue #24/#41/#52 e PR #55: i comandi del record sono riportati senza cambiare i flag
(controllato per confronto automatico delle righe); le correzioni di quella rilettura sono le discrepanze 10-13 e le note sui passi 1, 2, 3, 4, 5, 6, 8, 14, 16 e 17.

## 0. Prima di cominciare

**Cosa e un job unico.** Un solo job `stud` da 8 CPU e 32 GB (`--time=04:00:00`) ricalcola Pearson dai file tenuti da LM per 40 dataset
(C6 = `EXP-20261002-LM-1`, primario; C4 = `EXP-20260925-LM-2`, secondario) e non rifitta nulla. Stima del record: 1.5-2.5 CPU-ore, quindi meno di
un'ora di tempo reale una volta partito (mai misurato su dati reali). L'attesa in coda puo essere di ore.
Dopo il job, un secondo job leggero (1 CPU, 20 min) scrive il report. Mai due job insieme: QOS `stud` = 1 in esecuzione + 1 in coda.

**Chi lo esegue.** Tu, dalla tua macchina, VPN Bocconi attiva, nella root del repo, in bash (Git Bash su Windows).
Un agente puo digitare i comandi al posto tuo ma non tocca mai la passphrase e non sottomette il job (passo 12) senza che tu abbia superato i passi 1-11.

**Segreti.** ID universitario, passphrase e password non compaiono mai nel repo, in `.env.hpc`, nelle issue, nelle PR, nei commenti, nei log incollati
ne in questo file. L'ID sta solo nel campo `User` del tuo `~/.ssh/config`; la passphrase solo in ssh-agent (`ssh-add`, la digiti tu).
`hpc.sh` lavora solo con chiave + agent
(`BatchMode=yes`): una password di account non puo funzionare tramite `hpc.sh`.
Ogni output che incolli su GitHub: sostituisci prima qualsiasi percorso assoluto di home con `~` e non incollare mai l'output di `hpc.sh check`, `sync`,
`log` (la prima riga ha il nome del nodo di calcolo) e `status` (colonna dei nodi), ne quello di `du -sh <percorso>`, `ls -l` e `df` (percorso di home, proprietario).
Nei commenti scrivi solo le righe riassuntive date in questo file.

**Variabili** (da ripetere in ogni terminale nuovo; `L` e `R` servono in locale perche i comandi del record le espandono prima di spedirle al cluster):

```bash
cd <root del tuo clone di grnbench>
EXP=EXP-20261003-GB-1
A='source /software/miniconda3/etc/profile.d/conda.sh && conda activate grnbench'
L=/tmp/grnbench_artifacts_LM                 # cartella condivisa da LM sul nodo di login (issue #41)
R=/tmp/grnbench_registered_references_LM     # riferimenti registrati (issue #24)
FETCH=~/grnbench_fetch                       # cartella locale FUORI dal repo
D=benchmarks/results/$EXP                    # cartella dei risultati nel repo (passi 17-20)
: "${EXP:?}" "${A:?}" "${L:?}" "${R:?}" "${FETCH:?}" "${D:?}"   # nessuna variabile deve essere vuota
```

`JOB` e `JOB2` (passi 12 e 16) vivono solo nel terminale che li ha creati: se lo chiudi, rileggili dal commento di lancio su #52 e riassegnali a mano
(`JOB=<numero>`). Se il record usa `R` per "la tua copia dei riferimenti di #24": qui `R` e la cartella di LM (passo 7) perche la tua copia non esiste ancora;
se l'hai gia, puoi puntare `R` a quella, ma i due hash del passo 7 restano il controllo.

**Regole di stop valide per ogni comando `hpc.sh`** (non ripetute passo per passo):

| messaggio | cosa significa | cosa fai |
|---|---|---|
| `Could not resolve hostname` / timeout | VPN giu | esegui il comando VPN che `hpc.sh` stampa come prima riga, poi riprova UNA volta; se fallisce ancora, fermati e segnala |
| `Permission denied (publickey)` | chiave non nell'agent o `User` sbagliato | fermati; `ssh-add` lo fai tu; controlla `User` in `~/.ssh/config`; riprova UNA volta; mai in loop |
| `Not allowed at this time` | blocco dopo login falliti | fermati, aspetta 10 minuti, correggi la causa, riprova UNA volta; mai in loop |
| `cp: ... Permission denied` / `cannot open` letto da remoto | e un permesso su un file, non l'autenticazione, anche se `hpc.sh` stampa comunque il suggerimento "Auth failed" | fermati e segnala; non riprovare in loop; una sola `bash hpc/hpc.sh check` per escludere l'auth |
| qualunque altro errore non previsto dal passo | non indovinare | fermati e segnala su #52 (senza percorsi personali) |

Nessun calcolo sul nodo di login: ci girano solo copie, checksum, download, `sbatch`/`squeue`/`sacct`, il controllo `data_manifest.py` e il `plan`
(Python leggero, come prescrive il record; il record chiama "login node: copies and checksums only" anche questi due, UNVERIFIED se la policy del nodo li tollera).
Sul cluster non si cancella nulla: l'unica eccezione e un file `.part` di cisTarget lasciato da un download interrotto (passo 9). Mai una riga di `rows/` per far passare un report.

---

## PARTE A. Sulla tua macchina

### 1. Sei la persona giusta e nessun segreto e finito nell'albero del repo

```bash
gh issue view 52 --json assignees,labels,state --jq '{assegnatari:[.assignees[].login], label:[.labels[].name], stato:.state}'
git check-ignore .env.hpc
grep -c '^GRNBENCH_HPC_HOST=.' .env.hpc          # stampa solo il conteggio, non il valore
git status --short
ssh-add -l                                       # Windows/Git Bash: /c/Windows/System32/OpenSSH/ssh-add.exe -l
```

Atteso: assegnatario `Gabriele-Bambini`, label `["experiment"]` (senza `running`), stato `OPEN`; `.env.hpc` stampato dal secondo comando (e ignorato);
`1` dal terzo (il file esiste e ha l'alias; si crea con `cp .env.hpc.example .env.hpc` e si imposta `GRNBENCH_HPC_HOST` al TUO alias di `~/.ssh/config`, docs/HPC.md sez. 1);
`git status --short` senza output; `ssh-add -l` elenca almeno una chiave (e il solo controllo che la tua chiave e nell'agent: la passphrase l'hai digitata tu, mai un agente).
La VPN e il login vero si provano al passo 5 con `hpc.sh check`, prima di `sync`.

STOP e segnala se: non sei GB (la #52 e di GB: chiedi a GB/LM di riassegnarla o apri una issue di riproduzione `<EXP-ID> reproduction`, AGENTS.md sez. 6; UNVERIFIED chi sei);
`.env.hpc` non e ignorato o manca; `git status --short` elenca file (soprattutto uno con credenziali: spostalo fuori dall'albero, `sync` copia tutto l'albero sul cluster);
`ssh-add -l` dice `The agent has no identities` / `Could not open a connection to your authentication agent`: carica la chiave TU (`ssh-add ~/.ssh/<tua-chiave>`, la passphrase la digiti tu nel tuo terminale), poi riprova il controllo una volta.

### 2. Repo aggiornato e bacheca: nessuno sta gia correndo su questo

```bash
git fetch origin
PYTHONPATH=src python -m grnbench.board          # senza `pip install -e .` serve PYTHONPATH=src (AGENTS.md sez. 4); con il venv dal lock attivo
gh issue list --label running --state open
```

Atteso: "nothing running" (alla lettura del 2026-10-07: zero issue `running` aperte) e #52 tra gli esperimenti rivendicati.
Rileggi comunque `AGENTS.md` (sezioni 5-7, 10, 11) e le voci ACTIVE di `.projectos/DECISIONS.md` (D-020: riferimenti e cisTarget con hash registrati).

STOP se: una issue `running` cita `EXP-20261003-GB-1`, `~/gb_runs/EXP-20261003-GB-1` o `~/grnbench_artifacts`: non rilanciare, commenta quella issue.
La bacheca ha una riga che comincia con `- ?` (per esempio `gh not usable` / `gh failed`): non ha letto GitHub, quindi "nothing running" non prova nulla.
`gh auth login` lo fai tu, poi ripeti il passo; senza `gh` funzionante non c'e modo di sapere cosa gira e non si va avanti (la pagina web delle issue con label `running` e un controllo manuale accettabile).

### 3. Passo UMANO: la review di LM su #55 (condizione 2 della issue #52). Non si salta.

```bash
gh pr view 55 --json state,headRefOid,isDraft,reviewDecision,mergeStateStatus
gh pr view 55 --json reviews --jq '.reviews[] | [.author.login, .state, .commit.oid, .submittedAt] | @tsv'
```

Cosa c'e oggi su #55 (letto da GitHub): PR aperta, testa `661514b`. Sotto l'account di LM (`elle13eth`) ci sono due review: `COMMENTED` su `94dc34f`
e `APPROVED` su `661514b`, ma il testo di entrambe dice di essere una review di agente indipendente "for LM", "not a second-human review".
La descrizione della PR ha ancora la casella "review: LM ... LM's approval pending" non spuntata e `STATE.md` dice "LM approves after reading them".
Quindi: la review umana di LM su nota 2.1 e codice statistico non risulta ancora fatta (UNVERIFIED se LM l'ha letta fuori da GitHub).

Cosa fai tu, e solo tu: chiedi a LM una risposta scritta su #55. Testo suggerito (il repo e in inglese):

> @elle13eth: before I launch EXP-20261003-GB-1 (#52) I need your own decision on #55 as statistical reviewer. The two reviews under your account (`94dc34f`, `661514b`)
> are worded as independent agent reviews, not second-human reviews. Please either (a) read record note 2.1 and the statistical code (`benchmarks/motif_sensitivity.py`) yourself
> and approve, or (b) state explicitly that you have read them yourself and adopt the `661514b` review as your own, or (c) tell me you are unavailable.
> I am not launching until one of these is on the PR.

Atteso prima di andare avanti: una di (a), (b), (c) scritta da LM su #55, e `headRefOid` ancora `661514b3713b4937d1bef0aee599591dd2962d36`
(se la testa e cambiata, l'approvazione non copre il codice che lancerai: richiedila). Un "ok" generico, o la sola approvazione di agente gia presente, non basta:
(b) vale solo se LM dichiara di aver letto lui.
Se LM e indisponibile, AGENTS.md sez. 5 prevede un'altra strada (merge dopo review di agente indipendente, dichiarata come tale nella review e nel messaggio di squash,
e LM puo contestarla dopo con una issue). Le condizioni sono due e non le decidi tu a posteriori: LM ha detto di essere indisponibile (caso (c)), oppure non ha risposto
in 2 giorni lavorativi **alla richiesta qui sopra** (il conto parte da quando la posti, non dalle review del 2026-10-06: alla data in cui questo file e stato scritto, 2026-10-07, non era ancora scaduto nulla).
Se la prendi, scrivilo su #52 e #55 con il motivo. La review di agente della sez. 5 e una sessione nuova che non ha scritto la PR, con suite completa in venv dal lock
e output vero. L'approvazione di `661514b` dichiara `569 passed, 1 skipped` ("a fresh session"; UNVERIFIED chi l'ha avviata e se e quella della sez. 5): non e una review
umana, e puo valere come tale solo se una delle due condizioni sopra e vera e il messaggio di squash lo dice ("independent agent review, not a second-human review").

STOP: nessuna risposta scritta di LM, e nessuna delle due condizioni della strada sez. 5 soddisfatta e dichiarata per scritto = non si va avanti.
Non fare merge e non lanciare sull'onda della sola approvazione di agente.

### 4. Quale codice va sul cluster, e il record e quello pre-registrato

Due strade, entrambe ammesse; scegli una.

- **A (consigliata, dopo il passo 3): fai il merge squash di #55**, poi `git switch main && git pull --ff-only`.
  Il job gira sul commit di `main` che `sync` registra. I risultati si committano poi da `main`, senza PR impilate.
- **B: lanci dalla branch prima del merge.**
  ```bash
  git switch gb/motif-sensitivity && git pull --ff-only
  git rev-parse HEAD        # atteso 661514b3713b4937d1bef0aee599591dd2962d36, cioe la testa che la review copre
  ```
  Poi, per i risultati, aspetta il merge e parti da `main` aggiornato (niente PR impilate). Se fai il merge squash dopo il lancio:
  NON eseguire `hpc.sh sync` tra il job e un'eventuale risottomissione (vedi passo 15).

Controlli (in entrambe le strade, dalla root del repo):

```bash
git status --short                                                    # nessun output
ls experiments/$EXP.md experiments/$EXP.snapshot.tsv benchmarks/motif_sensitivity.py
git fetch origin pull/55/head                                         # porta in locale 9d3f049 anche dopo uno squash
git show HEAD:experiments/$EXP.md | awk '/^## 1\. Registration/{f=1} /^## 2\. How it runs/{f=0} f' | sha256sum
git diff --numstat 9d3f049 HEAD -- experiments/$EXP.md
git diff --stat 6ad729c HEAD -- benchmarks/real_bench.py src/grnbench/infer src/grnbench/eval src/grnbench/tfs.py src/grnbench/autoprep.py
for S in EXP-20261002-LM-1 EXP-20260925-LM-2; do ls benchmarks/results/$S/rows | wc -l; done
grep -c -F '| pearson_motif | pearson_signed | 20 | +0.203 | [+0.147, +0.259] | 1.62e-06 | 1 | pearson_motif |' benchmarks/results/EXP-20261002-LM-1/REPORT_REAL.md
grep -c -F '| pearson_motif | pearson_signed | 20 | +0.182 | [+0.138, +0.227] | 1.36e-07 | 1 | pearson_motif |' benchmarks/results/EXP-20260925-LM-2/REPORT_REAL.md
```

Atteso, riga per riga:
1. nessun output;
2. i tre file esistono;
3. (fetch senza errori);
4. `7823bf8ab0beec137fd8d66a5b9cb73b71b9c4e381af35f385f48a874a767608  -` (sezione 1 byte-identica alla pre-registrazione `9d3f049`, 19573 byte; verificato da me su `9d3f049` e `661514b`;
   si usa `git show` e non il file in chiaro perche su Windows `* text=auto` puo cambiare i fine riga e quindi l'hash);
5. `118	0	experiments/EXP-20261003-GB-1.md` su `main` `42c6fb2` (solo righe aggiunte: la nota 2.1 e la sezione 2.2 con le note 7 e 8; erano `86` a `661514b`, prima della sezione 2.2; informativo: l'hash del punto 4 e il controllo duro);
6. esattamente `src/grnbench/autoprep.py | 1 +` e `1 file changed, 1 insertion(+)` (il codice di punteggio e quello dei run C4/C6; la riga di `autoprep.py` e dichiarata nella nota 2.1);
7. `20` e `20`;
8. `1` e `1` (le due righe di intestazione che il gate 2 deve riprodurre).

Solo se sei su Windows: `git config --get core.longpaths` deve stampare `true` (altrimenti `git config core.longpaths true` prima del pull: la cartella
`benchmarks/results/EXP-20260925-LM-2/viz/` fa fallire il clone con "Filename too long").

STOP se: l'hash della sezione 1 e diverso (la sezione 1 e la pre-registrazione: un cambio e una nuova registrazione, non un ritocco); la diff del punto 6 mostra altri file (il codice che
calibra non e piu quello dei run: il gate 1 potrebbe fallire per ragioni non statistiche); i conteggi non sono 20 e 20; un grep da `0`.

Su `origin/main` letto il 2026-10-07 (`cfeee83`) e il 2026-10-10 (`42c6fb2`, dopo il merge di #55) la diff del punto 6 e esattamente quella attesa. Ma #40 e #69 (PR aperte alla lettura del
2026-10-07) toccano `src/grnbench/infer/classical.py` e `base.py`, cioe il codice di Pearson, e #86 e #88 (LM le tiene fino a dopo il job, nota 8) non devono entrare nel commit del job:
se una di queste viene unita su `main` prima del lancio, nella strada A il punto 6 elenchera altri file. Non e un rifiuto del codice di #55: lancia dal commit `42c6fb2`
(`git switch --detach 42c6fb2`, poi `git rev-parse HEAD` deve stampare `42c6fb211c70fe775512295ca68d2176b2fbd16a`) oppure aspetta il parere di LM.
(`git fetch origin pull/55/head` porta ora la testa finale della PR, `b3ec128`, non piu `661514b`; ha lo stesso albero di `42c6fb2`.)
Non giudicare tu se quelle modifiche cambiano i bit dei punteggi: lo dice solo il gate 1.

### 5. Sincronizza il codice (una volta sola)

Ordine del record: `git status --short`, poi `hpc.sh check`, poi `hpc.sh sync` (prima si prova il login, poi si spedisce l'albero).

```bash
git status --short                      # ancora vuoto
bash hpc/hpc.sh check
bash hpc/hpc.sh sync
git rev-parse HEAD
bash hpc/hpc.sh run 'cat GIT_COMMIT'
```

Atteso: `check` stampa `connected: ...`, la versione di Slurm e la cartella `~/grnbench` (o `no ~/grnbench yet: run sync` al primo uso; il nome della cartella e `GRNBENCH_HPC_DIR` di `.env.hpc`);
`sync` stampa `synced <40 caratteri esadecimali> to ...`; l'hash e uguale a `git rev-parse HEAD`; `cat GIT_COMMIT` lo ripete, senza `-dirty`.
Annotalo: e il commit che andra nel commento della issue e in `CODE_COMMIT`. Non incollare le righe di `check` e `sync` (contengono utente e nome del nodo): incolla solo l'hash.

STOP se: l'hash finisce con `-dirty` (risultati da un commit sporco non sono accettati dal record: torna al passo 4); `check` o `sync` falliscono (regole di stop sopra: VPN, `Permission denied`, `Not allowed at this time`).
Da qui in poi non fare un secondo `sync` fino alla chiusura della issue.

### 6. Coda libera, ambiente, spazio

```bash
bash hpc/hpc.sh status
bash hpc/hpc.sh run "$A && python --version"
bash hpc/hpc.sh run 'df -h ~ | tail -1'
```

Atteso: `status` mostra solo l'intestazione di `squeue` (nessun tuo job); Python `3.11.x`; `df` mostra molti GB liberi
(servono circa 1.3 GB di cisTarget + 0.4 GB di input + pochi MB di output + temporanei del job; la quota di home e UNVERIFIED, `df` mostra il disco, non la quota).
Niente import pesanti di prova sul nodo di login: se mancasse un pacchetto lo mostra il passo 11 (`plan`), che li importa comunque, prima di spendere il job.

STOP se: hai gia un job in coda o in corso (il limite `stud` e 1 + 1: aspetta che finisca, non sommare); l'ambiente `grnbench` manca
(`cannot activate`/`No module`): `bash hpc/hpc.sh run 'bash hpc/setup_env.sh'` (installazione lunga, lo permette la skill) e ripeti il controllo; spazio sotto ~5 GB.
Solo se questo account non ha mai lanciato un job: prima prova da 1 minuto
`bash hpc/hpc.sh submit "--partition=stud --qos=stud --time=2 --cpus-per-task=1 --mem=1G -o logs/hello_%j.out --wrap='hostname; nproc'"`, aspetta che finisca, poi continua.

---

## PARTE B. Input sul cluster (nodo di login: solo copie e checksum)

### 7. Controlla le sorgenti prima di copiare

`/tmp` del nodo di login e locale al nodo e puo essere ripulito. `hpc.sh run` deve atterrare sul nodo dove LM ha messo i file (nome in #41); se il tuo alias
ti porta su un altro nodo, le cartelle non ci saranno (UNVERIFIED a quale nodo arrivi).

```bash
bash hpc/hpc.sh run "L=$L; R=$R; ls -ld \$L/EXP-20261002-LM-1 \$L/EXP-20260925-LM-2 \$R >/dev/null && for X in EXP-20261002-LM-1 EXP-20260925-LM-2; do echo \$X rows=\$(find \$L/\$X/rows -type f | wc -l) task_dirs=\$(find \$L/\$X/artifacts -mindepth 1 -maxdepth 1 -type d | wc -l) processed=\$(find \$L/\$X/artifacts -name processed.h5ad | wc -l) size=\$(du -sh \$L/\$X | cut -f1); done; sha256sum \$R/collectri_human.tsv \$R/collectri_mouse.tsv"
bash hpc/hpc.sh run 'test -e ~/grnbench_artifacts && echo "ESISTE GIA" || echo "assente: ok"'
```

Atteso: per ciascuno dei due set `rows=20 task_dirs=20`, `processed` circa 90 in tutto tra i due (il record parla di "90 contexts in all"; almeno 80, due bracci per dataset) e una dimensione dell'ordine di 170-220 MB
(record: C6 circa 214 MB, C4 circa 179 MB). Le due righe di `sha256sum` devono essere esattamente

```
4e5597b16a0e795804e1ee0d47bbb420e144e9c0878f8b2ea05cac34be7e303f  .../collectri_human.tsv
514abd67c4c3b02c896e186cefe1c48070e77c3bafe71df871d0c9e86ba8dab4  .../collectri_mouse.tsv
```

(sono gli hash registrati in `experiments/datasets.tsv`, ripetuti da LM in #24). Il secondo comando deve dire `assente: ok`.
Non copiare `EXP-20260926-LM-1` (C5): e stage 2, fuori ambito.

STOP se: `No such file or directory` (cartella ripulita o nodo diverso): chiedi a LM di rimetterla sul tuo nodo, non ricostruire gli artifact e non rilanciare C4/C6; permesso negato su un file: LM deve rendere
leggibile (`chmod -R a+rX`), non riprovare in loop; conteggi diversi da 20/20; hash dei riferimenti diversi (non sono i byte registrati: D-020 vieta di sostituirli, chiedi a LM);
`ESISTE GIA`: una copia precedente, parziale o no, non si sovrascrive alla cieca: guarda cosa contiene e segnala. Eccezione: se e il resto di un passo 8 TUO
interrotto (connessione caduta a meta copia) salta questo controllo e rifai il passo 8 (e idempotente, vedi li).

### 8. Copia gli artifact di LM nel tuo account, con sha256 (comando del record, alla lettera)

```bash
: "${L:?}" "${R:?}"            # guardia locale: L e R devono essere valorizzate
for X in EXP-20261002-LM-1 EXP-20260925-LM-2; do
  bash hpc/hpc.sh run "mkdir -p ~/grnbench_artifacts/$X && (cd $L/$X && find rows artifacts -type f -print0 | sort -z | xargs -0 sha256sum) > ~/grnbench_artifacts/$X/INPUT_SHA256SUMS && cp -r $L/$X/rows $L/$X/artifacts ~/grnbench_artifacts/$X/ && cd ~/grnbench_artifacts/$X && sha256sum -c --quiet INPUT_SHA256SUMS && echo $X inputs OK"
done
bash hpc/hpc.sh run 'cd ~/grnbench_artifacts && for X in EXP-20261002-LM-1 EXP-20260925-LM-2; do echo $X listed=$(wc -l < $X/INPUT_SHA256SUMS) on_disk=$(find $X/rows $X/artifacts -type f | wc -l) size=$(du -sh $X | cut -f1); done'
```

Atteso: due righe `EXP-20261002-LM-1 inputs OK` e `EXP-20260925-LM-2 inputs OK`, nessuna riga `FAILED`; poi, per ciascun set, `listed` uguale a `on_disk` e la stessa dimensione del passo 7.
(`INPUT_SHA256SUMS` e calcolato sulla copia di LM e verificato sulla tua: prova che la copia e fedele, non che i file di LM siano quelli originali; quella prova e il gate 1 nel report.)

STOP se: compare `FAILED`, `WARNING: ... did NOT match`, `No such file` o `listed` != `on_disk`: non usare la copia, segnala; se `L` o `R` erano vuote la guardia ferma il passo prima di spedire nulla.
Connessione caduta durante la copia (ssh esce con errore, nessuna riga `inputs OK` per quel set): rilancia SOLO questo passo 8. Il comando ricalcola `INPUT_SHA256SUMS` dalla sorgente, ricopia
(`cp -r` sovrascrive file uguali) e riverifica; nessun altro passo ha ancora letto quella cartella. Se `FAILED` resta anche dopo il secondo giro, segnala: non cancellare nulla.

### 9. Riferimenti registrati e cisTarget nel data root dedicato `~/grnbench_data_gb1`

```bash
bash hpc/hpc.sh run "mkdir -p ~/grnbench_data_gb1/references && cp $R/collectri_human.tsv $R/collectri_mouse.tsv ~/grnbench_data_gb1/references/"
bash hpc/hpc.sh run 'mkdir -p logs && nohup bash hpc/fetch_cistarget.sh ~/grnbench_data_gb1/cistarget > logs/fetch_cistarget_gb1.log 2>&1 &'
bash hpc/hpc.sh run 'tail -3 logs/fetch_cistarget_gb1.log'        # ripeti fino a "done: 6 files in ..."
```

Il record scrive `bash hpc/fetch_cistarget.sh ~/grnbench_data_gb1/cistarget` in primo piano dentro un solo `hpc.sh run`; qui lo stesso script gira con `nohup` in sottofondo (come la skill per
`fetch_engine_datasets.sh`) perche sono circa 1.3 GB e una connessione caduta lascerebbe un `.part`. Stesso script, stessi file: non e una deviazione dell'esperimento.
Se questo account ha gia in cartella gli stessi sei file cisTarget (per esempio sotto un altro data root), puoi copiarli invece di riscaricarli: decide comunque il passo 10.
Non incollare la riga `done:` del log (contiene il percorso della home).

Atteso: il `tail` finisce con `done: 6 files in ...`; poi
```bash
bash hpc/hpc.sh run 'find ~/grnbench_data_gb1 -type f | wc -l; find ~/grnbench_data_gb1 -name "*.part" | wc -l'
```
stampa `8` e `0`.

STOP se: errore HTTP di `curl` (la rete del nodo di login verso il server cisTarget e UNVERIFIED), disco pieno, o un file `.part` rimasto (cancellalo tu, e solo questo, poi rilancia il fetch: lo script salta i file gia completi).

### 10. Controllo D-020 del data root (deve uscire con 0)

```bash
bash hpc/hpc.sh run "$A && python benchmarks/data_manifest.py ~/grnbench_data_gb1 --check experiments/$EXP.snapshot.tsv --require-all"; echo "exit=$?"
```

Atteso: l'ultima riga di `data_manifest` e `ok=8 ok-eol=0 MISMATCH=0 UNLISTED=0 MISSING=0`, poi `exit=0`.
STOP se: qualunque `MISMATCH`, `UNLISTED` (file in piu, `.part` compresi), `MISSING`, oppure `ok-eol` diverso da 0 (il tool lo accetta come "uguale a meno dei fine riga", ma il report poi stampera
"NOT the registered hash" perche il job hasha i byte: i riferimenti devono essere i byte registrati). Un file che non corrisponde non si ricostruisce con un altro fetch: D-020 vieta di sostituire un hash registrato.

### 11. Piano delle 40 unita

```bash
bash hpc/hpc.sh run "$A && mkdir -p ~/gb_runs/$EXP ~/gb_tmp && PYTHONPATH=src python benchmarks/motif_sensitivity.py plan --run EXP-20261002-LM-1=\$HOME/grnbench_artifacts/EXP-20261002-LM-1 --run EXP-20260925-LM-2=\$HOME/grnbench_artifacts/EXP-20260925-LM-2 --out \$HOME/gb_runs/$EXP/sens_tasks.json"
```

Atteso: `40 units -> .../sens_tasks.json; excluded by plan: 0; with notes: <n>`. Se compaiono righe con `note` (una riga della cartella di run diversa dal commit committed) sono note, non esclusioni: il gate 1 confronta
comunque con le righe committed. Copiale (con `~` al posto della home) nel commento di lancio.

STOP se: `excluded by plan` e maggiore di 0 (artifact mancanti: prima controlla che non sia la tua copia, confrontando i conteggi del passo 8 con la sorgente; se manca alla sorgente,
la regola di esclusione (a) del record lo escluderebbe e lo riporterebbe, ma decidete tu e GB/LM prima di spendere il job; non toccare gli artifact); un `ValueError`
su `rows/` o `REPORT_REAL.md` mancanti (l'albero sincronizzato non contiene i risultati committed: torna ai passi 4-5); qualunque traceback.

---

## PARTE C. Lancio, attesa, risultati

### 12. Sottometti il job unico (comando del record, alla lettera)

```bash
JOB=$(bash hpc/hpc.sh submit "--time=04:00:00 --export=ALL,OUT=\$HOME/gb_runs/$EXP,SCRIPT=benchmarks/motif_sensitivity.py,TASKS=\$HOME/gb_runs/$EXP/sens_tasks.json,REFS=\$HOME/grnbench_data_gb1/references,CHUNK=100000,GRNBENCH_CISTARGET_DIR=\$HOME/grnbench_data_gb1/cistarget,TMPDIR=\$HOME/gb_tmp hpc/engine_bench.sbatch")
echo "JOB=$JOB"
bash hpc/hpc.sh status "$JOB"
```

Atteso: un numero di job, poi `PENDING` o `RUNNING` con il nome `gb-engine`. (`--qos=stud`, 8 CPU, 32 GB vengono dall'intestazione dello sbatch; `CHUNK=100000` = tutte le 40 unita in un solo job.)
Per la stima di avvio: `bash hpc/hpc.sh run 'squeue --me --start'`.

STOP se: `QOSMaxJobsPerUserLimit`/`QOSMaxSubmitJobPerUserLimit` (hai gia un job: non sottomettere altro); sbatch rifiuta qualcosa; il numero di job non compare.
Non usare `--array`, non passare `--n-perm` (il registrato e 100: un altro valore il report lo segnala come DEVIATION), non passare mai `MODELS` con GRNBoost2 o RegDiffusion (lo script rifiuta lo stage 2).

### 13. Etichetta `running` e commento di lancio (subito, stesso minuto)

```bash
gh issue edit 52 --add-label running
gh issue comment 52 --body "launched: HPC job $JOB, QOS stud, 8 CPUs, 32 GB, limit 04:00:00. Code commit <hash del passo 5>, clean tree. Started: <data e ora>. Expected finish: <ora di avvio + circa 1 h (record: 1.5-2.5 CPU-h)>. Outputs: ~/gb_runs/EXP-20261003-GB-1 (run folder), ~/grnbench_artifacts/ (copied inputs). Check: bash hpc/hpc.sh status $JOB. Plan: 40 units, excluded by plan 0."
```

Atteso: label aggiunta (esiste: `running`, `done`, `failed` sono label del repo); commento con job, commit, ETA, output, comando di controllo. Niente nomi utente, nodi o percorsi di home.
STOP se: `gh` rifiuta di modificare le label (permessi, UNVERIFIED per il tuo account): chiedi a LM di aggiungere `running`, ma il job resta tracciato nel commento.

### 14. Sorveglia: stato, progresso, controllo precoce della calibrazione

Stato con intervallo crescente (dalla skill; da lasciare aperto o da fare a mano ogni 15 minuti):

```bash
id=$JOB; d=60
while :; do
  out=$(bash hpc/hpc.sh status $id) || { echo "ssh failed: stop and report"; break; }
  echo "$out" | grep -Eq 'COMPLETED|FAILED|CANCELLED|TIMEOUT|OUT_OF_ME|NODE_FAIL' && { echo "$out"; break; }
  echo "$out" | grep -q 'QOSMax' && { echo "$out"; echo "blocked by QOS limit: stop and report"; break; }
  sleep $d; d=$(( d < 900 ? d * 2 : 900 ))
done
```

Progresso (riga per unita: `<set>__<dataset>: <stato>`; una unita pulita stampa esattamente `completed`):

```bash
bash hpc/hpc.sh log "$JOB" 30
bash hpc/hpc.sh run "grep -hE '^EXP-[0-9A-Za-z-]+__.*: ' logs/gb-engine_${JOB}_*.out | sed 's/^[^:]*: //' | cut -c1-40 | sort | uniq -c"
```

Atteso: nel log iniziale `code <stesso hash del passo 5> python Python 3.11.x` (la riga prima, `node <nome>`, non va incollata); il secondo comando conta, per esempio, `40 completed`
(il nome del file di log viene da `#SBATCH --job-name=gb-engine` e `--output=logs/%x_%A_%a.out`; per un job senza array `%a` diventa un numero fisso, il `*` lo copre).
Prima di eseguire il secondo comando controlla che il file esista: `bash hpc/hpc.sh run "ls logs/gb-engine_${JOB}_*.out"`; se il nome e diverso, usa quello che `hpc.sh log` ha mostrato (`== logs/...`).
Ogni volta che controlli, commenta su #52 (regola della bacheca), per esempio:
`gh issue comment 52 --body "checked: RUNNING 0h40m, 17/40 units, all completed"`.

Controllo precoce (dopo i primi minuti, le prime righe sono di C6): `completed: GATE ...` significa gate 1 fallito su quel dataset (sara escluso e riportato, regola (b); non si ripara). `error: ...` e un errore di esecuzione
(regola (c): una sola risottomissione, passo 15). `check_failed: ...` e un errore di codice (insieme "non interpretabile", nessuna risottomissione).

STOP e segnala su #52 (decidete con GB/LM se annullare con `bash hpc/hpc.sh cancel $JOB`) se: tutte le prime unita di C6 (girato su cluster) mostrano `GATE` (probabile errore sistematico di input o codice, non rumore); compare un `check_failed`;
`log` mostra `cannot activate conda env`; lo stato e `FAILED`, `TIMEOUT`, `OUT_OF_ME`, `NODE_FAIL`, `CANCELLED`. Su C4, qualche `GATE` e atteso come rischio: C4 e girato su Windows, il record lo dichiara (nota di piattaforma).
Per `TIMEOUT`, `OUT_OF_ME`, `NODE_FAIL` o righe `error` si segnala su #52 e poi si segue il passo 15 (una sola risottomissione); per `check_failed`, `cannot activate conda env` o `CANCELLED` no: si ferma tutto fino alla risposta di GB/LM.

### 15. Job finito: esito, righe, eventuale risottomissione (una sola)

```bash
bash hpc/hpc.sh status "$JOB"
bash hpc/hpc.sh run "cd ~/gb_runs/$EXP && ls rows | wc -l && grep -L '\"status\": \"completed\"' rows/*.json; echo fine"
bash hpc/hpc.sh run "grep -hE '^EXP-[0-9A-Za-z-]+__.*: ' logs/gb-engine_${JOB}_*.out | sed 's/^[^:]*: //' | cut -c1-40 | sort | uniq -c"
```

Atteso: `sacct` `COMPLETED` e `ExitCode 0:0`; `40` righe e nessun file elencato prima di `fine`; il riepilogo non contiene `error` ne `check_failed`.
(`completed: GATE ...` non rende il job `FAILED`: il codice di uscita e 2 solo per `error`/`check_failed`. Quindi un job `COMPLETED` puo contenere esclusioni per gate 1: si leggono nel report.)
Dopo una risottomissione il log del job nuovo mostra `skipped (completed)` per le unita gia finite: la verita sono i file in `rows/` (secondo comando), non l'ultimo log.

- `TIMEOUT`, `OUT_OF_ME`, `NODE_FAIL`, righe `error` o file mancanti: risottometti UNA volta lo stesso identico comando del passo 12, con lo stesso commit.
  Il codice NON ti ferma alla seconda volta (rifa ogni riga `error` a ogni risottomissione: discrepanza 13): la regola (c) del record e una disciplina tua. Una seconda risottomissione e una deviazione da segnalare su #52, non una scorciatoia.
  Prima controlla `bash hpc/hpc.sh run 'cat GIT_COMMIT'` = hash del passo 5 e NON eseguire `hpc.sh sync` (lo script rifiuta di mischiare commit nella stessa cartella e rispedirebbe l'errore; se rifiuta, lo sbatch ha gia copiato il nuovo
  `GIT_COMMIT` in `OUT/CODE_COMMIT`: il messaggio dice come ripristinarlo, e il report altrimenti mostra `MIXED CODE COMMITS`). Se `OUT_OF_ME`, e un giudizio da riportare, non da improvvisare.
  Un'unita ancora in errore dopo la risottomissione e esclusa e riportata; leggi ogni riga `error` prima di accettare l'esclusione (un bug deterministico dentro `score_models` oggi diventa un'esclusione, non ferma il set).
- `check_failed`: nessuna risottomissione. L'insieme e "non interpretabile" finche il codice non e corretto; dopo la correzione si rifanno tutte le unita su un `OUT` nuovo e vuoto, con una deviazione datata che nomina i due commit.
  Fermati e segnala a GB/LM.

STOP se: righe che non sono ne `completed` ne un errore previsto; qualunque riga che sparisce; la tentazione di cancellare una riga per far passare il report (mai).

### 16. Report (job leggero, comando del record)

```bash
JOB2=$(bash hpc/hpc.sh submit "--partition=stud --qos=stud --time=00:20:00 --cpus-per-task=1 --mem=4G -o logs/sens_report_%j.out --wrap='$A && PYTHONPATH=src python benchmarks/motif_sensitivity.py report --out \$HOME/gb_runs/$EXP'")
echo "JOB2=$JOB2"
bash hpc/hpc.sh status "$JOB2"            # ripeti fino a COMPLETED
bash hpc/hpc.sh run "(cd ~/gb_runs/$EXP && wc -c SENSITIVITY.md SENSITIVITY.json) && tail -5 logs/sens_report_${JOB2}.out | cut -c1-150"
```

(`wc -c` e non `ls -l`: niente proprietario nell'output; `logs/` e relativo alla cartella del repo sul cluster, dove `hpc.sh run` entra da solo.)

Atteso: `COMPLETED 0:0`; esistono `SENSITIVITY.md` e `SENSITIVITY.json` (due conteggi di byte, nessun `No such file`). Il report e scritto da `motif_sensitivity.py`, nessun numero e digitato a mano.
Se il log dice `source: not found`, la shell di `--wrap` e `sh` senza `source` (UNVERIFIED): risottometti (stesso report, la shell cambia soltanto; verificato solo l'espansione dei virgolettati, in locale) con

```bash
JOB2=$(bash hpc/hpc.sh submit "--partition=stud --qos=stud --time=00:20:00 --cpus-per-task=1 --mem=4G -o logs/sens_report_%j.out --wrap='bash -c \"$A && PYTHONPATH=src python benchmarks/motif_sensitivity.py report --out \$HOME/gb_runs/$EXP\"'")
```

STOP se: il job del report finisce `FAILED`/`TIMEOUT`, o fallisce con un messaggio chiaro (per esempio righe con numeri di seed diversi, o cartella senza `sens_tasks.json`): non e una riga da ritoccare, segnala.
Il fallback `bash -c` qui sopra si usa solo per `source: not found`, una volta; qualunque altro errore del report si segnala e basta.

### 17. Checksum sul cluster, fetch, copia dei soli file piccoli

```bash
bash hpc/hpc.sh run "for X in EXP-20261002-LM-1 EXP-20260925-LM-2; do cp ~/grnbench_artifacts/\$X/INPUT_SHA256SUMS ~/gb_runs/$EXP/INPUT_SHA256SUMS.\$X; done && du -sh ~/gb_runs/$EXP"
bash hpc/hpc.sh run "cd ~/gb_runs/$EXP && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS && ls"
mkdir -p "$FETCH"
bash hpc/hpc.sh fetch "~/gb_runs/$EXP" "$FETCH"
(cd "$FETCH/$EXP" && sha256sum -c --quiet SHA256SUMS && echo "fetch OK")
```

Atteso: la cartella del run e piccola (MB, non GB: contiene `rows/`, `perm_seeds/`, `sens_tasks.json`, `RUN_INPUTS.json`, `CODE_COMMIT`, `SENSITIVITY.*`, `INPUT_SHA256SUMS.*`; gli artifact grandi restano in `~/grnbench_artifacts`); `fetch OK`.
Il primo comando e quello del record (il suo `du -sh` stampa il percorso assoluto di home: non incollarlo); il secondo e l'ordinario `SHA256SUMS` degli altri esperimenti (aggiunta mia).

Copia nella cartella dei risultati (e solo tua: mai scrivere nella cartella EXP di un altro, mai riusare un ID). Se `$D` esiste gia, la prima riga si ferma e non copia nulla:

```bash
[ ! -e "$D" ] && mkdir -p "$D" \
  && cp -r "$FETCH/$EXP/rows" "$D/" \
  && cp "$FETCH/$EXP"/{SENSITIVITY.md,SENSITIVITY.json,sens_tasks.json,RUN_INPUTS.json,CODE_COMMIT} "$D/" \
  && cp "$FETCH/$EXP"/INPUT_SHA256SUMS.* "$D/" \
  && du -sh "$FETCH/$EXP/perm_seeds" "$D" \
  && echo "decidi: perm_seeds entra solo se la somma resta sotto ~5 MB (il report non rilegge perm_seeds)"
```

Solo se la somma dei due `du` sta sotto ~5 MB:

```bash
cp -r "$FETCH/$EXP/perm_seeds" "$D/"
```

Poi, in ogni caso:

```bash
(cd "$D" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS) && du -sh "$D"
grep -rIlE '(/ho[m]e/|/Users/|C:.{1,2}Users)' "$D"; echo "grep percorsi: fatto"
```

Atteso: `du -sh $D` sotto circa 5 MB (se `perm_seeds` non ci sta, non copiarla e scrivi nel record che resta sul cluster in `~/gb_runs/EXP-20261003-GB-1/perm_seeds`);
il `grep` non stampa file prima di `fatto` (aggiungi al pattern il tuo ID e il tuo utente locale, senza scriverli qui).
STOP se: la prima riga non stampa i due `du` (cartella gia esistente o una copia fallita: non riusare una cartella a meta, guarda cosa c'e e segnala); un percorso personale o un ID compare in un file
(non committare; avvisa LM, AGENTS.md sez. 11).

### 18. Rigenera il report in locale dai file copiati e confrontalo (mai a mano)

```bash
T=$(mktemp -d) && cp -r $D/. $T/ && PYTHONPATH=src python benchmarks/motif_sensitivity.py report --out $T > /dev/null \
  && diff $D/SENSITIVITY.md $T/SENSITIVITY.md && diff $D/SENSITIVITY.json $T/SENSITIVITY.json && echo "report identico"
```

Serve il venv dal lock (AGENTS.md sez. 4) sulla tua macchina. Atteso: `report identico`. STOP se: le due versioni differiscono (cifre dell'ultima posizione tra piattaforme, UNVERIFIED): non scegliere a occhio e non ritoccare nessuna delle due;
allega il `diff` (senza percorsi personali) a #52 e aspetta GB/LM. Se serve ricalcolare dopo una correzione del codice di analisi: correzione + test nella PR dei risultati, report rigenerato dalle righe invariate
(`python benchmarks/motif_sensitivity.py report --out <cartella>`), e una deviazione datata nel record (AGENTS.md sez. 7).

### 19. Leggi la calibrazione PRIMA di guardare il verdetto

```bash
grep -nE 'units planned|units without a row|INCOMPLETE|DEVIATION|MIXED CODE COMMITS|CODE_COMMIT|NOT the registered|not in experiments|MISSING|not interpretable|FAIL' $D/SENSITIVITY.md
grep -nE '^- datasets:|auto pipeline is the baseline|built-in checks|Gate 2, per set|^- (equal|DIFFERENT)' $D/SENSITIVITY.md
grep -nE '^\| (pearson_signed|pearson_motif|motif_prior) \|' $D/SENSITIVITY.md
grep -nF '| pearson_motif | pearson_signed | 20 | +0.203 | [+0.147, +0.259] | 1.62e-06 | 1 | pearson_motif |' $D/SENSITIVITY.md
grep -nF '| pearson_motif | pearson_signed | 20 | +0.182 | [+0.138, +0.227] | 1.36e-07 | 1 | pearson_motif |' $D/SENSITIVITY.md
```

Atteso (tutto da `SENSITIVITY.md`, nulla a mano):

| controllo | atteso |
|---|---|
| testata | `units planned: 40; rows: 40; permutation seeds: [100]`, `units without a row: 0`, nessun `INCOMPLETE`, nessuna riga `DEVIATION` |
| commit | `CODE_COMMIT` = hash del passo 5 senza `-dirty`; nessun `MIXED CODE COMMITS` |
| input | ogni `input references/...` e `input cistarget/...` con `(matches experiments/datasets.tsv)`; nessun `NOT the registered hash`, nessun `MISSING from RUN_INPUTS.json` |
| per set | `datasets: 20; passed gate 1 and analysed: 20; excluded (reported, never silently): 0`; datasets con auto = baseline `13/20` (C6) e `12/20` (C4), come nelle righe committed |
| gate 1 | tabella di calibrazione: `max abs diff` <= 1e-12 per modello e braccio, `edge files matching` N/N; `built-in checks ... held on N/N contexts` |
| gate 2 | `**pass**` per ciascun set e nessuna riga `DIFFERENT`; le due righe di intestazione di C6 e C4 ritrovate dal `grep -F` (il report stampa la riga committed) |

Cosa fare se qualcosa non torna (dal record, non lo decidi a posteriori):
- un dataset fallisce il gate 1: e escluso e riportato con la ragione e la dimensione dello scarto (regola (b)), mai sostituito. Per C4 e un rischio dichiarato (C4 girato su Windows). I test usano i dataset che passano.
  Elenca gli id esclusi nel commento della issue: la riesecuzione supplementare su Windows (nota 2.1 punto 5) e di LM, solo se LM e d'accordo, su una cartella `OUT` separata; lo script non ha ancora `report --label` ne un modo di combinare due `OUT`,
  quindi quell'analisi non si riporta finche non esiste una piccola PR a parte. Non farla non e una deviazione.
- gate 2 fallisce con tutti i dataset passati dal gate 1, oppure c'e un `check_failed`: l'insieme e "non interpretabile". Non si legge nessun altro verdetto di quell'insieme finche la causa non e trovata, corretta e registrata come deviazione datata.
- `MIXED CODE COMMITS` o input non registrato: non accettare, segnala.

### 20. Regola di decisione (la applica lo script; tu leggi la riga, non la ricalcoli)

Per ogni insieme, braccio `auto`, famiglia di Holm di 6 ipotesi dirette (1 coppia x 3 test x 2 direzioni, alpha 0.05), unita = dataset.
T1: `pooled(pearson_motif) > pooled(perm_db)`. T2: `macro_R(pearson_motif) > macro_R(pearson_signed)`. T3: `macro_R(pearson_motif) > macro_R(gene_tf)`. Un test e supportato se p_Holm < 0.05 nella direzione scritta.
L = pooled(real) - pooled(base); D = pooled(real) - pooled(perm_db); L_art = pooled(tfmean_db) - pooled(base); L_gene = pooled(gene_tf) - pooled(base). Precedenza (la riga piu cautelativa vince):

| regola | condizione | verdetto stampato |
|---|---|---|
| 0 | un controllo di codice e fallito, oppure gate 2 fallito con tutti i dataset passati | non interpretabile (nessun'altra regola si legge) |
| 1 | L <= 0 | inconclusivo |
| 2 | il contrario di T1 o di T2 e supportato | coerente con soli effetti di design |
| 3 | T1 e T2 supportati | "supporta evidenza a livello di target" solo se T3 supportato, D >= 0.5 L, L_art < 0.5 L, L_gene < 0.5 L e p_perm < 0.05 (robustezza di T1 al database condiviso); altrimenti misto, con le condizioni mancanti nominate (e, se manca T3 o L_gene, "consistent with a TF-generic gene propensity") |
| 4 | uno solo tra T1 e T2 supportato | misto se L_art >= 0.5 L, altrimenti inconclusivo |
| 5 | nessuno dei due | soli effetti di design se L_art >= 0.5 L, altrimenti inconclusivo |

Cosa significa per il caveat su C4/C6 (nessuno stato di claim cambia; l'esperimento e esplorativo e non crea claim di testa; i numeri registrati restano, il gate 2 li riproduce):
- supporta evidenza a livello di target: il caveat puo essere ristretto per Pearson su quell'insieme; RegDiffusion e GRNBoost2 (C5, C6) lo tengono fino allo stage 2.
- misto: il caveat resta e nomina le componenti con le loro dimensioni dal report; se T3/L_gene falliscono dice anche che il lift e coerente con una propensione genica generica.
- soli effetti di design: il caveat si rafforza; il linguaggio meccanicistico in `docs/models_that_work.md` esce; un test confermativo con `macro_R` primario richiederebbe un D-xxx.
- inconclusivo: il caveat resta com'e (inconclusivo non e equivalenza). Non interpretabile: nulla cambia finche la causa non e risolta.
- C6 (primario) e C4 (secondario) discordi: si riportano entrambi e il caveat segue il piu cautelativo, nell'ordine della nota 2.1: soli effetti di design > inconclusivo = non interpretabile > misto > supporta target-level
  (vale solo per combinare i due insiemi; dentro un insieme vale la precedenza numerata).
- il braccio `baseline` e la sua famiglia sono secondari e danno solo accordo di segno; un `Secondary checks ... DISAGREE` si nomina e non cambia il verdetto.

```bash
grep -nE 'verdict \(|registered lift L|target-identity gap|within-TF lift|design-effect control|gene-propensity control|ratios to L|T1 permutation rank|Rule 2 fired|downgraded' $D/SENSITIVITY.md
```

Da non fare: leggere un insieme "non interpretabile"; citare un numero non presente nel report; usare `perm_db` o `p_perm` come prova sul database reale (le p dei t test sono condizionate all'unico database condiviso, lo dice il report);
rieseguire sotto questo ID con modelli stage 2; modificare la sezione 1.

### 21. Cosa registrare nel record e nella PR dei risultati

Sulla branch `gb/<argomento>` nuova da `main` (non impilata su #55), file e voci:
1. `experiments/EXP-20261003-GB-1.md`, sezione 3, riempita copiando dal report: `status` (done/failed); `produced by` = `python benchmarks/motif_sensitivity.py report --out ~/gb_runs/EXP-20261003-GB-1` (job `<JOB2>`, commit del passo 5)
   e la rigenerazione locale del passo 18; `report file` = `benchmarks/results/EXP-20261003-GB-1/SENSITIVITY.md`; `headline` e `all declared comparisons` = le tabelle di T1-T3 per insieme e braccio (con n per test e dimensione della famiglia);
   `verdict` = le etichette della sezione 1 come stampate (supporta evidenza a livello di target / misto / coerente con soli effetti di design / inconclusivo / non interpretabile) per C6 e per C4, e quella combinata secondo la nota 2.1;
   `figures` = nessuna in stage 1 (solo tabelle), con la ragione data dal record, e `viz/` e `latex/` non applicabili perche la cartella non e un run di `real_bench`;
   `deviations` = nessuna oppure l'elenco datato (risottomissione per regola (c), esclusioni da gate 1 con id e ragioni, eventuale rifiuto di LM alla riesecuzione Windows, hash e commit).
   Sezione 1 non si tocca. La sezione 4 (riproduzione) resta vuota: questo run non e una riproduzione e un agente non rende mai `reproduced`.
2. `experiments/README.md`: la riga di GB-1 passa da `registered | pending` a `done` (o `failed`) con il verdetto scritto dal report.
3. `.projectos/STATE.md` (dove siamo rimasti: verdetti, job, issue) e `.projectos/MAP.md` (la nuova cartella dei risultati); riga in `CHANGELOG.md` sotto `## Unreleased`;
   `paper/CONTRIBUTIONS.md` solo se registra un lavoro nuovo (ruoli CRediT).
4. Non toccare `paper/CLAIMS.md` ne `docs/models_that_work.md` in questa PR: il testo del caveat e un passo successivo che cita la PR e le righe del report (e commenta su #41 e #45).
5. PR con `.github/pull_request_template.md`: cosa, perche, ID esperimento, test (nessun codice cambiato: scrivilo, oppure il comando eseguito con l'output vero), rischi, passaggio di consegne
   (job, commit, cartelle sul cluster `~/gb_runs/EXP-20261003-GB-1` e `~/grnbench_artifacts/` da non cancellare finche la PR non e unita). Review: LM (cartella `benchmarks/results` + codice statistico; se tocchi `paper/`, obbligatoria).

### 22. Chiudi la issue

```bash
gh issue edit 52 --remove-label running --add-label done
gh issue close 52 --comment "result: <incolla le due righe di verdetto (C6, C4) da SENSITIVITY.md>; n = <per insieme, dal report>; job <JOB>/<JOB2>, commit <hash>; results PR #<n>"
```

Se il run e fallito o non interpretabile: `--add-label failed` al posto di `done` e nel commento la causa in una riga. Non chiudere con `done` un run senza risultati letti dal report.
Commento a parte su #41/#45 con il link ai risultati, perche il testo del caveat dipende da questo verdetto.

---

## Discrepanze trovate tra record, issue, PR e script

1. **Review di LM.** #52 e #55 indicano la review di LM come passo da fare prima del lancio; su #55 esiste una `APPROVED` sotto l'account di LM sulla testa `661514b`, ma e scritta come review di agente ("not a second-human review"),
   e la descrizione della PR dice ancora "LM's approval pending". Il passo 3 chiede a LM una conferma scritta.
2. **Stima del tempo.** #52 dice "well under the plan's 1.5-2.5 h"; il record parla di 1.5-2.5 CPU-ore, cioe meno di un'ora di tempo reale con 8 worker, con limite di 4 h. Stesso senso, unita diverse.
3. **Variabili `L` e `R` del record.** I comandi del record le espandono in locale (virgolette doppie), non sul cluster: senza `export` locale diventano stringhe vuote. Ho aggiunto la guardia `: "${L:?}" "${R:?}"` (blocco delle variabili e passo 8).
4. **Il fetch cisTarget** nel record e in primo piano; qui in sottofondo (stesso script), vedi passo 9.
5. **Skill `grnbench-hpc` sez. 7** parla di `COMPLETE.json`, `real_bench.py plan/report` e `RUN_CONFIG.json`: nulla di tutto questo esiste per `motif_sensitivity.py`. Qui il run produce `RUN_INPUTS.json` (hash dei riferimenti e cisTarget, piattaforma senza nome host),
   righe con il commit, `SENSITIVITY.md/json`. L'assenza di `COMPLETE.json` non e un errore. D-020 cita `RUN_CONFIG.json`: per questo run il record rimanda a `RUN_INPUTS.json`.
6. **Vocabolario del verdetto.** Il modello di sezione 3 dice `supports | refutes | inconclusive`; la sezione 1 del record definisce cinque etichette diverse. Usa quelle della sezione 1.
7. **AGENTS.md sez. 7** chiede figure di esplorazione e `viz/`/`latex/`; il record dichiara "figures: none in stage 1 (tables only)" e la cartella dei risultati non e un run di `real_bench` (`grnbench report`/`viz` non si applicano).
   La ragione sta solo nel record; il revisore puo contestarla.
8. **`hpc.sh` e "Permission denied".** Qualunque `permission denied` in stderr, anche di un `cp` remoto su un file, fa stampare il suggerimento "Auth failed: load your key": fuorviante per i permessi dei file di LM (passo 7).
9. **Nodo di login e `/tmp`.** #41 e #24 indicano `/tmp` di un nodo preciso (il nome e in #41; qui non e scritto); il record dice che `hpc.sh run` deve raggiungere quel nodo. L'alias personale potrebbe atterrare altrove (UNVERIFIED);
   LM in #24 si e offerto di rimettere i file sul nodo di chi atterra altrove.
10. **Ordine `check`/`sync`.** Il record fa `git status`, poi `hpc.sh check`, poi `hpc.sh sync`; il passo 5 segue questo ordine (un login fallito si vede prima di spedire l'albero), e non quello della skill (che mette `check` nella sezione 1 e `sync` nella 2, senza ripeterlo).
11. **`R` nel record.** Il record dice "`R` = GB's copy of the registered references from #24"; la copia di GB non esiste ancora (in #24 LM ha messo i file in `/tmp`). Qui `R` e la cartella di LM e il passo 9 la copia nel data root dedicato.
12. **PR aperte che toccano il codice di Pearson.** #40 e #69 modificano `src/grnbench/infer/classical.py` e `base.py` (letto da `refs/prs/40` e `refs/prs/69` contro `origin/main` del 2026-10-07). Se si uniscono prima di #55, il controllo di identita del passo 4 (punto 6) non torna piu
    sulla strada A. Il passo 4 dice cosa fare (strada B sul commit `661514b`).
13. **Skill `grnbench-hpc` sez. 8 e `MAX_ATTEMPTS`.** La skill dice che una riga che fallisce di nuovo diventa `permanent_error` a `MAX_ATTEMPTS`. `motif_sensitivity.py` accetta `--max-attempts` (lo passa lo sbatch di #55) ma non lo usa: non ha stati `transient_error`/`permanent_error`,
    una riga `error` viene rifatta a OGNI risottomissione. La regola (c) del record ("una sola risottomissione") e quindi una disciplina dell'operatore, non un limite del codice (passo 15).

## UNVERIFIED (non verificato da qui; nessun comando e stato eseguito su cluster o dati reali)

- chi sei (assegnatario di #52 e GB) e se la review umana di LM e avvenuta fuori da GitHub;
- che `/tmp/grnbench_artifacts_LM/...` e `/tmp/grnbench_registered_references_LM` esistano ancora oggi, sul nodo dove atterra il tuo alias, leggibili dal tuo account, e con struttura `rows/` + `artifacts/<dataset>/{auto,baseline}/ctx*/processed.h5ad`
  (dedotta dal codice di `plan`/`run`, non vista); i conteggi `processed` e le dimensioni sono quelli del record;
- l'ambiente conda `grnbench` e la quota/spazio nella tua home; la rete del nodo di login verso il server dei file cisTarget; i tempi di download e di coda;
- che `--wrap` usi una shell con `source` (fallback indicato al passo 16);
- la durata reale del job (stima del record, mai misurata), l'esattezza dei risultati tra Windows (C4) e Linux (rischio dichiarato nel record), le cifre identiche del report tra cluster e locale;
- i permessi del tuo account GitHub per modificare le label (le label `running`, `done`, `failed` ESISTONO nel repo: verificato);
- `ssh-add -l` su Git Bash/Windows (il percorso `/c/Windows/System32/OpenSSH/ssh-add.exe` e quello che `hpc.sh` usa per ssh; per `ssh-add` e dedotto), e il comportamento di `hpc.sh` con un alias che ha un `GRNBENCH_HPC_DIR` diverso da `grnbench`;
- il nome esatto del file di log (`logs/gb-engine_<JOB>_<n>.out`, dedotto da `#SBATCH --output` e dal glob di `hpc.sh log`): il passo 14 lo fa controllare;
- tutto quello che ho eseguito io e solo locale: lettura di record, script, test e issue; simulazione dell'espansione dei virgolettati dei comandi e delle righe di controllo (conteggi, `sha256sum -c`, riepilogo del log) su cartelle finte.
