# ProjectOS — memoria di progetto cross-agent

Un solo sistema per Claude Code, Codex, Antigravity (e qualsiasi altro agente o umano):
ogni progetto ha una cartella `.projectos/` in markdown con **scopo, mappa, stato, tentativi, decisioni, checkpoint**.
Nessun agente rilegge la codebase, nessuno ripete un tentativo fallito, nessun lavoro perso quando finisce l'usage.

## Cosa c'è nel pacchetto

| file | ruolo |
|---|---|
| `engine/pos.py` | CLI `pos` + hook (Python stdlib, un file). Scan, brief, checkpoint, log, usage, sessioni |
| `engine/PROTOCOL.md` | il protocollo unico (copiato in ogni progetto) |
| `skills/` | 5 skill: `projectos`, `pos-map`, `pos-checkpoint`, `pos-log`, `pos-resume` |
| `agents/claude/` | subagent `pos-mapper` (haiku) e `pos-architect` (sonnet) |
| `install.py` | installa tutto su questa macchina (idempotente, backup automatici) |

## Per progetto: `.projectos/`

| file | contenuto |
|---|---|
| `PURPOSE.md` | perché esiste, quando è "finito", criteri di successo misurabili, vincoli, non-obiettivi (tuo) |
| `MAP.md` | architettura: come si lancia, moduli, flusso dati, gotcha (scritta da modello economico) |
| `STATE.md` | focus, fatto, in corso, prossimi passi, blocchi, note per il prossimo agente |
| `ATTEMPTS.md` | log append-only dei tentativi con esito, metriche, evidenze, **lezione** |
| `DECISIONS.md` | decisioni (ADR-lite), con supersede |
| `checkpoints/` | snapshot immutabili (manuali + automatici) |
| `scan.md` | scan meccanico gratuito: layout, dipendenze, entry point, indice simboli, classe di grandezza |

Più: blocco breve in `AGENTS.md` (Codex), `CLAUDE.md` (Claude), `.agent/rules/projectos.md` (Antigravity) che punta al protocollo. Le istruzioni vivono in un solo posto → niente overlap tra agenti.

## Come si attiva (automatico)

1. Apri Claude/Codex in una cartella → hook **SessionStart**:
   - progetto già inizializzato → inietta un *brief* (~1.5k token): scopo, stato, ultimo checkpoint, tentativi falliti con lezioni, decisioni in vigore, freschezza mappa, altri agenti attivi, usage;
   - non inizializzato → `pos init` (scan tier-0, **zero token LLM**, + bootstrap graphify hub) e l'agente viene istruito a mappare con **haiku** in background;
   - cartella "contenitore" (≥2 sottoprogetti, es. `progetti GNN`) → non inizializza, suggerisce `pos discover`.
   - Ignorate: home, Desktop/Downloads/Documents, root dei drive, AppData, Temp, `ignore_paths` in `~/.projectos/config.json`, cartelle con `.projectos-ignore`.
2. Antigravity (niente hook): la regola globale in `~/.gemini/GEMINI.md` gli dice di fare `pos brief` / `pos init` come primo passo.

## Checkpoint e usage

| trigger | Claude Code | Codex | Antigravity |
|---|---|---|---|
| fine unità di lavoro | agente (protocollo + skill) | agente | agente |
| usage ≥75% | promemoria (UserPromptSubmit) | promemoria | — |
| usage ≥88% e ≥95% | **checkpoint forzato** (Stop hook blocca e ordina) | forzato (Stop hook) | — |
| 60 min con file modificati e nessun checkpoint | nudge (Stop hook) | nudge | — |
| prima della compattazione del contesto | auto-checkpoint meccanico (file modificati + richieste recenti) | idem | — |
| fine sessione | auto-checkpoint meccanico se ci sono modifiche | idem | — |

Fonti di usage: Claude CLI → statusline (`rate_limits` 5h/7d + contesto esatto); Claude desktop → solo token dal transcript
(percentuale solo se la finestra è nota: `context_windows` in config); Codex → rollout (`rate_limits` + contesto esatto).
Soglie in `~/.projectos/config.json` (`warn_pct`, `force_pct`, `periodic_checkpoint_min`).

## Economia modelli

| tier | chi | per cosa |
|---|---|---|
| T0 | script `pos` | scan, fingerprint, brief, checkpoint automatici — 0 token |
| T1 | haiku / mini | mappatura, riassunti (`pos-mapper`) |
| T2 | sonnet | sintesi mappa solo per progetti grandi (classe L, `pos-architect`) |
| T3 | opus / top | solo ragionamento nuovo; **mai** mappatura |

Piano di mappatura automatico dalla classe di grandezza: S → 1 haiku; M → 2-3 haiku in parallelo + merge haiku;
L → fino a 5 haiku + 1 sonnet che legge solo le parti. Rimappatura **incrementale**: solo i gruppi di file cambiati (fingerprint).

## Comandi utili

```
pos brief                     # contesto compatto
pos status                    # mappa fresca/stale, checkpoint, sessioni, usage
pos attempt list --grep gat   # cosa è già stato provato
pos attempt add -t "..." -s failed --result "..." --lesson "..."
pos decide -t "..." --decision "..." --why "..."
pos checkpoint -m "..."       # snapshot
pos usage                     # contesto / limiti
pos claim "task"              # evita overlap con altri agenti attivi
pos export --out HANDOFF.md   # passaggio a un agente senza accesso ai file
pos list                      # tutti i progetti (entrambi i PC, via hub)
pos discover "<drive>" --init   # inizializza tutto (tier-0), poi skill pos-map mode all
pos doctor                    # verifica installazione
```

## Installazione / altro PC / disinstallazione

- Installa o aggiorna: `python "<hub>\shared\projectos\install.py"` (anche STEP 8 di `install_master.ps1`).
- Dopo aver modificato il pacchetto: rilancia `install.py`, poi `pos upgrade --all` per aggiornare `PROTOCOL.md` nei progetti.
- Codex: al primo avvio approva gli hook con `/hooks` (richiesto da Codex per sicurezza).
- Disinstalla: `python install.py --uninstall` (tiene `.projectos/` dei progetti e `~/.projectos`).
- Backup di ogni config toccata: `~/.projectos/backups/<timestamp>/`.

## Limiti noti

- Claude desktop app non esegue la statusline → niente % di rate limit 5h/7d lì; restano: checkpoint per unità di lavoro, nudge periodico, auto-checkpoint pre-compattazione e fine sessione.
- Antigravity non ha hook → dipende dalla disciplina dell'agente (istruito a checkpointare dopo ogni unità).
- Subagent nuovi (`pos-mapper`) visibili a Claude solo dalla sessione successiva all'installazione.
