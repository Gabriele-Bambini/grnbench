# grnbench: audit scientifico e percorso di pubblicazione

Data: 13 settembre 2026. Snapshot locale v0.6.0. Audit del codice e della letteratura; nessuna riparazione del codice applicativo e nessun nuovo benchmark biologico eseguito.

**Verdetto: NO-GO per il paper nella forma attuale. GO condizionato a un pilot ristretto di validazione metodologica.** Il progetto contiene infrastruttura riutilizzabile, ma la narrativa di novità è già ampiamente coperta. Alcuni risultati inferenziali sono inoltre invalidati da controesempi sul codice corrente. Non è possibile promettere pubblicabilità: il pilot deve poter concludere che il contributo non basta.

## 1. Che progetto è realmente

È un harness Python per inferenza di reti regolatorie da scRNA-seq: ingestione → QC/preprocessing → selezione HVG e TF → adapter di inferenza → valutazione di edge list → report Markdown/HTML. Non contiene un nuovo modello GNN che costituisca il contributo centrale.

| Componente | Implementazione | Valore / limite |
|---|---|---|
| Input e QC | `ingest.py`, `hygiene.py`, `qc.py`, `pipeline.py` | Formati multipli, controlli counts, aggregazione duplicati, QC esplicito. Buona infrastruttura; restano percorsi reali non coperti. |
| Preprocessing | `strategies.py`, `features.py`, `tfs.py` | Strategie separabili, HVG ∪ TF, supporto pseudobulk. La scelta del gruppo e l'identità dei TF devono rispettare il task biologico. |
| Inferenza | `infer/`, `plugins.py`, esempi | Correlazioni, MI/CLR, regressioni, RF; adapter GENIE3/GRNBoost2/ARACNe3; plugin Python/CLI. Adapter presente ≠ algoritmo validato. |
| Valutazione | `eval/`, `cli_eval.py` | AP, AUROC, precision@k, EPR, null, overlap, resampling e tier. Qui si concentrano i blocker scientifici. |
| Esecuzione | CLI, `hpc/`, manifest | Workflow pratico; ambiente non congelato e provenienza dei backtest incompleta. |

La tesi dichiarata è: «quanto della classifica è reale?». Per sostenerla servono un estimando preciso, un disegno di campionamento e controllo dell'errore. Un elenco di p-value e intervalli non fornisce automaticamente queste proprietà.

## 2. Tabella delle evidenze locali

Livelli: **R** = riprodotto con fixture eseguibile; **C** = verificato nel codice; **D** = soltanto dichiarato nella documentazione. Le righe di sorgente si riferiscono allo snapshot auditato.

| Priorità | Evidenza | Risultato | Limite / conseguenza |
|---|---|---|---|
| Bloccante, R | `eval/significance.py:69-80,134-145`; `eval_findings.json` | Il default ricalcola `k=n_pos` dopo aver rimosso metà gold. EPR osservato A=4, B=2.5; resampling dichiara B migliore con intervallo A−B [−6.8,−3.59]. Con k=100 congelato: [+0.8,+2.2]. | Cambia la quantità valutata. Non è un'incertezza sul medesimo confronto. Il fix del cutoff da solo non convalida l'inferenza. |
| Bloccante, R/C | `eval/nulls.py:26-38` | A→x,A→y,B→x,B→y diventano A→x,A→x,B→y,B→y: grado unico 2→1. Un altro input senza diagonale genera 3 self-loop su 6 archi. | Conservati i margini delle righe di un multigrafo, non il grado della rete semplice ammessa dal benchmark. Il null altera supporto e duplicazioni. |
| Bloccante, R/C | `eval/metrics.py:326-340`; `cli_eval.py:275-293` | 100 archi con gold ai rank 1 e 100: AP=.51. Troncando al primo arco: AP=1; il positivo omesso sparisce dal denominatore. | Ogni metodo viene valutato sulla propria lista; i punteggi non costituiscono un confronto sullo stesso universo. La misura condizionata alla lista può avere senso, ma non identifica superiorità fra metodi. |
| Bloccante per queste baseline, R | `infer/classical.py:79-82`; `reproduce_mi_findings.json` | Matrice tutta zero 16×3: MI(A,B)=2.079442≈ln(8), atteso 0; CLR=1. Su Bernoulli indipendenti una permutazione comune delle righe cambia i punteggi. | Il binning spezza gli ex aequo con l'indice di cella condiviso tra geni. La fixture dimostra score dipendenti dall'ordine, non un'inversione di ranking in quel caso. |
| Alta, C | `eval/significance.py:110-168` | Sottocampionamento senza reinserimento di annotazioni positive; quelle escluse diventano negative. Percentili chiamati intervallo 95%; nessuna correzione della molteplicità. | È sensibilità alla rimozione di annotazioni note, non bootstrap della verità biologica. Coverage, Type-I e validità sotto bias di annotazione non dimostrati. |
| Alta, R | `eval/significance.py:187-214` | A/B e A/C non significativi, B/C significativo: output `[A,B,C]` nello stesso tier. | Controlla soltanto il leader del tier. Inoltre non significativo ≠ equivalente. Servono matrice/grafo dei confronti e separata definizione di equivalenza. |
| Alta, C | `eval/significance.py:226-230,274-277` | Il report dice che vincitori diversi contro gold diversi non possono essere entrambi corretti. | Binding ChIP, interazione fisica e risposta a perturbazione sono oggetti diversi. Discordanza tra task non dimostra errore statistico o causalità. |
| Alta, C | `examples/plugins/ridge_model.py:53-58` | La Ridge multi-output include il target fra i regressori quando il target è un TF. | Rimuovere il self-loop dall'output non rimuove il self-predictor dal fit. Penalizza i pesi TF→TF; non descriverlo automaticamente come inflazione dei punteggi. |
| Alta, C/D | `benchmarks/beeline_backtest.py:162-170` | Seleziona i primi N indici di GeneOrdering e tutti i TF presenti. Il protocollo ufficiale corrente filtra significatività e seleziona N non-TF variabili oltre ai TF variabili. | Non dimostrata replica dello stesso protocollo. Per confronto storico fissare commit, release del corpus e selezione esatta. |
| Alta, C/D | `benchmarks/README.md:105-117`; backtest:73-108,450-452 | Il pass usa margine arbitrario 3× sopra valori pubblicati. La tabella include EPR=3.069 pur chiamandoli tutti entro ceiling=1.20. | Un massimo empirico altrui non è un bound universale o un test di leakage. «Numeri vicini a random» non dimostra correttezza. |
| Alta, C/D | `OUT/run_manifest.json`, inventario iniziale | Solo manifest con DATA.h5ad/TFS.txt, stato `no edges`, versione .5.0. Nessun risultato BEELINE CSV/JSON o edge list reale nello snapshot. | Tabella README e p=.95 documentati, non riprodotti. Gli output potrebbero esistere altrove; non si afferma che non siano mai stati prodotti. |
| Media, C + probe | `tfs.py`; `infer/correlation.py:16`; `infer/classical.py:47-51` | TF caricati uppercase, matching in inferenza exact-case: Pou5f1/Sox2 contro POU5F1/SOX2 non trova regolatori. | Il backtest usa TF nel case dell'espressione e non copre questo percorso ingest→infer. |
| Media, C + probe | `report.py:19`; `pipeline.py:318` | Fingerprint usa `read_bytes()` anche sull'input directory; fallisce per directory 10x documentate. | Il fallimento avviene nella provenienza finale, dopo la scrittura del processed file. |

Ulteriori limiti: null privi di gold hit esclusi dal denominatore dei p-value; `random_null` conserva anch'esso i conteggi marginali, diversamente dalla descrizione; stabilità combina subsampling cellule e seed algoritmico, mentre Jaccard top-k taglia arbitrariamente gli ex aequo; gold loader perde segno/tipo di evidenza. La prevalenza è una baseline convenzionale per AP, non la sua aspettativa random esatta a campione finito: con un positivo e quattro candidati E[AP_random]=.520833, prevalenza=.25.

Parti corrette da preservare: AP ai confini degli score ex aequo; AUROC con rank medi; precision@k attesa sotto tie-breaking uniforme e relativi bounds; pseudocount Monte Carlo +1. La correttezza di questi componenti non ripara il disegno sperimentale.

## 3. Letteratura: cosa esiste già

Ricerca mirata di novità, non systematic review esaustiva. Ricerca svolta il 13 settembre 2026 su pagine di autori, editori, arXiv, bioRxiv/PMC, Research Square e repository ufficiali. Le cifre seguenti sono risultati dichiarati dagli autori: non sono stati rieseguiti i loro esperimenti.

| Fonte / stato | Entità e risultato pertinente | Limite / implicazione per grnbench | Confidenza nella sovrapposizione |
|---|---|---|---|
| [BEELINE, Pratapa et al., 2020](https://www.nature.com/articles/s41592-019-0690-6), peer-reviewed | Benchmark standardizzato scRNAseq con reti sintetiche, Boolean models e reference sperimentali; accuratezza, robustezza, efficienza. | Il workflow preprocess/infer/eval non è novità. Non risolve ogni problema di inferenza statistica. | Alta |
| [geneRNIB, repository ufficiale](https://github.com/openproblems-bio/task_grn_inference), [preprint dal 2025](https://doi.org/10.1101/2025.02.25.640181) | Benchmark modulare vivo, componenti per dati/metodi/metriche e leaderboard. Repository attuale elenca 11 dataset; aggiornamento autore riporta 12 metodi e 11 metriche. | Anche estensibilità e integrazione non costituiscono da sole un contributo scientifico distinto. Versioni del preprint hanno conteggi differenti. | Alta per architettura; medio-alta per conteggi correnti |
| [Kendiukhov, arXiv:2603.03493v1](https://arxiv.org/html/2603.03493v1), 3 marzo 2026, preprint | Diagnostica inversioni dei ranking GRN lungo assi di protocollo. Tre tessuti/sei metodi; inversioni 16.3% per candidati e 32.1% per reference; CI Wilson e permutazioni dei nomi dei metodi. | È un precedente diretto della tesi attuale. Non usa bootstrap né rewiring degree-preserving: non confondere i null. | Alta |
| [Three Classes of Confound…, Research Square](https://doi.org/10.21203/rs.3.rs-8997629/v1), 26 marzo 2026, preprint | Toolkit `grn_confound_audit`: batch/donor, prossimità, topologia; double-edge swaps diretti e report. | Anche audit topologico GRN ha precedente. Panel piccolo e 48 null; la validità dei suoi test non è verificata qui. Problemi sospetti non vanno presentati come errori dimostrati. | Alta sul precedente |
| [Partial-label metric ceilings…, BioSystems](https://doi.org/10.1016/j.biosystems.2026.105864), online 30 giugno 2026 | Ceiling e normalizzazioni sotto annotazioni parziali, coverage incerta e sensibilità a missingness non uniforme. | «Gold incompleto distorce le metriche» non è nuova tesi. Le assunzioni di osservabilità restano cruciali. | Alta |
| [challengeR, Wiesenfarth et al., 2021](https://www.nature.com/articles/s41598-021-82017-6), peer-reviewed | Bootstrap dei ranking, confronti pairwise, significance ranking e visualizzazione dell'incertezza dei benchmark. | Bootstrap e rinuncia al ranking totale sono già noti. Il dominio e l'unità campionaria differiscono dal GRN. | Alta |
| [Inference for Ranks, Mogstad et al.](https://doi.org/10.1093/restud/rdad006), online 2023, volume 2024; [csranks](https://search.r-project.org/CRAN/refmans/csranks/html/csranks.html) | Confidence sets simultanei per ranghi e confronti multipli; implementazioni bootstrap/stepdown. | Anche ordinamenti parziali con inferenza simultanea hanno teoria precedente. Applicarla non basta a rivendicare un nuovo metodo. | Alta |
| [FERRET, Eicher e Quackenbush](https://pmc.ncbi.nlm.nih.gov/articles/PMC13419529/), preprint 25 luglio 2026, DOI 10.64898/2026.07.21.739795 | Valuta similarità entro/fra tipi cellulari e comportamento attraverso cutoffs senza un gold completo. | «Stabilità biologica senza gold» non è un'area libera. Stabilità resta distinta dalla correttezza causale. | Alta sul tema; dettagli non auditati |
| [BEAR-GRN, repository ufficiale](https://github.com/UzunLab/BEAR-GRN) | Preprocessing, inferenza, valutazione, subsampling/Jaccard ed estensibilità. | Precedente software. Data/pubblicazione scientifica non stabilite; non usare il placeholder bibliografico del README come citazione valida. | Alta esistenza; bassa bibliografia |
| [CausalBench, 2025](https://www.nature.com/articles/s42003-025-07764-y), peer-reviewed | Inferenza di reti su grandi dati single-cell interventional, valutazione biologica e quantitativa. | Aggiungere genericamente Perturb-seq non crea novità. Offre un possibile substrato di validazione con task esplicito. | Alta |

Ulteriore segnalazione da verificare prima di fissare il contributo: il manoscritto autore [From ceilings to corrected estimates…](https://www.researchgate.net/publication/412127718_From_ceilings_to_corrected_estimates_capture-recapture_on_reference_networks_rescales_but_rarely_reorders_gene_regulatory_network_benchmarks), agosto 2026, tratta capture–recapture, correzioni di osservabilità e bootstrap per TF. Identità/versione archivistica non risolte qui: non usato come pilastro del verdetto. Anche l'[audit recente di perturbazioni e reti curate](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7336418) è un segnale di concorrenza, non evidenza sperimentale convalidata in questo audit.

La fonte dei valori 1.10/1.20 è [DAZZLE, PLOS Computational Biology 2025](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013603): sono massimi osservati in un confronto, non limiti matematici. Il [codice BEELINE](https://github.com/Murali-group/Beeline/blob/master/utils/generateExpInputs.py) e il [corpus Zenodo v4](https://zenodo.org/records/19009603) devono essere congelati insieme; il fetch locale usa v4, compatibile con BEELINE aggiornato nel marzo 2026. Una replica storica richiede una scelta diversa, esplicita e verificata.

## 4. Tre direzioni possibili

| Direzione | Giudizio | Cosa occorre |
|---|---|---|
| Pubblicare l'attuale «benchmark con plugin, null e ranking instabili» | **NO-GO** | Sovrapposizione diretta e fondamenta statistiche difettose. Più adapter o un dashboard non risolvono. |
| Rendere grnbench un software affidabile, possibilmente componente di geneRNIB | **GO ingegneristico** | Riparazioni, test avversariali, ambiente congelato, artefatti e adozione esterna. Un eventuale software paper dipende da utilità dimostrata, non dalla sola quantità di codice. |
| Studio di metascienza: quando un benchmark GRN supporta davvero una dichiarazione di superiorità | **GO condizionato; direzione raccomandata** | Problema esterno riprodotto, estimando fissato, calibrazione empirica e confronto con statistica già valida. Nuovo algoritmo opzionale: non va inventato se bastano metodi esistenti. |

Titolo di lavoro, non claim acquisito: **When do GRN benchmarks support superiority claims? Calibration under dependent and incomplete reference evidence.**

## 5. Disegno scientifico proposto

### Domanda e oggetto misurato

Domanda primaria: «Le procedure che dichiarano un metodo GRN superiore controllano la probabilità di una falsa dichiarazione quando dati, reference e archi hanno la struttura di dipendenza prevista dal task?»

Congelare task T, universo U di coppie TF×target ammissibili e budget k. Per esempio precision@100, se biologicamente giustificato dal budget di validazione; AP come endpoint secondario. Il valore 100 è una proposta di protocollo, non un optimum da scegliere dopo aver visto i risultati.

Definire Δ_ab = E[M(f_a(D),R;U,k) − M(f_b(D),R;U,k)], indicando esattamente quali acquisizioni di D e R potrebbero essere ripetute e quale popolazione rappresentino. Con un dataset e una reference fissati, il punteggio è deterministico. Non si ottiene incertezza biologica inventando un campionamento degli archi.

- Endpoint primario possibile: recupero di interazioni annotate da una specifica modalità sperimentale nello stesso contesto. È **silver-reference recovery**, non ricostruzione garantita della GRN vera.
- Binding ChIP, membership in PPI, risposta downstream a perturbazione e regolazione diretta sono task diversi. Non aggregarli in un unico gold per poi interpretare le differenze come contraddizioni.
- Sul sintetico, distinguere superiorità rispetto al grafo latente e rispetto alle annotazioni osservate. Missingness biased può cambiare il secondo estimando legittimamente.
- Split, unità indipendenti e iperparametri fissati prima della valutazione; seed algoritmici e subsampling delle cellule costituiscono fattori separati.

### Contratto di valutazione necessario

1. Identificatori di specie/versione espliciti; mapping verificato senza collisioni silenziose. U congelato prima delle predizioni; hash completo di U e gold.
2. Ogni metodo deve produrre uno score per U, oppure dichiarare score assenti. Policy prefissata per astensione/output sparso, ad esempio blocco ex aequo in fondo per candidati non proposti; riportare la sensitivity a tale convenzione. Vietato ridurre U separatamente per ogni metodo.
3. Duplicati, diagonale, score NaN e segno gestiti da un unico contratto. Top-k fissato anche durante i resample. AP non interpolata nominata come tale.
4. Le annotazioni non osservate non sono negativi biologici certificati. Gold thinning rinominato sensitivity analysis; nessuna promessa di correggere MNAR senza assunzioni identificanti.
5. Null distinti: permutazione di score su U per un'ipotesi di ranking casuale; eventuale null di hub con condizionamento esplicito. Edge-swap solo su reti semplici sparse con vincoli, mixing e supporto verificati. Il supporto bipartito completo non offre un rewiring topologico informativo.
6. Decisioni espresse come confronti supportati/non risolti. Equivalenza soltanto con margine pratico prefissato. Non trasformare non-significatività in tier di equivalenza.

### Esperimenti e comparator

**E0 — validità del motore.** Trasformare le fixture allegate in regressioni dopo la riparazione; aggiungere invarianza all'ordine, candidati mancanti, diagonale/duplicati, TF mouse e input 10x directory. Confrontare numericamente lo stesso output di un metodo con evaluator ufficiale e implementazione indipendente della stessa metrica, usando la stessa convenzione per i tie.

**E1 — pilot su un protocollo esterno.** Riprodurre fedelmente almeno una procedura pubblica che formula una concreta dichiarazione di superiorità. Documentare assunzioni, dati, versioni e claim. Non attribuire a BEELINE o geneRNIB i bug trovati qui. Se nessun claim esterno risulta problematico, il paper deve registrare questo risultato o cambiare obiettivo; non costruire un avversario artificiale.

**E2 — calibrazione controllata.** Prima simulazioni statistiche economiche con null globale e parziali, effetti quasi nulli, dipendenza entro/fra TF, target hub condivisi e annotazioni uniformemente o selettivamente mancanti. Variare 4/8/16/32 unità indipendenti per mettere in crisi le ipotesi. Poi due famiglie di simulatori meccanicistici: [SERGIO](https://doi.org/10.1016/j.cels.2020.08.003) e [BoolODE](https://github.com/murali-group/BoolODE), con topologie e scenari tenuti fuori dallo sviluppo. Sono famiglie diverse, non prove biologiche indipendenti.

Conoscere il grafo simulato non significa conoscere E[M_a]−E[M_b]. Stimare la superiorità vera con un Monte Carlo esterno indipendente e sufficientemente preciso; usare repliche interne distinte per misurare coverage e false superiority. Pilot di 200 repliche per scenario; conferma almeno 1.000 per gli scenari principali, aumentate in base all'errore Monte Carlo richiesto. Registrare intervalli binomiali per i tassi osservati.

Comparator minimi: ranking descrittivo; gold thinning attuale come diagnostica, non baseline statistica valida; test paired sulla stessa statistica con Holm; simultaneità max-T/stepdown e confidence sets per ranghi ove applicabili; challengeR come pratica di riferimento. Wilcoxon di challengeR non testa in generale la differenza media AP: non confrontare estimandi diversi. Confrontare potenza a errore empirico abbinato e costo simile. Holm è valido sotto dipendenza se i p-value individuali sono validi: non è il nemico da battere per definizione.

Metriche primarie di audit: probabilità di almeno una falsa superiorità nella famiglia dichiarata, coverage simultanea, potenza per un margine pratico prefissato, quota di decisioni non risolte e costo. Una procedura che si astiene sempre non è una soluzione utile. Una procedura calibrata sul proprio modello può fallire quando TF condividono target o i cluster non sono indipendenti: dimostrarne i limiti.

**E3 — verifica reale.** Primo perimetro: 6–8 metodi computazionalmente accessibili, includendo Spearman, Pearson, Ridge corretta, GENIE3, GRNBoost2 e metodi con struttura differente; MI/CLR solo dopo correzione. Track con prior cromatinici o modelli preaddestrati separata per informazione disponibile e rischio di contaminazione. Usare BEELINE/geneRNIB per riuso e confronto, senza moltiplicare piattaforme.

Per la validazione sperimentale scegliere contesti con guide/reagenti indipendenti e, dove disponibili, repliche biologiche o donatori. Separare scoperta e conferma prima di sviluppare i criteri. Due guide contro lo stesso gene non sono due donatori, e lo split casuale delle cellule non è una replica sperimentale. Efficienza di perturbazione, batch e potenza del test di replica devono essere misurati. Una decisione non replicata può essere underpowered o riguardare un task diverso; non chiamarla automaticamente falso positivo. Le risposte a CRISPR possono essere indirette: nessuna promessa di identificare archi causali diretti.

**E4 — utilità reale del risultato.** Il contributo deve mostrare quali decisioni cambiano rispetto ai comparator validi e perché conta per un budget di validazione. La sola percentuale di tier o di ranking invertiti è insufficiente. Se una baseline standard risolve il problema, il risultato onesto è un audit con raccomandazioni riproducibili, non una nuova teoria.

### Gate per continuare o fermarsi

| Gate | Passa se | Se fallisce |
|---|---|---|
| G0, software | Fixture riparate; test passano; stesso output/stesso evaluator producono risultati concordanti. | Fermare i benchmark scientifici. |
| G1, novità | Problema documentato in almeno un protocollo esterno, distinto dai soli bug locali e dai precedenti sopra. | Mantenere software/portfolio; nessun claim di nuovo metodo. |
| G2, identificabilità | Task e unità indipendenti consentono di stimare l'errore rivendicato. | Dichiarare non stimabile, ridurre claim o cambiare dati. |
| G3, calibrazione | Tolleranza ε e precisione Monte Carlo prefissate: limite superiore unilaterale dell'errore ≤ α+ε e limite inferiore della coverage ≥ 1−α−ε negli scenari confermativi, con potenza utile. | Intervalli troppo ampi = inconclusivo, aumentare repliche; errore eccessivo = correggere assunzioni/metodo o riportare il limite. |
| G4, valore aggiunto | Beneficio rispetto a Holm/max-T/csranks correttamente applicati, oppure importante evidenza empirica nuova su protocolli reali. | Non vendere un wrapper statistico come innovazione. |
| G5, conferma | Decisioni e limiti verificati su scenari/dati separati dallo sviluppo; nessuna selezione post hoc dei soli casi favorevoli. | Risultato esplorativo, non confermativo. |

Il pilot dovrebbe precedere qualsiasi grande campagna HPC. Ordine di lavoro: correzioni e contratto → replica esterna → simulazioni economiche → decisione G1/G2 → espansione. Tempi e RAM del confronto completo vanno stimati dal pilot; non sono stati misurati qui. Aggiungere subito una nuova GNN aumenterebbe il costo senza risolvere il problema identificato.

## 6. Artefatti necessari prima di chiamarlo pubblicabile

- Repository con storia/versione identificabile; questo snapshot locale non contiene `.git`.
- Ambiente congelato, CI e risultato dei test; niente fallback silenzioso da installazione completa a minima presentato come successo.
- Manifest completi: hash di input anche multi-GB, inventario ricorsivo delle directory, feature/TF/gold effettivi, seed, parametri, commit, versioni di tutte le dipendenze e container digest.
- Edge list persistite per ogni metodo e seed. Il backtest corrente scarta `_prov` e conserva gli archi solo in memoria.
- Dataset/release/licenze, mapping, split e policy di accesso tracciati. Verificare le licenze prima della redistribuzione.
- Tabella/figure generate esclusivamente dagli artefatti, risultati negativi inclusi; separazione esplicita tra dati, inferenze e limiti.
- Protocollo preregistrato o congelato prima degli esperimenti confermativi; revisione statistica indipendente.

Una versione software solida sarebbe già un buon progetto. Un paper convincente richiede un risultato verificabile che cambi ciò che si può concludere dai benchmark. Il lavoro necessario è soprattutto nel disegno sperimentale e nell'evidenza, non nell'aggiunta di funzionalità.

## 7. Riproduzione dell'audit e limiti

Comandi dal root del progetto:

```powershell
rtk proxy python audit/2026-09-13/reproduce_eval_findings.py
rtk proxy python audit/2026-09-13/reproduce_mi_findings.py
```

Eseguiti anche dal coordinatore: sei finding di valutazione riprodotti; MI costante e violazione dell'invarianza all'ordine riprodotte. Questi script caricano i moduli reali con importlib, evitando soltanto gli initializer che richiedono scanpy. Sono fixture dell'audit, non certificazione dell'intero pacchetto né test della riparazione.

L'esecuzione completa di pytest e l'ambiente osservato sono registrati in `pytest.log`, `pytest.xml`, `environment-observed.txt` e `VERIFICA.md`. L'ambiente di audit è esterno al progetto e usa system-site-packages; non equivale a una build ermetica. I dati BEELINE e i risultati numerici del README non sono stati rieseguiti in questo audit.

File applicativi invariati. Creati solo dossier/fixture/evidenze e collegamento Graphify richiesto dalle regole globali. Il grafo condiviso iniziale non conteneva il progetto; il rapporto del parent era precedente a questo snapshot. L'evidenza scientifica locale deriva dal codice e dalle esecuzioni, non dal grafo.

## 8. Post-mortem compatto

Problema ricorrente: test che controllano coerenza interna e documentazione persuasiva hanno fatto apparire validato un estimando diverso da quello dichiarato. Azione proposta: generare tabelle README da output versionati e introdurre proprietà scientifiche avversariali — universo fisso, cutoff fisso, costanti con MI zero, invarianza all'ordine, null con supporto valido — accanto ai test del percorso ordinario. Nessuna memoria globale aggiornata.
