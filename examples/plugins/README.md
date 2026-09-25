# Example plugins

A grnbench plugin is one object with one method. `adata` arrives **already
normalized** (see `grnbench.infer.prepare_adata`), so `obs`, `var` and
`layers["counts"]` are reachable.

```python
class MyGRNModel:
    name = "my-model"
    def infer(self, adata, tf_list, **params) -> "pd.DataFrame":
        # -> columns [source, target, weight]  (+ optional direction ∈ {-1,0,1})
        ...
```

Return whatever column names you like — `regulator`/`tf`/`score`/`MI`/`sign` and
friends are matched case-insensitively; otherwise the first three columns are
taken positionally. Self-loops are dropped and the list is sorted by descending
weight for you.

## Files here

| file | what it is |
|---|---|
| `ridge_model.py` | a real method: per-target ridge regression on the TF panel (scikit-learn). Not a stub — usable as a third opinion next to GENIE3/GRNBoost2. |
| `ridge.json` | spec file wrapping `ridge_model.py`, so the plugin is one `--plugin` argument |
| `docker_model.json` | template for any containerized model |

## Try it

```bash
grnbench infer processed.h5ad --out runs/infer --tf-list tfs.txt \
  --plugin examples/plugins/ridge.json

# equivalently, straight at the file:
grnbench infer processed.h5ad --out runs/infer --tf-list tfs.txt \
  --plugin examples/plugins/ridge_model.py:RidgeGRN

# and in the full loop, scored against a gold standard next to the built-ins:
grnbench run --data d.h5ad --out runs/x --tf-list tfs.txt --gold dorotheaABC \
  --plugin examples/plugins/ridge.json --stability
```

`grnbench methods` then lists `ridge` alongside the built-ins, with its origin.

## Shipping it properly

Declare an entry point in your own package and grnbench discovers it with no
flags at all:

```toml
[project.entry-points."grnbench.methods"]
my-model = "my_pkg.models:MyGRNModel"
```

## Containers

grnbench writes `expr.parquet` (cells × genes), `tfs.csv` (one column, `tf`) and
`params.json` into a work directory, runs your command, then reads back
`edges.parquet`. Placeholders available in `command`:

| placeholder | value |
|---|---|
| `{workdir}` `{expr}` `{tfs}` `{out}` `{params}` | absolute host paths |
| `{expr_name}` `{tfs_name}` `{out_name}` `{params_name}` | basenames, for use under a container mount |
| `{threads}` `{seed}` | forwarded from the CLI |

Set `"format": "tsv"` if pyarrow is unavailable in your environment. A non-zero
exit code, a timeout, or a missing output file all fail loudly with the
command's stderr tail attached — a plugin never fails silently.

If the image or binary is missing, the method is *listed as unavailable with the
reason* rather than crashing the run: unavailability is data.
