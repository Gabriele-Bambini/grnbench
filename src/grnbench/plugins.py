"""Custom-model plugin contract — bring your own GRN method.

Two ways in, both landing in the same registry as the built-in adapters and
therefore in the same edge schema, the same manifests and the same evaluation:

1. **Python object** (in-process). Anything exposing ``infer``:

   .. code-block:: python

       class MyGRNModel:
           name = "my-model"
           description = "..."                       # optional
           def infer(self, adata, tf_list, **params) -> pd.DataFrame:
               # -> columns [source, target, weight] (+ optional direction)
               ...

   ``adata`` is the *normalized* AnnData (see
   :func:`grnbench.infer.base.prepare_adata`), so ``obs``/``var``/
   ``layers["counts"]`` are all reachable. A bare function
   ``f(adata, tf_list, **params)`` works too. Register it either by declaring
   an entry point in your own package::

       [project.entry-points."grnbench.methods"]
       my-model = "my_pkg.models:MyGRNModel"

   (discovered automatically, no flags needed) or ad hoc with
   ``--plugin my_pkg.models:MyGRNModel``.

2. **External command / Docker image** (out-of-process), described by a JSON or
   YAML file::

       {
         "name": "my-docker-model",
         "description": "...",
         "command": ["docker", "run", "--rm",
                     "-v", "{workdir}:/data", "myimage:1.0",
                     "--expr", "/data/{expr_name}",
                     "--tfs", "/data/{tfs_name}",
                     "--out", "/data/{out_name}"],
         "format": "parquet",
         "timeout_s": 3600
       }

   grnbench writes ``expr.parquet`` (cells x genes), ``tfs.csv`` and
   ``params.json`` into a work directory, runs the command, then reads back
   ``edges.parquet``. Placeholders: ``{workdir} {expr} {tfs} {out} {params}``
   for absolute host paths, ``{expr_name} {tfs_name} {out_name}
   {params_name}`` for the basenames a container sees under its own mount, plus
   ``{threads}`` and ``{seed}``.

Availability is data, never an exception at import time: a plugin whose image,
binary or import is missing stays listed with its reason and is skipped loudly.
"""
from __future__ import annotations

import json
import os
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from grnbench.infer.base import MethodError, MethodSpec, register_method

ENTRY_POINT_GROUP = "grnbench.methods"

#: canonical output columns a plugin must (eventually) produce
_EDGE_ALIASES = {
    "source": {"source", "regulator", "tf", "tf1", "gene1", "from",
               "transcription_factor"},
    "target": {"target", "target_gene", "gene2", "to"},
    "weight": {"weight", "score", "importance", "mi", "value", "coefficient"},
    "direction": {"direction", "sign", "mode", "effect"},
}

_loaded_entry_points = False


class PluginError(MethodError):
    """A plugin could not be loaded, or violated the adapter contract."""


# --------------------------------------------------------------------------- #
# discovery
# --------------------------------------------------------------------------- #

def ensure_plugins_loaded() -> None:
    """Load ``grnbench.methods`` entry points once per process (idempotent)."""
    global _loaded_entry_points
    if _loaded_entry_points:
        return
    _loaded_entry_points = True  # set first: a failure must not retry-loop
    for spec in load_entry_point_methods():
        try:
            register_method(spec)
        except MethodError:
            # a third-party package must never be able to shadow a built-in
            continue


def load_entry_point_methods() -> list[MethodSpec]:
    """Discover installed packages advertising ``grnbench.methods`` entry points.

    A plugin that fails to import is returned as an *unavailable* spec carrying
    the import error, so ``grnbench-infer --methods list`` shows the breakage
    instead of the method silently vanishing.
    """
    from importlib.metadata import entry_points

    specs: list[MethodSpec] = []
    try:
        eps = list(entry_points(group=ENTRY_POINT_GROUP))
    except Exception:  # noqa: BLE001 - metadata backends vary across environments
        return specs
    for ep in eps:
        try:
            obj = ep.load()
        except Exception as e:  # noqa: BLE001 - a broken plugin is data
            specs.append(MethodSpec(
                name=ep.name,
                description=f"plugin from entry point '{ep.value}'",
                runner=None,
                unavailable=f"import failed: {type(e).__name__}: {e}",
                origin="entry_point",
            ))
            continue
        try:
            specs.append(python_method_spec(obj, name=ep.name, origin="entry_point"))
        except PluginError as e:
            specs.append(MethodSpec(
                name=ep.name,
                description=f"plugin from entry point '{ep.value}'",
                runner=None,
                unavailable=str(e),
                origin="entry_point",
            ))
    return specs


# --------------------------------------------------------------------------- #
# python plugins
# --------------------------------------------------------------------------- #

def python_method_spec(obj,
                       name: str | None = None,
                       description: str | None = None,
                       origin: str = "python") -> MethodSpec:
    """Wrap a class / instance / function into a :class:`MethodSpec`.

    Accepted shapes, in priority order: a class with ``infer`` (instantiated
    with no arguments), an object with ``infer``, or a plain callable. The
    resulting runner is called as ``runner(adata, tf_list, **params)``.
    """
    if isinstance(obj, type):
        if not hasattr(obj, "infer"):
            raise PluginError(
                f"plugin class {obj.__qualname__} has no .infer(adata, tf_list, **params)"
            )
        try:
            obj = obj()
        except Exception as e:  # noqa: BLE001
            raise PluginError(
                f"plugin class {obj.__qualname__} could not be instantiated with no "
                f"arguments: {type(e).__name__}: {e}"
            ) from e

    if hasattr(obj, "infer"):
        runner = obj.infer
    elif callable(obj):
        runner = obj
    else:
        raise PluginError(
            f"plugin object {obj!r} is neither callable nor exposes .infer()"
        )

    resolved = name or getattr(obj, "name", None) or getattr(obj, "__name__", None)
    if not resolved:
        raise PluginError("plugin has no name: set a `name` attribute or pass one")
    desc = (description or getattr(obj, "description", None)
            or (getattr(obj, "__doc__", None) or "").strip().splitlines()[:1])
    if isinstance(desc, list):
        desc = desc[0] if desc else ""
    return MethodSpec(
        name=str(resolved),
        description=str(desc) or f"custom {origin} plugin",
        runner=runner,
        takes="adata",
        origin=origin,
    )


def split_object_ref(spec: str) -> tuple[str | None, str | None]:
    """Split ``"head:Attribute"`` on the LAST colon, tolerating Windows drives.

    ``C:\\models\\x.py`` (no attribute) returns ``(None, None)`` because the
    trailing part is a path, not an identifier.
    """
    head, sep, attr = spec.rpartition(":")
    if not sep or not head or not attr or "/" in attr or "\\" in attr:
        return None, None
    return head, attr


def _object_from_file(py_path: str, attr: str):
    """Import ``attr`` out of a .py file that is not on ``sys.path``."""
    import importlib.util

    p = Path(py_path)
    if not p.exists():
        raise PluginError(f"plugin file not found: {p}")
    mod_spec = importlib.util.spec_from_file_location(f"grnbench_plugin_{p.stem}", p)
    if mod_spec is None or mod_spec.loader is None:  # pragma: no cover
        raise PluginError(f"cannot load a Python module from {p}")
    module = importlib.util.module_from_spec(mod_spec)
    try:
        mod_spec.loader.exec_module(module)
    except Exception as e:  # noqa: BLE001
        raise PluginError(f"error executing {p}: {type(e).__name__}: {e}") from e
    try:
        return getattr(module, attr)
    except AttributeError as e:
        raise PluginError(f"{p.name} has no attribute '{attr}'") from e


def import_object(spec: str):
    """Import ``"package.module:Attribute"`` and return the attribute."""
    if ":" not in spec:
        raise PluginError(
            f"'{spec}' is not a module path: expected 'package.module:Attribute'"
        )
    mod_name, _, attr = spec.rpartition(":")
    import importlib
    try:
        mod = importlib.import_module(mod_name)
    except Exception as e:  # noqa: BLE001
        raise PluginError(f"cannot import '{mod_name}': {type(e).__name__}: {e}") from e
    try:
        return getattr(mod, attr)
    except AttributeError as e:
        raise PluginError(f"module '{mod_name}' has no attribute '{attr}'") from e


def resolve_object(ref: str, relative_to: Path | None = None):
    """Resolve a Python object reference in either supported spelling.

    ``package.module:Attribute`` (importable) or ``path/to/file.py:Attribute``
    (loaded straight off disk, so a one-off model needs no packaging). Relative
    file paths inside a spec file resolve against that file's directory.
    """
    head, attr = split_object_ref(ref)
    if head is None:
        raise PluginError(
            f"'{ref}' is not an object reference: expected "
            "'package.module:Attribute' or 'path/to/file.py:Attribute'"
        )
    if head.lower().endswith(".py"):
        p = Path(head)
        if not p.is_absolute() and relative_to is not None:
            candidate = relative_to / p
            if candidate.exists():
                p = candidate
        return _object_from_file(str(p), attr)
    return import_object(f"{head}:{attr}")


# --------------------------------------------------------------------------- #
# external command / docker plugins
# --------------------------------------------------------------------------- #

def _write_expression(expr: pd.DataFrame, path: Path, fmt: str) -> None:
    if fmt == "parquet":
        try:
            expr.to_parquet(path)
        except ImportError as e:  # pragma: no cover - guarded by the probe
            raise PluginError(
                "parquet exchange needs pyarrow: pip install 'grnbench[plugins]', "
                "or set \"format\": \"tsv\" in the plugin spec"
            ) from e
    else:
        expr.to_csv(path, sep="\t")


def _read_edges(path: Path, fmt: str) -> pd.DataFrame:
    if not path.exists():
        raise PluginError(
            f"plugin command finished successfully but wrote no {path.name}; "
            "the command must write its edge table to the --out path"
        )
    df = pd.read_parquet(path) if fmt == "parquet" else pd.read_csv(path, sep="\t")
    return canonical_edge_frame(df)


def canonical_edge_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Project an arbitrary edge table onto [source, target, weight(, direction)].

    Column names are matched case-insensitively against known aliases; when no
    alias matches, the first three positional columns are used (and the failure
    mode is loud if there are fewer than three).
    """
    if df.shape[1] < 3:
        raise PluginError(
            f"plugin edge table needs >= 3 columns [source, target, weight]; "
            f"got {list(df.columns)}"
        )
    found: dict[str, str] = {}
    # exact canonical names win over aliases: a table carrying BOTH "source" and
    # "regulator" must not resolve by column order.
    for canon in _EDGE_ALIASES:
        for col in df.columns:
            if str(col).strip().lower() == canon:
                found[canon] = col
                break
    for col in df.columns:
        low = str(col).strip().lower()
        for canon, aliases in _EDGE_ALIASES.items():
            if low in aliases and canon not in found:
                found[canon] = col
    if "source" not in found or "target" not in found or "weight" not in found:
        cols = list(df.columns)
        found = {"source": cols[0], "target": cols[1], "weight": cols[2]}
    out = pd.DataFrame({
        "source": df[found["source"]].astype(str),
        "target": df[found["target"]].astype(str),
        "weight": pd.to_numeric(df[found["weight"]], errors="coerce").astype(float),
    })
    if "direction" in found:
        out["direction"] = pd.to_numeric(
            df[found["direction"]], errors="coerce").fillna(0).astype(int)
    bad = int(out["weight"].isna().sum())
    if bad:
        raise PluginError(
            f"plugin edge table has {bad} non-numeric weight(s) in column "
            f"'{found['weight']}'"
        )
    return out


def _command_probe(cmd: list[str], fmt: str, name: str):
    def probe() -> tuple[bool, str]:
        exe = cmd[0]
        if os.sep in exe or (os.altsep and os.altsep in exe):
            if not Path(exe).exists():
                return False, f"plugin '{name}': executable '{exe}' not found"
        elif shutil.which(exe) is None:
            return False, f"plugin '{name}': '{exe}' not on PATH"
        if fmt == "parquet":
            import importlib.util
            if importlib.util.find_spec("pyarrow") is None:
                return False, (f"plugin '{name}': parquet exchange needs pyarrow "
                               "(pip install 'grnbench[plugins]')")
        return True, ""
    return probe


def command_method_spec(cfg: dict, source: str = "<inline>") -> MethodSpec:
    """Build a :class:`MethodSpec` that shells out to an external command."""
    name = cfg.get("name")
    if not name:
        raise PluginError(f"{source}: plugin spec needs a 'name'")
    raw_cmd = cfg.get("command")
    if not raw_cmd:
        raise PluginError(f"{source}: plugin '{name}' needs a 'command'")
    cmd = shlex.split(raw_cmd) if isinstance(raw_cmd, str) else [str(c) for c in raw_cmd]
    fmt = str(cfg.get("format", "parquet")).lower()
    if fmt not in {"parquet", "tsv"}:
        raise PluginError(f"{source}: plugin '{name}': format must be parquet or tsv")
    timeout_s = int(cfg.get("timeout_s", 3600))
    extra_env = {str(k): str(v) for k, v in (cfg.get("env") or {}).items()}
    keep_workdir = cfg.get("workdir")

    def runner(adata, tf_list=None, *, seed: int = 42, threads: int = 4,
               **params) -> pd.DataFrame:
        from grnbench.infer.base import _to_frame
        expr = _to_frame(adata)
        ext = "parquet" if fmt == "parquet" else "tsv"

        def _execute(workdir: Path) -> pd.DataFrame:
            expr_p = workdir / f"expr.{ext}"
            tfs_p = workdir / "tfs.csv"
            out_p = workdir / f"edges.{ext}"
            params_p = workdir / "params.json"

            _write_expression(expr, expr_p, fmt)
            pd.DataFrame({"tf": [str(t) for t in (tf_list or [])]}).to_csv(
                tfs_p, index=False)
            params_p.write_text(
                json.dumps({"seed": seed, "threads": threads,
                            **{k: v for k, v in params.items()}}, default=str, indent=2),
                encoding="utf-8")

            subs = {
                "workdir": str(workdir), "expr": str(expr_p), "tfs": str(tfs_p),
                "out": str(out_p), "params": str(params_p),
                "expr_name": expr_p.name, "tfs_name": tfs_p.name,
                "out_name": out_p.name, "params_name": params_p.name,
                "threads": str(threads), "seed": str(seed),
            }
            argv = [tok.format_map(_SafeSubs(subs)) for tok in cmd]
            env = {**os.environ, **extra_env}
            try:
                proc = subprocess.run(argv, capture_output=True, text=True,
                                      timeout=timeout_s, env=env)  # noqa: S603
            except FileNotFoundError as e:
                raise PluginError(f"plugin '{name}': cannot execute {argv[0]!r}: {e}") from e
            except subprocess.TimeoutExpired as e:
                raise PluginError(
                    f"plugin '{name}' exceeded timeout_s={timeout_s}") from e
            if proc.returncode != 0:
                tail = (proc.stderr or proc.stdout or "").strip()[-1500:]
                raise PluginError(
                    f"plugin '{name}' failed (rc={proc.returncode}): {tail}")
            return _read_edges(out_p, fmt)

        if keep_workdir:
            wd = Path(keep_workdir)
            wd.mkdir(parents=True, exist_ok=True)
            return _execute(wd)
        with tempfile.TemporaryDirectory(prefix=f"grnbench_{name}_") as td:
            return _execute(Path(td))

    return MethodSpec(
        name=str(name),
        description=str(cfg.get("description")
                        or f"external command plugin ({cmd[0]})"),
        runner=runner,
        takes="adata",
        origin="command",
        probe=_command_probe(cmd, fmt, str(name)),
    )


class _SafeSubs(dict):
    """format_map helper: unknown ``{placeholders}`` survive untouched."""

    def __missing__(self, key):  # pragma: no cover - trivial
        return "{" + key + "}"


# --------------------------------------------------------------------------- #
# spec files + the user-facing entry point
# --------------------------------------------------------------------------- #

def _read_spec_file(path: Path) -> list[dict]:
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() in {".yaml", ".yml"}:
        try:
            import yaml
        except ImportError as e:
            raise PluginError(
                f"{path.name}: YAML plugin specs need pyyaml (pip install pyyaml); "
                "JSON specs work with no extra dependency"
            ) from e
        data = yaml.safe_load(text)
    else:
        try:
            data = json.loads(text)
        except json.JSONDecodeError as e:
            raise PluginError(f"{path.name}: invalid JSON: {e}") from e
    if isinstance(data, dict) and "methods" in data:
        data = data["methods"]
    if isinstance(data, dict):
        data = [data]
    if not isinstance(data, list):
        raise PluginError(
            f"{path.name}: expected a plugin object or a list of them")
    return data


def load_plugin_file(path: str | Path) -> list[MethodSpec]:
    """Load one or more plugin definitions from a JSON/YAML file."""
    p = Path(path)
    if not p.exists():
        raise PluginError(f"plugin spec file not found: {p}")
    specs = []
    for cfg in _read_spec_file(p):
        if not isinstance(cfg, dict):
            raise PluginError(f"{p.name}: each plugin must be an object, got {cfg!r}")
        if cfg.get("python"):
            obj = resolve_object(str(cfg["python"]), relative_to=p.parent)
            specs.append(python_method_spec(
                obj, name=cfg.get("name"), description=cfg.get("description"),
                origin="python"))
        else:
            specs.append(command_method_spec(cfg, source=str(p)))
    return specs


def load_plugin(spec: str) -> list[MethodSpec]:
    """Resolve one ``--plugin`` argument.

    Accepted: a ``.json``/``.yaml`` spec file, ``package.module:Attribute``, or
    ``path/to/file.py:Attribute``.
    """
    head, attr = split_object_ref(spec)
    if head is not None and head.lower().endswith(".py"):
        return [python_method_spec(resolve_object(spec), origin="python")]
    p = Path(spec)
    if p.exists():
        return load_plugin_file(p)
    if head is not None:
        return [python_method_spec(resolve_object(spec), origin="python")]
    raise PluginError(
        f"--plugin '{spec}': not an existing spec file and not a "
        "'package.module:Attribute' or 'path/to/file.py:Attribute' reference"
    )


#: name -> the --plugin argument that registered it, for idempotent re-registration
_registered_refs: dict[str, str] = {}


def register_plugins(specs: list[str] | None, *, replace: bool = False) -> list[dict]:
    """Load and register every ``--plugin`` argument; return a provenance list.

    Registering the *same* argument twice in one process is a no-op — ``grnbench
    run`` passes ``--plugin`` to both its inference and its evaluation stage, and
    that must not look like a name collision. A *different* plugin claiming an
    already-taken name still fails loudly.
    """
    ensure_plugins_loaded()
    registered: list[dict] = []
    for raw in specs or []:
        for spec in load_plugin(raw):
            if _registered_refs.get(spec.name) != raw:
                try:
                    register_method(spec, replace=replace)
                except MethodError as e:
                    raise PluginError(str(e)) from e
                _registered_refs[spec.name] = raw
            registered.append({"name": spec.name, "origin": spec.origin,
                               "requested_as": raw})
    return registered
