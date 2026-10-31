"""What each registered model is, for the interface: family, rough speed, who ships it.

Pure data and two small functions, free of Streamlit, so the model picker is testable.
Speeds are indicative (about 1,000 genes x 1,500 cells on a 4-core machine) and only meant
to tell a user which models to leave out of a first quick look; a method this table does not
know (a plugin, a model registered later) is listed as "unknown speed" and never hidden.
"""
from __future__ import annotations

import pandas as pd

# method -> (family, speed class). 0 = seconds, 1 = about a minute, 2 = minutes or more.
MODEL_INFO: dict[str, tuple[str, int]] = {
    "pearson_signed": ("Correlation", 0), "spearman_signed": ("Correlation", 0),
    "bicor": ("Correlation", 0), "propr_rho": ("Correlation", 0), "wgcna_tom": ("Correlation", 0),
    "mutual_info": ("Information theory", 0), "mrnet": ("Information theory", 0),
    "aracne_dpi": ("Information theory", 0), "clr": ("Information theory", 1),
    "aracne3": ("Information theory", 2),
    "ppcor": ("Partial correlation", 1), "pcnet": ("Partial correlation", 1),
    "lasso": ("Regression", 0), "bbsr": ("Regression", 0), "tigress": ("Regression", 1),
    "random_forest": ("Trees", 2), "extra_trees": ("Trees", 2), "gradient_boosting": ("Trees", 2),
    "genie3": ("Trees", 2), "grnboost2": ("Trees", 2),
    "regdiffusion": ("Deep learning", 2),
    "motif_prior": ("Sequence prior", 0),
    "pearson_motif": ("Model + motif prior", 0), "grnboost2_motif": ("Model + motif prior", 2),
    "regdiffusion_motif": ("Model + motif prior", 2),
    "pearson_cistarget": ("Model + motif prior", 0), "grnboost2_cistarget": ("Model + motif prior", 2),
    "community": ("Ensemble", 1),
    "random_baseline": ("Negative control", 0),
}
# Per-target L1 regression: its class above holds for a short TF list (the demo's 30; lasso 3 s). With no TF
# list every gene is a regulator and lasso takes ~30 min at 1,400 cells x 1,000 genes (measured).
REGULATOR_BOUND = {"lasso", "tigress", "community"}      # community runs lasso as a member
SPEED_LABEL = {0: "seconds", 1: "about a minute", 2: "minutes or more", None: "unknown speed"}
PRESETS = ("Fast (seconds)", "Fast, including extras", "Recommended by the engine", "Everything available", "None")
TIERS = ("default", "extra", "yours")
# A negative control (D-037) is ticked only by hand: no preset puts it among the models it compares.
CONTROL_FAMILY = "Negative control"


def speed_of(name: str) -> int | None:
    info = MODEL_INFO.get(name)
    return None if info is None else info[1]


def tier_of(entry: dict, user_methods=()) -> str:
    """default (the engine's default set), extra (registered but off by default), yours (a plugin).

    Same vocabulary as ``list_methods()["default"]``: ``--methods available`` runs the default tier.
    """
    if entry["name"] in set(user_methods) or entry.get("origin") not in (None, "builtin"):
        return "yours"
    return "default" if entry.get("default", True) else "extra"


def model_table(methods: list[dict], recommended=(), user_methods=(), *, all_regulators: bool = False) -> pd.DataFrame:
    """One row per registered method, ready for a data editor with a ``run`` checkbox."""
    recommended = set(recommended)
    rows = []
    for m in methods:
        family, speed = MODEL_INFO.get(m["name"], ("Your model" if tier_of(m, user_methods) == "yours" else "Other", None))
        if all_regulators and m["name"] in REGULATOR_BOUND:
            speed = 2
        rows.append({"run": bool(m["available"] and m["name"] in recommended),
                     "method": m["name"], "family": family, "tier": tier_of(m, user_methods),
                     "speed": SPEED_LABEL[speed], "available": bool(m["available"]),
                     "what it is": m.get("description", ""),
                     "why not available": m.get("unavailable_reason") or ""})
    table = pd.DataFrame(rows)
    order = {t: i for i, t in enumerate(TIERS)}
    return table.assign(_o=table["tier"].map(order)).sort_values(
        ["_o", "family", "method"], kind="stable").drop(columns="_o").reset_index(drop=True)


def preset_selection(table: pd.DataFrame, preset: str, recommended=()) -> list[str]:
    """Method names a preset ticks (always only the available ones)."""
    usable = table[table["available"]]
    if preset in ("Fast (seconds)", "Fast, including extras"):
        keep = usable[usable["speed"] == SPEED_LABEL[0]]
        if preset == "Fast (seconds)":
            keep = keep[keep["tier"] == "default"]
        return keep.loc[~keep["family"].isin(["Model + motif prior", "Sequence prior", CONTROL_FAMILY]), "method"].tolist()
    if preset == "Recommended by the engine":
        return [m for m in usable["method"] if m in set(recommended)]
    if preset == "Everything available":
        return usable.loc[usable["family"] != CONTROL_FAMILY, "method"].tolist()
    return []
