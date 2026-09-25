"""Plot recorded calibration JSON without rerunning or selecting experiments."""
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", type=Path)
    args = parser.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    data = json.loads(args.result.read_text(encoding="utf-8"))
    rows = data["scenarios"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8), sharey=True)
    for y, row in enumerate(rows):
        color = "#256e99" if row["independence_assumption_holds"] else "#bd4d29"
        for ax, key in zip(axes, ["familywise_false_superiority", "simultaneous_coverage"]):
            bounds = row["gate"][key]
            x = bounds["rate"]
            ax.errorbar(x, y, xerr=[[x - bounds["lower"]], [bounds["upper"] - x]],
                        fmt="o", color=color, capsize=3, markersize=6)
    axes[0].set_yticks(range(len(rows)), [r["name"].replace("_", " ") for r in rows])
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Probability of at least one false superiority call")
    axes[1].set_xlabel("Probability all pairwise intervals cover the truth")
    axes[0].axvline(rows[0]["gate"]["error_ceiling"], ls="--", color="#777777", label="gate tolerance")
    axes[1].axvline(rows[0]["gate"]["coverage_floor"], ls="--", color="#777777")
    axes[0].set_xlim(0, max(.15, max(r["gate"]["familywise_false_superiority"]["upper"] for r in rows) + .04))
    axes[1].set_xlim(max(0, min(r["gate"]["simultaneous_coverage"]["lower"] for r in rows) - .04), 1.01)
    for ax in axes:
        ax.grid(axis="x", alpha=.2)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Score-level calibration: independent units and pseudoreplication", fontsize=14)
    reps = rows[0]["gate"]["familywise_false_superiority"]["trials"]
    confidence = rows[0]["gate"]["familywise_false_superiority"]["one_sided_confidence"]
    fig.text(.5, .015, f"{reps:,} experiments per scenario; whiskers: separate one-sided {confidence:.0%} exact Monte Carlo bounds.\n"
             "Blue: independent units. Orange: duplicated units. Dashed lines: prespecified tolerance. Not biological validation.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .08, 1, .94))
    for extension in ["png", "svg"]:
        fig.savefig(args.result.with_name(f"calibration.{extension}"), dpi=170, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
