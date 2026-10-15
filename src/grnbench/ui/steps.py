"""The screens of the interface: names, order, groups and the plain-language text of each header.

Pure data, free of Streamlit, so the app, the sidebar and the tests share one list. ``key`` is the label the
sidebar shows and the value kept in ``st.session_state["step"]``.
"""
from __future__ import annotations

from dataclasses import dataclass

QUICK = "Quick run"
DATA = "1 · Data"
PROFILE = "2 · Profile"
PIPELINE = "3 · Pipeline"
COMPARE = "4 · Compare pipelines"
CLEAN = "5 · Clean"
MODELS = "6 · Models"
BENCHMARK = "7 · Benchmark"
EXPLORE_DATA = "Explore data"
EXPLORE_RESULTS = "Explore results"

START, WORKFLOW, EXPLORE = "Start here", "Step by step", "Explore"


@dataclass(frozen=True)
class Step:
    key: str
    group: str
    icon: str            # a Material symbol name
    title: str           # the headline of the page
    lede: str            # one or two sentences under it
    optional: bool = False

    @property
    def short(self) -> str:
        """The name in the menu without its number or the word Explore: ``Pipeline``, ``Benchmark``, ``Results``."""
        name = self.key.split(" · ", 1)[-1]
        return name.removeprefix("Explore ").capitalize() if self.group == EXPLORE else name


STEPS: tuple[Step, ...] = (
    Step(QUICK, START, "bolt", "From raw counts to compared network models",
         "Pick a dataset and some models, then press Run. grnbench reads and cleans the data with a fixed, "
         "well-tested pipeline, runs the models and, if you give a reference network, scores them. "
         "Every step can be redone with full control in the workflow."),
    Step(DATA, WORKFLOW, "upload_file", "Load your data",
         "Any common format (h5ad, 10x, loom, CSV/TSV). README files and download pages are read to recognise "
         "the protocol, the organism and the experimental design."),
    Step(PROFILE, WORKFLOW, "monitoring", "What your data looks like",
         "Size, protocol, batches, cell types and quality, measured from the counts and read from your notes. "
         "These facts drive the cleaning choices in the next step."),
    Step(PIPELINE, WORKFLOW, "tune", "How your data will be cleaned",
         "For each decision the engine picked one option from a cited knowledge base, using your profile. "
         "Look at what it is unsure about and change anything you disagree with."),
    Step(COMPARE, WORKFLOW, "compare_arrows", "Compare cleaning pipelines",
         "Run the recommended pipeline and one-change alternatives on your data and measure them. With a "
         "reference network the pipeline is chosen on one half of it and judged on the other half, so it is "
         "never graded on the links it was picked on.", optional=True),
    Step(CLEAN, WORKFLOW, "cleaning_services", "Clean the data",
         "The pipeline is applied once per cellular context. Your original file is never changed."),
    Step(MODELS, WORKFLOW, "hub", "Run network models",
         "Every model receives exactly the same cleaned data and returns a ranked list of regulator → gene links."),
    Step(BENCHMARK, WORKFLOW, "leaderboard", "Score the models",
         "Compare each model's links with a reference network of known regulation."),
    Step(EXPLORE_DATA, EXPLORE, "scatter_plot", "Explore your data",
         "Look at the cells and genes before and after cleaning."),
    Step(EXPLORE_RESULTS, EXPLORE, "share", "Explore the networks",
         "Browse each model's top links, draw the network, and see how much the models agree."),
)
LABELS = [s.key for s in STEPS]
BY_KEY = {s.key: s for s in STEPS}
GROUPS = (START, WORKFLOW, EXPLORE)
WORKFLOW_KEYS = [s.key for s in STEPS if s.group == WORKFLOW]


def position(key: str) -> tuple[int, int] | None:
    """(index, total) of a step inside the guided workflow, or None for the other screens."""
    return (WORKFLOW_KEYS.index(key) + 1, len(WORKFLOW_KEYS)) if key in WORKFLOW_KEYS else None
