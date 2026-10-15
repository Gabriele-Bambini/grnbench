"""Colours and Streamlit options of the interface, in one place.

``grnbench ui`` hands these to Streamlit as command-line options (a user's own ``--theme.*`` option wins), so
the look needs no config file next to the data. There is a light and a dark variant and Streamlit picks the one
the browser asks for; only the accent, the radii, the heading sizes and the chart colours are set.
"""
from __future__ import annotations

import json

ACCENT = "#0F766E"            # teal 700: white text on it is 5.5:1, so it is the fill of buttons and badges in both modes
ACCENT_ON_DARK = "#5EEAD4"    # teal 300: the accent as text or outline on a dark page (the fill would be too dim there)

# Colour-blind safe (Okabe and Ito), the same palette the embedding maps use.
CATEGORICAL = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#999999"]
SEQUENTIAL = ["#F0FDFA", "#CCFBF1", "#99F6E4", "#5EEAD4", "#2DD4BF", "#14B8A6", "#0D9488", "#0F766E", "#115E59", "#134E4A"]

# Shared by both variants.
THEME = {
    "baseRadius": "0.6rem",
    "buttonRadius": "0.5rem",
    "showSidebarBorder": "true",
    "headingFontSizes": json.dumps(["2rem", "1.5rem", "1.2rem"]),
    "headingFontWeights": json.dumps([700, 650, 600]),
    "chartCategoricalColors": json.dumps(CATEGORICAL),
    "chartSequentialColors": json.dumps(SEQUENTIAL),
}
# Per variant. Naming either variant (or ``--theme.base``) in a section of its own is what makes Streamlit follow
# the browser's light or dark setting; a plain ``--theme.primaryColor`` would pin one look for everyone.
LIGHT = {"primaryColor": ACCENT, "greenColor": "#15803D", "orangeColor": "#B45309", "redColor": "#B91C1C",
         "blueColor": "#1D4ED8"}
DARK = {"primaryColor": ACCENT}
# This is a tool for the data of a lab: no usage statistics, and no "Deploy" button.
SERVER = {"browser.gatherUsageStats": "false", "client.toolbarMode": "minimal"}


def streamlit_flags(given=()) -> list[str]:
    """Command-line options for ``streamlit run`` that the caller has not already given.

    ``--theme.primaryColor`` (without a variant) replaces our colour in both variants; naming one variant replaces
    that variant's only. ``--theme.base`` means the caller brings their own theme: none of our colours is added.
    """
    taken = {str(a).split("=", 1)[0].lstrip("-") for a in given if str(a).startswith("--")}
    options = [(f"theme.{k}", v) for k, v in THEME.items()]
    if "theme.base" not in taken:
        for variant, colours in (("light", LIGHT), ("dark", DARK)):
            options += [(f"theme.{variant}.{k}", v) for k, v in colours.items() if f"theme.{k}" not in taken]
    options += list(SERVER.items())
    flags: list[str] = []
    for name, value in options:
        if name not in taken:
            flags += [f"--{name}", value]
    return flags
