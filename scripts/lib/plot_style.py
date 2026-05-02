"""Common matplotlib/seaborn plotting setup."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


def setup_mpl(style="whitegrid", context="paper", font_scale=1.0):
    """Set up matplotlib for publication-quality figures."""
    matplotlib.use("Agg")
    sns.set_style(style)
    sns.set_context(context, font_scale=font_scale)


def savefig(fig, path, dpi=300, bbox_inches="tight"):
    """Save figure and close it."""
    fig.savefig(path, dpi=dpi, bbox_inches=bbox_inches)
    plt.close(fig)


# Default palette for WT vs mutant comparisons
WT_COLOR = "#1f77b4"       # blue
MUTANT_COLOR = "#d62728"   # red
CGAS_COLOR = "#2ca02c"     # green
TRIM_COLOR = "#ff7f0e"     # orange
