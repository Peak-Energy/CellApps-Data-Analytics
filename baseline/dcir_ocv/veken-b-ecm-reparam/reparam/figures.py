"""Figures written to outputs/<date>/figures/.
Colours: team datapack = orange (dashed), new-OCV set = blue, team-OCV set = green, new-OCV 2RC set = indigo,
team-OCV 2RC set = magenta (2RC sets dash-dotted). Palette checked for colour-vision separation."""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import load

INK, MUTED, SURF, SHADE = "#0b0b0b", "#52514e", "#fcfcfb", "#efeee9"
NEW, ORIG, TEAMOCV, NEW2, TEAMOCV2 = "#2a78d6", "#eb6834", "#1baf7a", "#5a4b8a", "#d13c8f"
SET_COLORS = {"team datapack": ORIG, "new OCV set": NEW, "team OCV set": TEAMOCV,
              "new OCV 2RC set": NEW2, "team OCV 2RC set": TEAMOCV2}
SET_STYLES = {"team datapack": (0, (3, 2)), "new OCV 2RC set": (0, (5, 1.5, 1, 1.5)), "team OCV 2RC set": (0, (5, 1.5, 1, 1.5))}


def style(name):
    return SET_STYLES.get(name, "-")
PHASE_SHADES = ["#f4f3ef", "#e8eef7", "#f7efe8"]
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK, "axes.grid": True,
                     "grid.color": "#e4e3df", "grid.linewidth": 0.6, "lines.linewidth": 1.4, "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})


def ocv(new, orig, path):
    fig, ax = plt.subplots(3, 3, figsize=(13, 9), sharex=True)
    for j, T in enumerate([15, 25, 45]):
        a, b = new[new.Temperature == T], orig[orig.Temperature == T]
        for i, c in enumerate(["OCV-Charge", "OCV-Discharge", "Average"]):
            ax[i, j].plot(b.SOC, b[c], color=ORIG, ls=(0, (3, 2)), label="team OCV")
            ax[i, j].plot(a.SOC, a[c], color=NEW, label="new OCV")
            tw = ax[i, j].twinx()
            tw.plot(a.SOC, 1e3 * (a[c].to_numpy() - b[c].to_numpy()), color=MUTED, lw=0.9)
            tw.set_ylim(-150, 150); tw.grid(False); tw.tick_params(colors=MUTED, labelsize=7)
            tw.spines["top"].set_visible(False)
            if j == 2:
                tw.set_ylabel("new − team (mV, grey line)", color=MUTED, fontsize=8)
            ax[i, j].set_title(f"{c}, {T} °C", loc="left")
        ax[2, j].set_xlabel("SOC")
    for i in range(3):
        ax[i, 0].set_ylabel("V")
    ax[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def entropy(new, orig, path):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    for j, c in enumerate(["Charge", "Discharge"]):
        ax[j].axvspan(0, 0.1, color=SHADE, lw=0); ax[j].axvspan(0.9, 1, color=SHADE, lw=0)
        ax[j].plot(orig.SOC, 1e3 * orig[c], "o--", color=ORIG, ms=4, label="team datapack")
        ax[j].plot(new.SOC, 1e3 * new[c], "o-", color=NEW, ms=4, label="candidate from the new OCV (not in either set)")
        ax[j].set_title(f"{c} branch (shaded = low confidence)", loc="left"); ax[j].set_xlabel("SOC")
    ax[0].set_ylabel("dU/dT (mV/K)"); ax[0].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def rc(tables, col, ylabel, path):
    """tables: {set name: datapack table}; one panel per temperature."""
    fig, ax = plt.subplots(2, 4, figsize=(15, 7))
    e = col.replace("Mean", "Error")
    for k, T in enumerate(load.R_T_GRID):
        a_ = ax.flat[k]
        measured = T in (15, 25, 45)
        if not measured:
            a_.set_facecolor(SHADE)
        for name, t in tables.items():
            g = t[t.Temperature == T]
            c = SET_COLORS[name]
            a_.fill_between(g.SOC, g[col] - g[e], g[col] + g[e], color=c, alpha=0.12, lw=0)
            a_.plot(g.SOC, g[col], color=c, ls=style(name), label=name)
        a_.set_yscale("log")
        a_.set_title(f"{T} °C — {'measured' if measured else 'Arrhenius-extrapolated'}", loc="left",
                     color=INK if measured else MUTED)
        a_.set_xlabel("SOC")
    ax[0, 0].set_ylabel(ylabel); ax[1, 0].set_ylabel(ylabel); ax[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"{ylabel} (bands = ± Error column; grey panels are extrapolated rows)", x=0.01, ha="left")
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def _shade_phases(ax, t, phase, label=False):
    edges = np.flatnonzero(np.r_[True, phase[1:] != phase[:-1], True])
    for k in range(len(edges) - 1):
        a, b = t[edges[k]], t[edges[k + 1] - 1]
        ax.axvspan(a, b, color=PHASE_SHADES[k % 3], lw=0, zorder=0)
        if label:
            ax.text(0.5 * (a + b), 1.02, phase[edges[k]], transform=ax.get_xaxis_transform(), ha="center",
                    va="bottom", fontsize=8, color=MUTED)


def whole_file(tr, sets, title, path, stats):
    """Measured vs simulated voltage, error and current over the entire file, for each parameter set."""
    t, V, I, ph, first = tr["t_h"], tr["V"], tr["I"], tr["phase"], tr["first"]
    fig, ax = plt.subplots(3, 1, figsize=(15, 9), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1.6, 0.8]))
    for a_ in ax:
        _shade_phases(a_, t, ph, label=a_ is ax[0])
    ax[0].plot(t, V, color=MUTED, lw=2.0, label="measured")
    ax[1].axhspan(-50, 50, color="#dedcd6", alpha=0.5, lw=0)
    for name in sets:
        v, _ = tr[(name, "SOC 1 at full charge")]
        c = SET_COLORS[name]
        ax[0].plot(t, v, color=c, lw=0.9, ls=style(name), label=name)
        e = 1e3 * (v - V); e[first] = np.nan
        ax[1].plot(t, e, color=c, lw=0.7, ls=style(name), label=f"{name}: RMSE {stats[name]:.0f} mV")
    ax[0].set_ylabel("V"); ax[0].legend(frameon=False, fontsize=8, loc="lower left")
    ax[1].set_ylim(-250, 250); ax[1].set_ylabel("model − measured (mV)\nband = ±50 mV")
    ax[1].legend(frameon=False, fontsize=8, loc="lower left")
    ax[2].plot(t, I, color=INK, lw=0.6); ax[2].set_ylabel("current (A)\ndischarge +"); ax[2].set_xlabel("time (h)")
    fig.suptitle(title, x=0.01, ha="left", y=1.0)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def hppc_overview(traces, sets, stats, path):
    """Model − measured over each entire HPPC file (one panel per file), for each parameter set."""
    fig, ax = plt.subplots(len(traces), 1, figsize=(15, 2.3 * len(traces)))
    for a_, (name, tr) in zip(np.atleast_1d(ax), traces.items()):
        t, V, first = tr["t_h"], tr["V"], tr["first"]
        _shade_phases(a_, t, tr["phase"], label=True)
        a_.axhspan(-50, 50, color="#dedcd6", alpha=0.5, lw=0)
        for s in sets:
            v, _ = tr[(s, "SOC 1 at full charge")]
            e = 1e3 * (v - V); e[first] = np.nan
            a_.plot(t, e, color=SET_COLORS[s], lw=0.6, ls=style(s), label=f"{s}: RMSE {stats[name][s]:.0f} mV")
        a_.set_ylim(-250, 250); a_.set_ylabel("mV")
        a_.set_title(name, loc="left", fontsize=9, pad=14)
        a_.legend(frameon=False, fontsize=8, loc="lower left")
    np.atleast_1d(ax)[-1].set_xlabel("time (h)")
    fig.suptitle("HPPC records, model − measured voltage over each whole file (band = ±50 mV)", x=0.01, ha="left")
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def rpt_summary(whole, order, path):
    """Whole-file RMSE and 99th-percentile |error| on the held-out RPT cell, by temperature, one bar per set."""
    temps = [15, 25, 45, 60]
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.8))
    w = 0.8 / len(order)
    for a, col, title in [(ax[0], "rmse_mV", "Whole-file RMSE (mV)"), (ax[1], "p99_mV", "99th percentile |error| (mV)")]:
        for k, s in enumerate(order):
            d = whole[whole.set == s].set_index("T")[col].reindex(temps)
            x = np.arange(len(temps)) + (k - (len(order) - 1) / 2) * w
            bars = a.bar(x, d.to_numpy(), w - 0.02, color=SET_COLORS[s], label=s, edgecolor=SURF, linewidth=1.0)
            for b_, v in zip(bars, d.to_numpy()):
                a.text(b_.get_x() + b_.get_width() / 2, v, f"{v:.0f}" if col == "p99_mV" else f"{v:.1f}",
                       ha="center", va="bottom", fontsize=8, color=INK)
        a.set_xticks(np.arange(len(temps)), [f"{T} °C" + (" (outside fit)" if T == 60 else "") for T in temps])
        a.set_title(title, loc="left"); a.grid(axis="x", visible=False); a.set_axisbelow(True)
    ax[0].legend(frameon=False, ncol=3, loc="upper left", bbox_to_anchor=(0, -0.1), fontsize=8)
    fig.suptitle("Held-out cell CS2D7485, entire RPT files simulated", x=0.01, ha="left", fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)
