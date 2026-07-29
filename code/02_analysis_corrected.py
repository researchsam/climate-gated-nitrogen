from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from linearmodels.panel import PanelOLS
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.model_selection import GroupKFold

warnings.filterwarnings("ignore")

REPO_ROOT = Path(__file__).resolve().parents[1]

parser = argparse.ArgumentParser(
    description="Reproduce the climate-gated nitrogen analysis and manuscript figures."
)
parser.add_argument(
    "--data",
    type=Path,
    default=REPO_ROOT / "data" / "county_year_model_frame.csv",
    help="County-year analytical panel.",
)
parser.add_argument(
    "--output",
    type=Path,
    default=REPO_ROOT / "results",
    help="Directory for figures and result tables.",
)
parser.add_argument(
    "--bootstrap-replicates",
    type=int,
    default=300,
    help="Number of coefficient-refitting bootstrap replicates.",
)
args = parser.parse_args()

DATA = args.data.resolve()
OUT = args.output.resolve()
OUT.mkdir(parents=True, exist_ok=True)

SEED = 42
CUT_PCT = 15.0
WET_FRAC = 1 / 3
CAP = 0.60
BOOT_B = args.bootstrap_replicates

BLUE = "#216B8C"
BLUE_LIGHT = "#BBD7E5"
CORAL = "#D45345"
GRAY = "#777777"
GREEN = "#3F8F5F"


def gof(obs, pred):
    obs, pred = np.asarray(obs, float), np.asarray(pred, float)
    m = np.isfinite(obs) & np.isfinite(pred)
    obs, pred = obs[m], pred[m]
    nse = 1 - np.sum((obs - pred) ** 2) / np.sum((obs - obs.mean()) ** 2)
    r = np.corrcoef(obs, pred)[0, 1]
    kge = 1 - np.sqrt(
        (r - 1) ** 2
        + (pred.std() / obs.std() - 1) ** 2
        + (pred.mean() / obs.mean() - 1) ** 2
    )
    rmse = np.sqrt(np.mean((obs - pred) ** 2))
    return dict(NSE=float(nse), KGE=float(kge), RMSE=float(rmse))


raw = pd.read_csv(DATA)
df = raw.rename(
    columns={
        "loading_kgha_lag1": "loading_lag1",
        "loading_kgha_lag2": "loading_lag2",
        "loading_kgha_lag3": "loading_lag3",
        "precipmm_lag1": "precip_lag1",
        "precipmm_lag2": "precip_lag2",
        "precipmm_lag3": "precip_lag3",
    }
).copy()
df["NS_precip"] = df["NS"] * df["precipmm"]
df["area_ha"] = df["county_area_ha"]
df["Value"] = df["tile_ha"]
df["rowcrop_frac"] = (df["corng_ha"] + df["soy_ha"]) / df["county_area_ha"]
df["rowcrop_precip"] = df["rowcrop_frac"] * df["precipmm"]
df = df.dropna(
    subset=[
        "CountyName",
        "Year",
        "loading_kgha",
        "NS",
        "precipmm",
        "NS_precip",
        "loading_lag1",
        "area_ha",
    ]
).copy()
df["Year"] = df["Year"].astype(int)


def fit_fe(frame, dep="loading_kgha", lag_col="loading_lag1", extra=None):
    cols = ["NS", "precipmm", "NS_precip", lag_col] + (extra or [])
    p = frame.dropna(subset=[dep] + cols).set_index(["CountyName", "Year"]).sort_index()
    exog = p[cols].assign(const=1.0)
    return PanelOLS(
        p[dep],
        exog,
        entity_effects=True,
        time_effects=True,
        drop_absorbed=True,
    ).fit(cov_type="kernel", kernel="bartlett", bandwidth=3)


fe = fit_fe(df)
bNS = float(fe.params["NS"])
bNSP = float(fe.params["NS_precip"])
zero_cross = -bNS / bNSP


def leakfree_z(frame, train_idx, col="precipmm"):
    tr = frame.loc[train_idx]
    mu = tr.groupby("CountyName")[col].mean()
    sd = tr.groupby("CountyName")[col].std()
    gmu, gsd = tr[col].mean(), tr[col].std()
    m = frame["CountyName"].map(mu).fillna(gmu)
    s = frame["CountyName"].map(sd).replace(0, np.nan).fillna(gsd)
    return (frame[col] - m) / s


def rf_new(seed=SEED):
    return RandomForestRegressor(
        n_estimators=400,
        min_samples_leaf=3,
        random_state=seed,
        n_jobs=1,
    )


d = df.reset_index(drop=True).copy()
y = d["loading_kgha"].to_numpy()
rf_base = ["loading_lag1", "precipmm", "Value", "precip_lag1", "NS", "NS_lag1"]
yr_cut = d["Year"].quantile(0.75)
tr = d.index[d["Year"] <= yr_cut]
te = d.index[d["Year"] > yr_cut]
d["precip_z"] = leakfree_z(d, tr)
features = rf_base + ["precip_z"]
rf = rf_new().fit(d.loc[tr, features], y[tr])
pred_te = rf.predict(d.loc[te, features])
rf_temporal = gof(y[te], pred_te)

gkf = GroupKFold(n_splits=5)
group_rows = []
for k, (a, b) in enumerate(gkf.split(d, y, groups=d["CountyName"]), 1):
    dz = d.copy()
    dz["precip_z"] = leakfree_z(dz, d.index[a])
    model = rf_new(SEED + k).fit(dz.iloc[a][features], y[a])
    group_rows.append(gof(y[b], model.predict(dz.iloc[b][features])))
rf_group = pd.DataFrame(group_rows)

perm = permutation_importance(
    rf,
    d.loc[te, features],
    y[te],
    n_repeats=30,
    random_state=SEED,
    n_jobs=1,
)
importance = (
    pd.DataFrame(
        {
            "feature": features,
            "importance": perm.importances_mean,
            "sd": perm.importances_std,
        }
    )
    .sort_values("importance", ascending=False)
    .reset_index(drop=True)
)


def build_calendar(frame, b0, b1, clip_response=True):
    x = frame[
        ["CountyName", "Year", "NS", "precipmm", "loading_kgha", "area_ha"]
    ].dropna().copy()
    x["slope_raw"] = b0 + b1 * x["precipmm"]
    x["slope_row"] = (
        x["slope_raw"].clip(lower=0.0) if clip_response else x["slope_raw"]
    )
    # Both scenarios act only on positive, physically reducible surplus.
    x["NS_reducible"] = x["NS"].clip(lower=0.0)
    x["dkg_ha_CAL"] = x["slope_row"] * (
        x["NS_reducible"] * (-CUT_PCT / 100.0)
    )
    x["tons_CAL"] = x["dkg_ha_CAL"] * x["area_ha"] / 1000.0
    return x


def apply_wet(frame, cap=CAP, wet_frac=WET_FRAC, miss=0.0, seed=SEED):
    x = frame.copy()
    x["wet_rank"] = x.groupby("CountyName")["precipmm"].rank(pct=True)
    sel = x["wet_rank"] >= (1 - wet_frac)
    if miss > 0:
        rng = np.random.default_rng(seed)

        def flip(g):
            idx = g.index
            s = sel.loc[idx].to_numpy().copy()
            nt, nf = int(miss * s.sum()), int(miss * (~s).sum())
            if nt:
                s[
                    rng.choice(np.where(s)[0], min(nt, s.sum()), replace=False)
                ] = False
            if nf:
                s[
                    rng.choice(np.where(~s)[0], min(nf, (~s).sum()), replace=False)
                ] = True
            return pd.Series(s, index=idx)

        sel = x.groupby("CountyName", group_keys=False).apply(flip).astype(bool)
    x["NS_reducible"] = x["NS"].clip(lower=0.0)
    all_ns = x.groupby("CountyName")["NS_reducible"].sum().rename("all")
    wet_ns = x.loc[sel].groupby("CountyName")["NS_reducible"].sum().rename("sel")
    effort = pd.concat([all_ns, wet_ns], axis=1).fillna(0.0)
    effort["cut_wet"] = (
        -np.minimum(
            cap,
            (CUT_PCT / 100.0)
            * effort["all"]
            / effort["sel"].replace(0, np.nan),
        )
    ).fillna(-CUT_PCT / 100.0)
    x["cut_wet"] = x["CountyName"].map(effort["cut_wet"])
    x["dkg_ha_WET"] = 0.0
    x.loc[sel, "dkg_ha_WET"] = x.loc[sel, "slope_row"] * (
        x.loc[sel, "NS_reducible"] * x.loc[sel, "cut_wet"]
    )
    x["tons_WET"] = x["dkg_ha_WET"] * x["area_ha"] / 1000.0
    return x


def annualized(x):
    years = x.groupby("Year")[["tons_CAL", "tons_WET"]].sum()
    cal, wet = -years["tons_CAL"].mean(), -years["tons_WET"].mean()
    return dict(
        calendar=float(cal),
        wet=float(wet),
        gain=float(wet - cal),
        pct=float(100 * (wet - cal) / cal) if cal else np.nan,
    )


scenario_rows = []
for clip_response, label in [
    (True, "Nonnegative marginal response"),
    (False, "Model-implied unconstrained response"),
]:
    cal = build_calendar(df, bNS, bNSP, clip_response=clip_response)
    row = annualized(apply_wet(cal))
    row["scenario"] = label
    scenario_rows.append(row)
scenarios = pd.DataFrame(scenario_rows)
primary = scenarios.iloc[0].to_dict()


def bootstrap_refit(frame, B=BOOT_B, seed=123):
    rng = np.random.default_rng(seed)
    counties = frame["CountyName"].unique()
    gains, betas, zeroes = [], [], []
    for _ in range(B):
        sample = rng.choice(counties, size=len(counties), replace=True)
        sub = pd.concat(
            [
                frame[frame.CountyName == county].assign(
                    CountyName=f"{county}__{k}"
                )
                for k, county in enumerate(sample)
            ],
            ignore_index=True,
        )
        fitted = fit_fe(sub)
        b0 = float(fitted.params["NS"])
        b1 = float(fitted.params["NS_precip"])
        gains.append(
            annualized(
                apply_wet(build_calendar(sub, b0, b1, clip_response=True))
            )["gain"]
        )
        betas.append(b1)
        zeroes.append(-b0 / b1)
    return np.asarray(gains), np.asarray(betas), np.asarray(zeroes)


boot_gain, boot_beta, boot_zero = bootstrap_refit(df)
boot_ci = np.percentile(boot_gain, [2.5, 97.5])

cal_primary = build_calendar(df, bNS, bNSP, clip_response=True)
miss20 = annualized(apply_wet(cal_primary, miss=0.20))["gain"]
late = df[df["Year"] > df["Year"].quantile(0.6)]
holdout = annualized(
    apply_wet(build_calendar(late, bNS, bNSP, clip_response=True))
)["gain"]


def apply_mask(frame, mask, cap=CAP):
    x = frame.copy()
    sel = mask.reindex(x.index).fillna(False).astype(bool)
    x["NS_reducible"] = x["NS"].clip(lower=0.0)
    all_ns = x.groupby("CountyName")["NS_reducible"].sum().rename("all")
    wet_ns = x.loc[sel].groupby("CountyName")["NS_reducible"].sum().rename("sel")
    effort = pd.concat([all_ns, wet_ns], axis=1).fillna(0.0)
    effort["cut_wet"] = (
        -np.minimum(
            cap,
            (CUT_PCT / 100.0)
            * effort["all"]
            / effort["sel"].replace(0, np.nan),
        )
    ).fillna(-CUT_PCT / 100.0)
    x["cut_wet"] = x["CountyName"].map(effort["cut_wet"])
    x["dkg_ha_WET"] = 0.0
    x.loc[sel, "dkg_ha_WET"] = x.loc[sel, "slope_row"] * (
        x.loc[sel, "NS_reducible"] * x.loc[sel, "cut_wet"]
    )
    x["tons_WET"] = x["dkg_ha_WET"] * x["area_ha"] / 1000.0
    return x


rng = np.random.default_rng(SEED)
n_wet = max(1, int(round(WET_FRAC * df.Year.nunique())))
idx_by_county = {
    county: cal_primary.index[cal_primary.CountyName == county].values
    for county in cal_primary.CountyName.unique()
}
placebo = []
for _ in range(200):
    mask = pd.Series(False, index=cal_primary.index)
    for idxs in idx_by_county.values():
        mask.loc[rng.choice(idxs, size=min(n_wet, len(idxs)), replace=False)] = True
    placebo.append(annualized(apply_mask(cal_primary, mask))["gain"])
placebo = np.asarray(placebo)

# Dynamic three-lag specification.
dyn = fit_fe(
    df,
    extra=["loading_lag2", "loading_lag3"],
)
lag_joint = dyn.wald_test(
    formula="loading_lag2 = 0, loading_lag3 = 0"
)

# Leave-one-county-out stability.
loo_rows = []
for county in sorted(df["CountyName"].unique()):
    sub = df[df["CountyName"] != county]
    fitted = fit_fe(sub)
    b0, b1 = float(fitted.params["NS"]), float(fitted.params["NS_precip"])
    gain = annualized(apply_wet(build_calendar(sub, b0, b1)))["gain"]
    loo_rows.append(
        dict(
            omitted=county,
            beta=b1,
            zero_cross=-b0 / b1,
            gain=gain,
        )
    )
loo = pd.DataFrame(loo_rows)

# Remove exactly ten counties with the largest composite mean-load/positive-surplus rank.
county_summary = df.groupby("CountyName").agg(
    mean_load=("loading_kgha", "mean"),
    mean_positive_surplus=("NS", lambda s: s.clip(lower=0).mean()),
)
county_summary["leverage_score"] = (
    county_summary["mean_load"].rank(pct=True)
    + county_summary["mean_positive_surplus"].rank(pct=True)
)
top10 = county_summary.nlargest(10, "leverage_score").index.tolist()
sub10 = df[~df["CountyName"].isin(top10)]
fe10 = fit_fe(sub10)
b10_ns = float(fe10.params["NS"])
b10_nsp = float(fe10.params["NS_precip"])
gain10 = annualized(
    apply_wet(build_calendar(sub10, b10_ns, b10_nsp))
)["gain"]

# Centering and nonlinear checks.
centered = df.copy()
centered["NS_c"] = centered["NS"] - centered["NS"].mean()
centered["P_c"] = centered["precipmm"] - centered["precipmm"].mean()
centered["NSP_c"] = centered["NS_c"] * centered["P_c"]
pcenter = centered.set_index(["CountyName", "Year"]).sort_index()
fe_center = PanelOLS(
    pcenter["loading_kgha"],
    pcenter[["NS_c", "P_c", "NSP_c", "loading_lag1"]].assign(const=1.0),
    entity_effects=True,
    time_effects=True,
).fit(cov_type="kernel", kernel="bartlett", bandwidth=3)

nonlinear = df.copy()
nonlinear["NS2"] = nonlinear["NS"] ** 2
nonlinear["P2"] = nonlinear["precipmm"] ** 2
fe_nonlinear = fit_fe(nonlinear, extra=["NS2", "P2"])

# Row-crop sensitivity.
tile_row_corr = float(df[["tile_frac", "rowcrop_frac"]].corr().iloc[0, 1])
fe_rowcrop = fit_fe(
    df,
    extra=["rowcrop_frac", "rowcrop_precip"],
)

dr = d.copy()
dr["rowcrop_frac"] = df.reset_index(drop=True)["rowcrop_frac"]
rf_row_features = features + ["rowcrop_frac"]
rf_row = rf_new().fit(dr.loc[tr, rf_row_features], y[tr])
perm_row = permutation_importance(
    rf_row,
    dr.loc[te, rf_row_features],
    y[te],
    n_repeats=30,
    random_state=SEED,
    n_jobs=1,
)
importance_row = (
    pd.DataFrame(
        {
            "feature": rf_row_features,
            "importance": perm_row.importances_mean,
            "sd": perm_row.importances_std,
        }
    )
    .sort_values("importance", ascending=False)
    .reset_index(drop=True)
)

# Coverage and load-denominator sensitivities.
coverage_rows = []
for threshold in [0.25, 0.50, 0.75]:
    sub = df[df["mean_coverage_frac"] >= threshold]
    fitted = fit_fe(sub)
    coverage_rows.append(
        dict(
            threshold=threshold,
            n=int(fitted.nobs),
            counties=int(sub.CountyName.nunique()),
            beta=float(fitted.params["NS_precip"]),
            p=float(fitted.pvalues["NS_precip"]),
            zero_cross=float(
                -fitted.params["NS"] / fitted.params["NS_precip"]
            ),
        )
    )
coverage = pd.DataFrame(coverage_rows)

fullarea = df.sort_values(["CountyName", "Year"]).copy()
fullarea["fullarea_lag1"] = fullarea.groupby("CountyName")[
    "loading_kgha_fullarea"
].shift(1)
fe_fullarea = fit_fe(
    fullarea,
    dep="loading_kgha_fullarea",
    lag_col="fullarea_lag1",
)

fwmc = df.sort_values(["CountyName", "Year"]).copy()
fwmc["fwmc_lag1"] = fwmc.groupby("CountyName")["FWMC_mgL"].shift(1)
fe_fwmc = fit_fe(fwmc, dep="FWMC_mgL", lag_col="fwmc_lag1")

# Manuscript figures.
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 180,
        "savefig.dpi": 300,
    }
)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
ax = axes[0]
ax.scatter(y[te], pred_te, s=18, alpha=0.62, color=BLUE, edgecolor="none")
lim = max(float(np.nanmax(y[te])), float(np.nanmax(pred_te)))
ax.plot([0, lim], [0, lim], "--", lw=1.2, color=GRAY)
ax.set_xlim(0, lim)
ax.set_ylim(0, lim)
ax.set_xlabel("Estimated load (kg N ha$^{-1}$ yr$^{-1}$)")
ax.set_ylabel("RF-predicted load (kg N ha$^{-1}$ yr$^{-1}$)")
ax.set_title(
    f"(a) Temporal hold-out: NSE={rf_temporal['NSE']:.2f}, "
    f"KGE={rf_temporal['KGE']:.2f}"
)

label_map = {
    "precip_z": "Precipitation anomaly",
    "Value": "Tile-drained area",
    "loading_lag1": "Antecedent load",
    "precip_lag1": "Antecedent precipitation",
    "precipmm": "Annual precipitation",
    "NS": "Nitrogen surplus",
    "NS_lag1": "Antecedent surplus",
}
plot_imp = importance.sort_values("importance")
ax = axes[1]
ax.barh(
    [label_map[x] for x in plot_imp["feature"]],
    plot_imp["importance"],
    xerr=plot_imp["sd"],
    color=BLUE,
    alpha=0.95,
    error_kw=dict(ecolor=GRAY, capsize=2, lw=1),
)
ax.axvline(0, color="#BBBBBB", lw=0.8)
ax.set_xlabel("Permutation importance")
ax.set_title("(b) Leakage-free predictor importance")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Fig1_RF_corrected.{ext}", bbox_inches="tight")
plt.close(fig)

p_grid = np.linspace(
    df["precipmm"].quantile(0.01),
    df["precipmm"].quantile(0.99),
    300,
)
marginal = bNS + bNSP * p_grid
cov = fe.cov
variance = (
    cov.loc["NS", "NS"]
    + p_grid**2 * cov.loc["NS_precip", "NS_precip"]
    + 2 * p_grid * cov.loc["NS", "NS_precip"]
)
se = np.sqrt(np.maximum(variance, 0))

fig, ax = plt.subplots(figsize=(8.2, 5.0))
ax.fill_between(
    p_grid,
    marginal - 1.96 * se,
    marginal + 1.96 * se,
    color=BLUE_LIGHT,
    alpha=0.65,
    label="95% confidence interval",
)
ax.plot(p_grid, marginal, color=BLUE, lw=2.5, label="Marginal NS effect")
ax.axhline(0, color=GRAY, ls="--", lw=1)
ax.axvline(zero_cross, color=CORAL, ls=":", lw=1.5)
ax.annotate(
    f"Zero-crossing\n{zero_cross:.0f} mm",
    xy=(zero_cross, 0),
    xytext=(zero_cross + 55, marginal.max() * 0.62),
    color=CORAL,
    arrowprops=dict(arrowstyle="-", color=CORAL, lw=1),
)
ax.set_xlabel("Annual precipitation (mm)")
ax.set_ylabel(
    "Marginal effect of N surplus on load\n"
    "$\\partial$Load / $\\partial$NS"
)
ax.set_title("Climate gating: surplus–load association strengthens in wetter years")
ax.legend(frameon=False, loc="lower right")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Fig2_marginal_effect_corrected.{ext}", bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), sharey=True)
for ax, row, title in zip(
    axes,
    scenarios.to_dict("records"),
    [
        "(a) Primary: nonnegative marginal response",
        "(b) Sensitivity: model-implied response",
    ],
):
    vals = [row["calendar"], row["wet"]]
    bars = ax.bar(
        ["Calendar\nuniform effort", "Wet-year\nsame total effort"],
        vals,
        color=["#9A9A9A", BLUE],
        width=0.62,
    )
    for bar, val in zip(bars, vals):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            val + max(vals) * 0.035,
            f"{val:,.0f}",
            ha="center",
            va="bottom",
            fontweight="bold",
        )
    ax.set_title(title)
    ax.set_ylim(0, 2150)
    ax.text(
        0.5,
        2025,
        f"Timing dividend: {row['gain']:,.0f} Mg N yr$^{{-1}}$",
        ha="center",
        color=CORAL,
        fontweight="bold",
    )
axes[0].set_ylabel("Modeled annual N-load removal (Mg N yr$^{-1}$)")
axes[0].text(
    0.5,
    1880,
    f"Refit-bootstrap 95% CI: {boot_ci[0]:,.0f}–{boot_ci[1]:,.0f}",
    ha="center",
    color=CORAL,
)
fig.suptitle("Equal total surplus-reduction effort, different timing", y=1.02)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Fig3_timing_corrected.{ext}", bbox_inches="tight")
plt.close(fig)

labels = [
    "Base",
    "Bootstrap\nmedian",
    "20%\nmisclassification",
    "Temporal\nholdout",
    "Random-year\nplacebo",
]
values = [
    primary["gain"],
    float(np.median(boot_gain)),
    miss20,
    holdout,
    float(placebo.mean()),
]
fig, ax = plt.subplots(figsize=(9, 5.1))
bars = ax.bar(labels, values, color=[BLUE, BLUE, GREEN, GREEN, "#B0B0B0"], width=0.62)
err_low = np.median(boot_gain) - boot_ci[0]
err_high = boot_ci[1] - np.median(boot_gain)
ax.errorbar(
    1,
    np.median(boot_gain),
    yerr=np.array([[err_low], [err_high]]),
    fmt="none",
    ecolor=CORAL,
    capsize=5,
    lw=1.8,
)
for bar, val in zip(bars, values):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        val + 35,
        f"{val:,.0f}",
        ha="center",
        va="bottom",
        fontweight="bold",
    )
ax.axhline(0, color="#BBBBBB", lw=0.8)
ax.set_ylabel("Added modeled removal from wet targeting (Mg N yr$^{-1}$)")
ax.set_title("Timing dividend across robustness tests")
ax.text(
    4,
    175,
    f"Placebo 95th percentile: {np.percentile(placebo, 95):,.0f}",
    ha="center",
    color=GRAY,
    fontsize=9,
)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Fig4_robustness_corrected.{ext}", bbox_inches="tight")
plt.close(fig)

results = {
    "sample": {
        "rows": len(df),
        "counties": int(df.CountyName.nunique()),
        "years": [int(df.Year.min()), int(df.Year.max())],
        "mean_load": float(df.loading_kgha.mean()),
        "range_load": [float(df.loading_kgha.min()), float(df.loading_kgha.max())],
        "mean_ns": float(df.NS.mean()),
        "range_ns": [float(df.NS.min()), float(df.NS.max())],
        "mean_precip": float(df.precipmm.mean()),
        "range_precip": [float(df.precipmm.min()), float(df.precipmm.max())],
    },
    "fixed_effects": {
        "params": {k: float(v) for k, v in fe.params.items()},
        "std_errors": {k: float(v) for k, v in fe.std_errors.items()},
        "pvalues": {k: float(v) for k, v in fe.pvalues.items()},
        "within_r2": float(fe.rsquared_within),
        "zero_cross_mm": float(zero_cross),
    },
    "rf_temporal": rf_temporal,
    "rf_group_mean": rf_group.mean().to_dict(),
    "rf_group_sd": rf_group.std().to_dict(),
    "dynamic": {
        "lag1": float(dyn.params["loading_lag1"]),
        "lag1_p": float(dyn.pvalues["loading_lag1"]),
        "lag2": float(dyn.params["loading_lag2"]),
        "lag2_p": float(dyn.pvalues["loading_lag2"]),
        "lag3": float(dyn.params["loading_lag3"]),
        "lag3_p": float(dyn.pvalues["loading_lag3"]),
        "joint_stat": float(lag_joint.stat),
        "joint_p": float(lag_joint.pval),
        "interaction": float(dyn.params["NS_precip"]),
        "interaction_p": float(dyn.pvalues["NS_precip"]),
    },
    "scenarios": scenarios.to_dict("records"),
    "bootstrap": {
        "n": len(boot_gain),
        "gain_median": float(np.median(boot_gain)),
        "gain_ci": boot_ci.tolist(),
        "beta_ci": np.percentile(boot_beta, [2.5, 97.5]).tolist(),
        "zero_cross_ci": np.percentile(boot_zero, [2.5, 97.5]).tolist(),
    },
    "robustness": {
        "miss20": float(miss20),
        "holdout": float(holdout),
        "placebo_mean": float(placebo.mean()),
        "placebo_p95": float(np.percentile(placebo, 95)),
    },
    "leave_one_out": {
        "beta_range": [float(loo.beta.min()), float(loo.beta.max())],
        "zero_cross_range": [
            float(loo.zero_cross.min()),
            float(loo.zero_cross.max()),
        ],
        "gain_range": [float(loo.gain.min()), float(loo.gain.max())],
        "max_abs_gain_shift": float((loo.gain - primary["gain"]).abs().max()),
    },
    "top10_removed": {
        "counties": top10,
        "beta": b10_nsp,
        "zero_cross": float(-b10_ns / b10_nsp),
        "gain": gain10,
    },
    "rowcrop": {
        "tile_rowcrop_corr": tile_row_corr,
        "ns_precip_beta": float(fe_rowcrop.params["NS_precip"]),
        "ns_precip_p": float(fe_rowcrop.pvalues["NS_precip"]),
        "rowcrop_precip_beta": float(fe_rowcrop.params["rowcrop_precip"]),
        "rowcrop_precip_p": float(fe_rowcrop.pvalues["rowcrop_precip"]),
    },
    "nonlinear": {
        "interaction_beta": float(fe_nonlinear.params["NS_precip"]),
        "interaction_p": float(fe_nonlinear.pvalues["NS_precip"]),
        "ns2_p": float(fe_nonlinear.pvalues["NS2"]),
        "p2_p": float(fe_nonlinear.pvalues["P2"]),
    },
    "fullarea": {
        "beta": float(fe_fullarea.params["NS_precip"]),
        "p": float(fe_fullarea.pvalues["NS_precip"]),
        "zero_cross": float(
            -fe_fullarea.params["NS"] / fe_fullarea.params["NS_precip"]
        ),
    },
    "fwmc": {
        "beta": float(fe_fwmc.params["NS_precip"]),
        "p": float(fe_fwmc.pvalues["NS_precip"]),
        "within_r2": float(fe_fwmc.rsquared_within),
    },
}

(OUT / "final_results.json").write_text(
    json.dumps(results, indent=2),
    encoding="utf-8",
)
importance.to_csv(OUT / "rf_permutation_importance.csv", index=False)
importance_row.to_csv(OUT / "rf_rowcrop_importance.csv", index=False)
rf_group.to_csv(OUT / "rf_grouped_cv.csv", index=False)
scenarios.to_csv(OUT / "timing_scenarios.csv", index=False)
coverage.to_csv(OUT / "coverage_sensitivity.csv", index=False)
loo.to_csv(OUT / "leave_one_county_out.csv", index=False)

print(json.dumps(results, indent=2))
