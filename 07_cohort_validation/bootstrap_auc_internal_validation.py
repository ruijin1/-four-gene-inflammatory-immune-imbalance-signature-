# -*- coding: utf-8 -*-
# ==========================================================================
#  Four-gene inflammatory signature (CD3G/PRTN3/CD6/THBS1) for MS-BPD
#  File   : 07_cohort_validation/bootstrap_auc_internal_validation.py
#  Status : ORIGINAL (as used in the analysis)
#  Paths  : the authors' working copy used absolute paths on a local Windows
#           machine. They have been replaced here by P("<relative path>")
#           calls. Set BPD_REPO to the repository root (or run the script from
#           the repository root) and place inputs as described in docs/DATA.md.
#  Note   : the statistical code itself is unchanged.
# ==========================================================================
import os as _os
def P(rel):
    """Resolve a repository-relative path (placeholder for local absolute paths)."""
    root = _os.environ.get("BPD_REPO") or _os.path.abspath(
        _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))
    return _os.path.join(root, *rel.replace("\\", "/").split("/"))

"""
Training-set four-gene signature: 1000-replicate bootstrap internal validation (optimism-corrected AUC + 95% CI)
Data: 预测值-系数.xlsx / sheet 'glm模型'  (0 = normal, 1 = MS-BPD)
"""
import numpy as np, pandas as pd, json
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve

RNG_SEED = 20260909
B = 1000
XL = P("data/raw/GSE32472/GSE32472_series_matrix.txt/正常组与中重度BPD的比较/12-17d/炎症并集/Nomogram/预测值-系数.xlsx")

df = pd.read_excel(XL, sheet_name='glm模型')
coef = pd.read_excel(XL, sheet_name='多因素模型情况')
GENES = ['CD3G', 'PRTN3', 'CD6', 'THBS1']

y = df['outcome'].astype(int).values
X = df[GENES].values
lp = df['linear_predictors'].values
pp = df['predicted_probability'].values
n = len(y)
print(f'n = {n}   MS-BPD = {int(y.sum())}   controls = {int((1-y).sum())}')

# ---- 0. Internal consistency check ----
c = dict(zip(coef['多因素模型变量'], coef['系数']))
lp_recalc = c['(Intercept)'] + sum(c[g] * df[g].values for g in GENES)
print(f"\n[check] maximum deviation of the linear predictor from the coefficients = {np.max(np.abs(lp_recalc - lp)):.2e}")
pp_recalc = 1 / (1 + np.exp(-lp_recalc))
print(f"[check] maximum deviation of predicted_probability = {np.max(np.abs(pp_recalc - pp)):.2e}")

# ---- 1. Apparent AUC (using the supplied linear predictor directly) ----
auc_app = roc_auc_score(y, lp)
# Refit once to confirm agreement with the manuscript model
mdl = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
mdl.fit(X, y)
auc_refit = roc_auc_score(y, mdl.decision_function(X))
print(f"\n[apparent] AUC (given linear predictor) = {auc_app:.4f}")
print(f"[apparent] AUC (refitted four-gene logistic regression)  = {auc_refit:.4f}")
print(f"[coef] intercept={mdl.intercept_[0]:.6f}  " +
      "  ".join(f"{g}={v:.6f}" for g, v in zip(GENES, mdl.coef_[0])))

# ---- 2. DeLong 95% CI (CI of the apparent AUC) ----
def delong_ci(y, score, alpha=0.05):
    from scipy.stats import norm
    y = np.asarray(y); score = np.asarray(score)
    pos, neg = score[y == 1], score[y == 0]
    n1, n0 = len(pos), len(neg)
    auc = (np.subtract.outer(pos, neg) > 0).mean() + 0.5 * (np.subtract.outer(pos, neg) == 0).mean()
    v10 = np.array([(pos[i] > neg).mean() + 0.5 * (pos[i] == neg).mean() for i in range(n1)])
    v01 = np.array([(neg[j] < pos).mean() + 0.5 * (neg[j] == pos).mean() for j in range(n0)])
    s10, s01 = v10.var(ddof=1), v01.var(ddof=1)
    se = np.sqrt(s10 / n1 + s01 / n0)
    z = norm.ppf(1 - alpha / 2)
    return auc, max(0.0, auc - z * se), min(1.0, auc + z * se), se

a_dl, lo_dl, hi_dl, se_dl = delong_ci(y, lp)
print(f"[DeLong ] apparent AUC = {a_dl:.4f}  (95% CI {lo_dl:.4f}-{hi_dl:.4f}), SE = {se_dl:.4f}")

# ---- 3. 1000-replicate bootstrap internal validation (Harrell optimism correction) ----
rng = np.random.default_rng(RNG_SEED)
opt, auc_boot, auc_test, corrected_b = [], [], [], []
skipped = 0
for b in range(B):
    idx = rng.integers(0, n, n)
    Xb, yb = X[idx], y[idx]
    if len(np.unique(yb)) < 2:
        skipped += 1
        continue
    try:
        m = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
        m.fit(Xb, yb)
    except Exception:
        skipped += 1
        continue
    sb = m.decision_function(Xb)      # performance within the bootstrap sample
    so = m.decision_function(X)       # performance of the bootstrap model on the original sample
    ab = roc_auc_score(yb, sb)
    ao = roc_auc_score(y, so)
    auc_boot.append(ab); auc_test.append(ao)
    opt.append(ab - ao)               # optimism
    corrected_b.append(auc_app - (ab - ao))

opt = np.array(opt); corrected_b = np.array(corrected_b)
mean_opt = opt.mean()
auc_corr = auc_app - mean_opt
ci_corr = np.percentile(corrected_b, [2.5, 97.5])
ci_opt = np.percentile(opt, [2.5, 97.5])

print(f"\n[bootstrap B={B}, succeeded: {len(opt)} runs, skipped {skipped} runs]")
print(f"  mean optimism                 = {mean_opt:.4f}  (95% CI {ci_opt[0]:.4f}-{ci_opt[1]:.4f})")
print(f"  optimism-corrected AUC        = {auc_corr:.4f}")
print(f"  corrected AUC 95% CI (percentile)  = {ci_corr[0]:.4f}-{ci_corr[1]:.4f}")
print(f"  mean AUC within bootstrap         = {np.mean(auc_boot):.4f}")
print(f"  mean AUC of the bootstrap models on the original sample = {np.mean(auc_test):.4f}")

# ---- 4. 0.632+ estimate (supplementary; frequently requested by reviewers) ----
auc_632 = 0.368 * auc_app + 0.632 * np.mean(auc_test)
g = np.mean(auc_test)
R = abs(mean_opt) / (abs(auc_app - g) + 1e-12) if abs(auc_app - g) > 0 else 0
w = 0.632 / (1 - 0.368 * R)
auc_632p = (1 - w) * auc_app + w * g
print(f"\n[0.632 ]  = {auc_632:.4f}")
print(f"[0.632+] (R={R:.4f}, w={w:.4f}) = {auc_632p:.4f}")

# ---- 5. Sensitivity / specificity at the Youden cut-off ----
fpr, tpr, thr = roc_curve(y, lp)
j = tpr - fpr
k = int(np.argmax(j))
print(f"\n[Youden] threshold={thr[k]:.4f}  sensitivity={tpr[k]:.3f}  specificity={1-fpr[k]:.3f}")

# ---- 6. Save ----
out = {
    "n": int(n), "MS_BPD": int(y.sum()), "controls": int((1 - y).sum()),
    "apparent_AUC": round(float(auc_app), 4),
    "apparent_AUC_DeLong_95CI": [round(float(lo_dl), 4), round(float(hi_dl), 4)],
    "bootstrap_iterations": int(len(opt)),
    "random_seed": RNG_SEED,
    "mean_optimism": round(float(mean_opt), 4),
    "optimism_95CI": [round(float(ci_opt[0]), 4), round(float(ci_opt[1]), 4)],
    "optimism_corrected_AUC": round(float(auc_corr), 4),
    "corrected_AUC_95CI_percentile": [round(float(ci_corr[0]), 4), round(float(ci_corr[1]), 4)],
    "auc_0.632": round(float(auc_632), 4),
    "auc_0.632plus": round(float(auc_632p), 4),
    "youden_threshold": round(float(thr[k]), 4),
    "sensitivity": round(float(tpr[k]), 3),
    "specificity": round(float(1 - fpr[k]), 3),
}
with open('bootstrap_auc_result.json', 'w', encoding='utf-8') as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
pd.DataFrame({
    "bootstrap_id": np.arange(1, len(opt) + 1),
    "auc_in_bootstrap": np.round(auc_boot, 6),
    "auc_on_original": np.round(auc_test, 6),
    "optimism": np.round(opt, 6),
    "corrected_auc": np.round(corrected_b, 6),
}).to_csv('bootstrap_auc_iterations.csv', index=False, encoding='utf-8-sig')
print("\n[saved] bootstrap_auc_result.json / bootstrap_auc_iterations.csv")
