import numpy as np, pandas as pd
from sklearn.linear_model import lasso_path
from grnbench.infer.classical import lasso, pearson_signed

rng = np.random.default_rng(0)
n = 300
true_tf = rng.normal(0, 1, n)                 # drives the target, mean 0-ish after shift
high_mean_tf = 8.0 + rng.normal(0, 0.05, n)    # independent of target, high mean (log-normalized scale)
target = 4.0 + 0.8 * true_tf + rng.normal(0, 0.3, n)
X = np.column_stack([true_tf + 2.0, high_mean_tf, target])   # all positive like log1p data
expr = pd.DataFrame(X, columns=["TRUE_TF", "HIGHMEAN_TF", "TARGET"])
print("corr(TRUE_TF,TARGET)=%.3f corr(HIGHMEAN_TF,TARGET)=%.3f" % (
    np.corrcoef(X[:,0], X[:,2])[0,1], np.corrcoef(X[:,1], X[:,2])[0,1]))
e = lasso(expr, tf_names=["TRUE_TF", "HIGHMEAN_TF"])
print(e[e.target == "TARGET"].to_string())
# same fit with centering (what an intercept does)
Xr = X[:, :2]; y = X[:, 2]
_, coefs, _ = lasso_path(Xr - Xr.mean(0), y - y.mean(), n_alphas=20)
print("centered coefs at knee:", np.abs(coefs[:, coefs.shape[1] // 2]))
_, coefs, _ = lasso_path(Xr, y, n_alphas=20)
print("uncentered (as shipped) coefs at knee:", np.abs(coefs[:, coefs.shape[1] // 2]))
