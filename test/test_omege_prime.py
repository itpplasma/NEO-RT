# %%
"""
Compares gradients of canonical frequencies with the ones computed by linear regression.
"""

import numpy as np
from numpy.testing import assert_allclose
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

df = pd.read_csv("/tmp/test_omega_prime.dat", delim_whitespace=True, comment='#')

# %% Omth derivatives
model = LinearRegression()
fit = model.fit(df[["s", "v", "eta"]], df["Omth"])
dOmth_num = fit.coef_
print("dOmth/d(s, v, eta) from regression:")
print(dOmth_num)
print("dOmth/d(s, v, eta) from d_Om_ds:")
print(df[["dOmthds", "dOmthdv", "dOmthdeta"]].iloc[0].values)

assert_allclose(dOmth_num[0], df["dOmthds"].iloc[0], rtol=1e-2)
print("assert_allclose dOmth/ds OK")
assert_allclose(dOmth_num[1], df["dOmthdv"].iloc[0], rtol=1e-3)
print("assert_allclose dOmth/dv OK")
assert_allclose(dOmth_num[2], df["dOmthdeta"].iloc[0], rtol=1e-3)
print("assert_allclose dOmth/deta OK")

# %% Om_tE derivatives
model = LinearRegression()
fit = model.fit(df[["s", "v", "eta"]], df["Om_tE"])
dOmtE_num = fit.coef_
print("Om_tE/d(s, v, eta) from regression:")
print(dOmtE_num)
print("Om_tE/ds from Om_tE_ds:")
print(df["dOm_tEds"][0])

assert_allclose(dOmtE_num[0], df["dOm_tEds"][0], rtol=1e-3)
print("assert_allclose dOm_tE/ds OK")

# %% Omph derivatives
model = LinearRegression()
fit = model.fit(df[["s", "v", "eta"]], df["Omph"])
dOmph_num = fit.coef_
print("dOmph/d(s, v, eta) from regression:")
print(dOmph_num)
print("dOmph/d(s, v, eta) from d_Om_ds:")
print(df[["dOmphds", "dOmphdv", "dOmphdeta"]].iloc[0].values)

assert_allclose(dOmph_num[0], df["dOmphds"].iloc[0], rtol=1e-2)
print("assert_allclose dOmph/ds OK")
assert_allclose(dOmph_num[1], df["dOmphdv"].iloc[0], rtol=1e-2)
print("assert_allclose dOmph/dv OK")
assert_allclose(dOmph_num[2], df["dOmphdeta"].iloc[0], rtol=1e-2)
print("assert_allclose dOmph/deta OK")

# %% Om derivatives
model = LinearRegression()
fit = model.fit(df[["s", "v", "eta"]], df["Om"])
dOm = fit.coef_
print("dOm/d(s, v, eta) from regression:")
print(dOm)
print("dOm/d(s, v, eta) from d_Om_ds:")
print(df[["dOmds", "dOmdv", "dOmdeta"]].iloc[0].values)

assert_allclose(dOmth_num[0], df["dOmthds"].iloc[0], rtol=1e-2)
print ("assert_allclose dOmth/ds OK")
assert_allclose(dOm[1], df["dOmdv"][0], rtol=1e-3)
print ("assert_allclose dOm/dv OK")
assert_allclose(dOm[2], df["dOmdeta"][0], rtol=1e-3)
print ("assert_allclose dOm/deta OK")

# %%
model = LinearRegression()
fit = model.fit(df[["Jbar1", "Jbar2", "Jbar3"]], df["Om"])
dOm = fit.coef_
print("dOm/d(Jbar1, Jbar2, Jbar3) from regression:")
print(dOm)

print(f"Omega prime old: {np.mean(df['Ompr_old'][0])} +/- {np.std(df['Ompr_old'])}")
print(f"Omega prime new: {np.mean(df['Ompr_new'][0])} +/- {np.std(df['Ompr_new'])}")

# %%
plt.figure()
plt.plot(df["Jbar2"], df["Om"], "x")


# %%
