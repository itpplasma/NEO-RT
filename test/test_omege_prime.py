"""
Compares gradients of canonical frequencies with the ones computed by linear regression.
"""

#%%
import numpy as np
from numpy.testing import assert_allclose
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

df = pd.read_csv("/tmp/test_omega_prime.dat", delim_whitespace=True, comment='#')

# %%
model = LinearRegression()
fit = model.fit(df[["s", "v", "eta"]], df["Omth"])
dOmth_num = fit.coef_
print("dOmth/d(s, v, eta) from numerical calculation:")
print(dOmth_num)
print("dOmth/d(s, v, eta) from code:")
print(df[["dOmthds", "dOmthdv", "dOmthdeta"]].iloc[0].values)

#assert_allclose(dOmth_num[0], df["dOmthds"].iloc[0], rtol=1e-6)
#print ("assert_allclose dOmth/ds OK")
assert_allclose(dOmth_num[1], df["dOmthdv"].iloc[0], rtol=1e-6)
print ("assert_allclose dOmth/dv OK")
assert_allclose(dOmth_num[2], df["dOmthdeta"].iloc[0], rtol=1e-6)
print ("assert_allclose dOmth/deta OK")

# %%
model = LinearRegression()
fit = model.fit(df[["s", "v", "eta"]], df["Om"])
dOm = fit.coef_
print("dOm/d(s, v, eta) from numerical calculation:")
print(dOm)
print("dOm/d(s, v, eta) from code:")
print(df[["dOmds", "dOmdv", "dOmdeta"]].iloc[0].values)

#assert_allclose(dOmth_num[0], df["dOmthds"].iloc[0], rtol=1e-6)
#print ("assert_allclose dOmth/ds OK")
assert_allclose(dOm[1], df["dOmdv"].iloc[0], rtol=1e-6)
print ("assert_allclose dOm/dv OK")
assert_allclose(dOm[2], df["dOmdeta"].iloc[0], rtol=1e-6)
print ("assert_allclose dOm/deta OK")

# %%
print(f"Omega prime old: {np.mean(df['Ompr_old'][0])} +/- {np.std(df['Ompr_old'])}")
print(f"Omega prime new: {np.mean(df['Ompr_new'][0])} +/- {np.std(df['Ompr_new'])}")

# %%
plt.figure()
plt.plot(Jbar3, Om, "x")

# %% Compare partial derivatives of frequencies

#vars = np.empty((len(Om), 3))
# vars[:, 0] = v
# vars[:, 1] = eta
# vars[:, 2] = J3  # pphi
vars = freq_data[:, 0:3]

fit = LinearRegression().fit(vars, Om)
print(f"Omega derivatives from numerical calculation:\n {fit.coef_}")

# %%
