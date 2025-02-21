#%%
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

freq_data = np.loadtxt('/tmp/test_omega_prime.dat')

s = freq_data[:, 0]
v = freq_data[:, 1]
eta = freq_data[:, 2]
J1 = freq_data[:, 3]
J2 = freq_data[:, 4]
J3 = freq_data[:, 5]
Jbar1 = freq_data[:, 6]
Jbar2 = freq_data[:, 7]
Jbar3 = freq_data[:, 8]
Om = freq_data[:, 9]
Ompr_old = freq_data[:, 10]
Ompr_new = freq_data[:, 11]

fit = LinearRegression().fit(freq_data[:, 6:9], Om)

Ompr_num = fit.coef_[2]

print(f"Omega prime from numerical calculation: {Ompr_num}")
print(f"Omega prime old: {np.mean(Ompr_old[0])} +/- {np.std(Ompr_old)}")
print(f"Omega prime new: {np.mean(Ompr_new[0])} +/- {np.std(Ompr_new)}")

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
