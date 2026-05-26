import numpy as np

mu = 0.6 # hr^-1
kappa = 5e-5 # hr ^-1
# doubling time
tau = np.log(2)/ mu
kappa_gen = kappa * tau
CN = 1 - np.log2(kappa_gen)

print(f"doubling time: {tau:.2f} hours")
print(f"genetic loss rate: {kappa_gen:.2e}")
print(f"PCN: {CN:.2f}") 
print(f"{1/1.1:.2f}")