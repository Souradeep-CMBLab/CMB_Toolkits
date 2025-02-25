# Using CoNIGS

1. Enable unlimited stack size

```bash
ulimit -s unlimited
```

1. Put the spectra file in spectra folder
2. Run following codes

```bash
./cholesky_new_unbin_L_1_2048_planck_new.exe
./gen_alm_real_L_1_final_new_3.exe
./gen_alm_img_L_1_final_new_3.exe
```

1. Cnvert the output in the format of healpy

```python
import numpy as np
import healpy as hp
#import matplotlib.pyplot as plt

almr1 = np.loadtxt("alm_real_beta10_lmax2000_nside2048_1.dat")
almi1 = np.loadtxt("alm_img_beta10_lmax2000_nside2048_1.dat")
almr2 = np.loadtxt("alm_real_beta10_lmax2000_nside2048_2.dat")
almi2 = np.loadtxt("alm_img_beta10_lmax2000_nside2048_2.dat")

N_side = 2048
elmax = 2000
mmax = 2000
alm_ar = np.ndarray(shape = (elmax+1, elmax+1), dtype = np.complex)

index = 0
for i in np.arange(0, elmax+1, 1):
    for j in np.arange(i, -1, -1):
        alm_ar[i, j] = np.complex(almr1[index], almi1[index])
        index += 1

alm1 = np.ndarray(shape = (len(almr1)), dtype = np.complex)

index = 0
for j in np.arange(0, elmax+1, 1):
    for i in np.arange(j, elmax + 1, 1):
        alm1[index] = alm_ar[i, j]
        index += 1

map1 = hp.alm2map(alm1, lmax = elmax, mmax=elmax, nside = N_side)

N_pix = 12*(N_side**2)
mu, sigma = 0.0, 5.0 # mean and standard deviation
np.random.seed(seed=10)
Noise = np.random.normal(mu, sigma, N_pix)
map_2048_withnoise_1st = map1 + Noise

np.savetxt("map_beta10_nside2048_WithNoise_1st.dat", map_2048_withnoise_1st)

alm_ar = np.ndarray(shape = (elmax+1, elmax+1), dtype = np.complex)

index = 0
for i in np.arange(0, elmax+1, 1):
    for j in np.arange(i, -1, -1):
        alm_ar[i, j] = np.complex(almr2[index], almi2[index])
        index += 1

alm1 = np.ndarray(shape = (len(almr2)), dtype = np.complex)

index = 0
for j in np.arange(0, elmax+1, 1):
    for i in np.arange(j, elmax + 1, 1):
        alm1[index] = alm_ar[i, j]
        index += 1

map2 = hp.alm2map(alm1, lmax = elmax, mmax=elmax, nside = N_side)

#N_pix = 12*(N_side**2)
mu, sigma = 0.0, 5.0 # mean and standard deviation
np.random.seed(seed=18)
Noise = np.random.normal(mu, sigma, N_pix)
map_2048_withnoise_2nd = map2 + Noise

np.savetxt("map_beta10_nside2048_WithNoise_2nd.dat", map_2048_withnoise_2nd)
```
