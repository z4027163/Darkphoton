# Low-mass dark photon scouting search
## 2HDM Scalar model interpretation 

*******************************

### SOME COORDINATES for the calculation:

- Mass Points: the list of mass points used in the limit calculation is included in the notebook and it is reported here
```
1.176, 1.187, 1.199, 1.211, 1.223, 1.235, 1.248, 1.26, 1.273, 1.286, 1.299, 1.312, 1.325, 1.338, 1.351, 1.365, 1.378, 1.392, 1.406, 1.42, 1.434, 1.449, 1.463,
1.478, 1.493, 1.508, 1.523, 1.538, 1.553, 1.569, 1.584, 1.6, 1.616, 1.632, 1.649, 1.665, 1.682, 1.699, 1.716, 1.733, 1.75, 1.768, 1.785, 1.803, 1.821, 1.839,
1.858, 1.876, 1.895, 1.914, 1.933, 1.953, 1.972, 1.992, 2.012, 2.032, 2.052, 2.073, 2.094, 2.114, 2.136, 2.157, 2.179, 2.2, 2.222, 2.245, 2.267, 2.29, 2.313,
2.336, 2.359, 2.383, 2.406, 2.43, 2.455, 2.479, 2.504, 2.529, 2.554, 2.58, 4.201, 4.243, 4.286, 4.328, 4.372, 4.415, 4.46, 4.504, 4.549, 4.595, 4.641, 4.687,
4.734, 4.781, 4.829, 4.877, 4.926, 4.975, 5.025, 5.075, 5.126, 5.177, 5.229, 5.282, 5.334, 5.388, 5.442, 5.496, 5.551, 5.606, 5.663, 5.719, 5.776, 5.834, 5.892,
5.951, 6.011, 6.071, 6.132, 6.193, 6.255, 6.318, 6.381, 6.445, 6.509, 6.574, 6.64, 6.706, 6.773, 6.841, 6.909, 6.978, 7.048, 7.119, 7.19, 7.262, 7.334, 7.408,
7.482, 7.557, 7.632, 7.709, 7.786, 7.864
```

- Acceptance:
A polynomial fit is obtained starting from the acceptances for `mass = [1., 2., 4., 5., 6., 7., 8., 20.] GeV` (retrieved using the Pythia generator).
The fit is then used in the notebook to evaluate the acceptance for each mass of interest.

- Branching ratio:
The list of branching ratios for each mass point under consideration is taken from `BR_IV_tgBeta05_allLimits.txt`.

- Cross section:
The list of cross sections for each mass point under is obtained by using a TGraph fit (see `allMassPoints_fit.C`) and is available in `xsec_muR05_muF05_allLimits.txt`.

- Obs limit:
The list of observed limits for each mass point is contained in `CMS_modelIndependentLimits_red_bothYears_twoIDs_90CL_allLimits_20220401.txt`.
(The previous set of limits is also provided for easy comparison and as a sanity check of the updates of the notebook: `CMS_modelIndependentLimits_red_bothYears_twoIDs_90CL_20210723_doubleCheck.txt`)


### USEFUL REFERENCES:

- Collider constraints on light pseudoscalars - JHEP 03 (2018) 178:
https://link.springer.com/article/10.1007/JHEP03(2018)178

- LHCB PAPER - JHEP 10 (2020) 156:
https://link.springer.com/article/10.1007%2FJHEP10%282020%29156

- CMS EXO-19-018 - Phys. Rev. Lett. 124, 131802:
https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.124.131802

- CMS at 7 TeV:
https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.109.121801
