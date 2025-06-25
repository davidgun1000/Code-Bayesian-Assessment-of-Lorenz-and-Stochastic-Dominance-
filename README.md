# Code-Bayesian-Assessment-of-Lorenz-and-Stochastic-Dominance-
This repository contains the code associated with the paper:
"Bayesian Assessment of Lorenz and Stochastic Dominance"
by David Gunawan, William Griffiths, and Duangkamon Chotikapanich,
Canadian Journal of Economics, 2020, Vol. 53(2), pp. 767–799.

The primary script is dominance_prog_gamma.m, which computes the posterior probabilities of Lorenz and stochastic dominance.
For further details on the methodology and interpretation, please refer to the published paper.

### 1. Prepare Your Data

- Provide your own dataset in the required format.

### 2. Generate Posterior Draws of Gamma Mixture Parameters

- Use 'gamma_mixture4.m' to estimate the parameters of a 4-component gamma mixture model for each year or subgroup.
- Example output files: 'gamma4_1999.mat', 'gamma4_2002.mat', etc.

### 3. Compute Posterior Probabilities of Dominance

- Use the posterior draws from Step 2 as input into 'dominance_prog_gamma.m'.
- The script will compute the posterior probabilities of:
  - Lorenz dominance
  - First-order stochastic dominance
  - Second-order stochastic dominance
