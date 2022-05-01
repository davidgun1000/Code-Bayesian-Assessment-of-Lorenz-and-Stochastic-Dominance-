# Code-Bayesian-Assessment-of-Lorenz-and-Stochastic-Dominance-
Code for Bayesian Assessment of Lorenz and stochastic dominance

This is the code for the paper "Bayesian assessement of Lorenz and stochastic dominance"
David Gunawan, William Griffiths, and Duangkamon Chotikapanich (2020) Canadian Journal of Economics, 53 (2), 767-799.
The main code is dominance_prog_gamma.m.
The main code will generate the posterior probabilities of Lorenz and stochastic dominance. 
Please read the paper for more information.
The user needs to supply their own data.

data_dominance_paper.mat is an example of the required data for this package. 
First, for each year, you need to use gamma_mixture4.m to give the posterior draws of the parameters of mixture of gamma. The files gamma4_1999.mat, gamma4_2002.mat,
gamma4_2005.mat, and gamma4_2008.mat are the examples.
Second, using the posterior draws obtained from the first step, you then use dominance_prog_gamma.m to obtain the posterior 
probabilities of Lorenz and stochastic dominance.
