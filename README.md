# Stability-of-Quasi-Biennial-Oscillations
The Holton--Lindzen--Plumb is a model that mdoels the Quasi-Biennial-Oscillation, the derivation of the model utilizes quasi-linear and WKB approximation.
The Code creates simulations of the model, the first code uses the symmetric version of the model, and the stability of the model uses eigenvalues.
The code for asymmetric model as been labelled with an ending "_a", the simulations and the stability is virtually similar to the symmetric version.

A quick outline of the code:
The model is controlled using Reynolds Number, so you can change this factor to get different results.

The code uses ODE45 to do the simulations.

It calculates the period of oscillations and give you a Hovmöller diagram to visualize the oscillations.
