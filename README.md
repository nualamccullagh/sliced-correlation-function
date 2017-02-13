# sliced-correlation-function
The 2-point correlation function of matter measures how clustered matter is as a function of scale. In a Gaussian Random Field, the clustering in regions below the mean density of the field (under-dense regions) is the same as the clustering in regions above the mean density (over-dense regions). However, in a non-Gaussian field, such as matter in the universe, this is not necessarily the case. This density-dependent clustering, also known as the sliced correlation function, can provide us more insight into the properties of the universe than the usual 2-point correlation function alone. 

This Python code can be used to model and plot the sliced correlation function (see https://arxiv.org/abs/1610.06215 for details). Two different models are considered for the matter density field: Gaussian and lognormal. An option for modelling galaxy bias is also included.

Files:
-sliced_corr_bias.py: Main code to compute the sliced correlation function for a given input power spectrum and given bias parameters.
-sliced_corr_plots.py: Shows how to run the code for both the Gaussian and lognormal models.

