This zip file contains a Matlab implementation of the opponent-channel
model described in:

Briley PM & Summerfield AQ. Age-related deterioration of the
representation of space in human auditory cortex.  Neurobiology of
Aging.

This function can be used to reproduce the key aspects of figures 4
and 5 of our paper. It includes the fitted parameter values for our
three participant groups (young, younger-old and older-old).

Specific details of model parameters and outputs are given at the
start of the function, and the function is annotated throughout. The
model is derived from that of Briley et al. (2013),
J. Assoc. Res. Otolaryngol. 14, 83-101, which is also available on
ModelDB (146050). At the core of the model is the description of the
spatial tuning of two channels, one maximally responsive to the left,
and one to the right, auditory hemispace. Our model can be used to
describe the tuning of these two channels, and consequently to predict
psychophysical spatial acuity ("minimum audible angles", MAAs) for
different azimuthal positions, as well as electroencephalographic
(EEG) responses to abrupt shifts in sound-source location.

The function has been re-written, but one key change from the previous
version is that the tuning of each channel is now described by the pdf
of the generalized Gaussian distribution, also known as the
exponential power distribution, which has both a width parameter
(related to the standard deviation of the conventional Gaussian
distribution) and a shape parameter, which allows the function to
adopt the conventional Gaussian shape, or to have a sharper or flatter
peak.


out = oppchanmodel(in,figs);

in: set to '' (default parameters), 'young' or 'youngold' or 'oldold'
(parameters fitted to the EEG data from one of our three participant
groups). Alternatively, in can be a structure with fields chans, ac
and pred. These fields specify the model parameters (see defaults for
more information).

figs: set to 0 (no figures), 1 (channel figures, including MAA
predictions), 2 (auditory-cortex-response figures), 3 (all figures).

out: structure containing the updated chans (channel parameters,
tuning curves and gradients of tuning curves), ac (auditory cortex
parameters and EEG response curves) and pred (MAA predictions) fields.
