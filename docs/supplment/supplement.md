



# Binarization of Solution

## Kernel Density Estimation

### Gaussian Kernel

Gaussian Kernel is form distribution itself and consists dirac sequence.

## Amplitude/Frequency Modulation

This method treat the original solution as 
amplitude modulated signal for spatial dimension.


1. Generate message signal, $m(x)$, from solution
2. Calculate Frequency modulated signal which represents the location of sources of which its maximum peak points.

$$o(x)= (m(x)+1) \cos(2 \pi \nu_b x + \phi_b)$$
$$f(x) = \cos(2 \pi (\nu_b + k\cdot m(x)) x + \phi_f)$$

About base frequency, $f_b$

* DS: center peak period, $T = \frac{1}{f_b}$
* NNLS: Same with DS but for large number of discretization it becomes the average values of peak terms. 

NNLS peak calculation: Wavelet based peak estimation
