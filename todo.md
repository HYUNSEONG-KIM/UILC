
---------------------------------------------------------------

2023-02-03
1. Implement Rectangular ESC n_max array search routine
2. OP: bc_expansion
3. OP: linear solve -> power weight name change
4. OP: binarization routine - (1: KDE), (2: FM modulation)
5. Uniformity factor estimation routinesimplementation

---------------------------------------------------------------

How to find exact distribution of LEDs from solution of linear equaton.

SOl: Non-negative solution, but it is not normalized
    Permit power allocation: 
        Just distribute sources to the peak points and allocate powers to its peak values.
    Prohibit power allocation(Only location distribution):
        Approximating given different powers to dense of sources.


Method 1:
    KDE(Kernel distribution estimating)

    Using kernel distribution: Gaussian
        or Wavelet(It can be negative so that expand the square root of probability distribution)
    
Method 2:
    Apply frequency modulation
        Frequency: superposition of gaussians(same shape, rho)

        ESC: Method result( Center uniform array ) => Base frequency
        Using solution of linear equation:- > Frequency modulate the ESC array to boundary,

$$x(t) = A_c \cos(2 \pi f(t)t + \phi_c)\\ f(t) = f_b + sol(t)$$ 