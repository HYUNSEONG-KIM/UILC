<table style="border-radius:8px;width:100%;">
<th style="text-align:center;background-color:rgb(0, 0, 0); color:white; border-top-left-radius: 10px;width:20%;">
Thm</th>
<th style="text-align:left; width:80%;">
Cauchy interlace </th>
<tr style="text-align:center;">
<td colspan="2">


For Hermit matrix $\mathbf{A}$ of order $n$ and its principle sub-matrix $\mathbf{B}$ of order $n-1$, about their ordered eigen values $\{ 
\lambda _i \}_{i=1}^n$ and $\{\mu_k \}_{k=1}^{n-1}$ respectively,

next inequality is hold.

$$\lambda_1 \leq \mu_1 \leq \lambda_2 \leq \mu_2 \leq \cdots \leq \mu_{n-2} \leq \lambda_{n-1} \leq \mu_{n-1} \leq \lambda_{n}$$


</td>
</tr>

</table>


Let $T_n$ be a Toeplitz matrix of sequence $\{a_k\}_{k=-(n-1)}^{n-1}$ of non-sigular matrix.

$$T := \begin{pmatrix}
a_0 & a_1 & a_2 & \cdots & a_{n-1}\\
a_1 & a_0 & a_1 & \cdots & a_{n-2}\\
a_2 & a_1 & a_0 & \cdots & a_{n-3}\\
{}  & \vdots & {} & \ddots & \vdots\\
a_{n-1} & a_{n-2} & a_{n-3} & \cdots & a_{0}\\
\end{pmatrix}$$

Thus, $\det(T) = \Pi_{i=1}^n \mu_i \neq 0$.



Now, the arise problem of raidation pattern convolution matrix is that discretized matrix is non-singular by the dimension of matrix or not. By the definition of radiation pattern function on $xy$ plane. Such function $R(x,y)$ always positive for all function, thus generated matrix always positive Hermit Toeplitz matix. There are various radiation parttern by the manufacturers and models of them, in spite of the such variance, those patterns can be approximated with linear combination of angular shifted cosine power(imperfect Lambertian) or Gaussian. Thus, showing these two functions properties are enough to use the method in practical cases.

$$R_g(a_1; a_2; a_3; \theta) = a_1 \exp(- \ln_2 (\frac{|\theta| - a_1}{a_2})^2)\\
R_l(s; \theta) = a_1 \cos^s(\theta)\\
\\
R_g(a_1; a_2; a_3; r) = a_1 \exp(- \ln_2 (\frac{|\arctan(\frac{r}{H})| - a_1}{a_2})^2)\\
R_l(s; r) = \frac{H^s}{(H^2 + (r^2))^{s/2 +1}}\\
R_{g2}(r) = \frac{s}{2 \pi H^2} \exp(- \frac{s}{2} (\frac{d}{H})^2)
$$



>Ivan Moreno and Ching-Cherng Sun, "Modeling the radiation pattern of LEDs," Opt. Express 16, 1808-1819 (2008) 

> Hongming Yang, JanW. M. Bergmans, Tim C. W. Schenk, Jean-Paul M. G. Linnartz, and Ronald Rietman, "An analytical model for the illuminance distribution of a power LED," Opt. Express 16, 21641-21646 (2008) 

> Hongming Yang, Bergmans, J. W. M., Schenk, T. C. W., Linnartz, J.-P. M. G., & Rietman, R. (2009). Uniform Illumination Rendering Using an Array of LEDs: A Signal Processing Perspective. IEEE Transactions on Signal Processing, 57(3), 1044–1057. doi:10.1109/tsp.2008.2009269 

Convinence of calculation, assumed that $R(r, \theta, \phi)$ has a axial symmetry then, on $xy$ plane the radiation function can be reduced to 1-dim function, $T(t) = R(|t|)$. 

$$t_{k} = T(\Delta_n \cdot k)$$

$$\mathbf{T}_n = \begin{bmatrix} 
T(0)& T(\Delta_n)& T(2\Delta_n)& \cdots & T((n-1)\Delta_n)\\
T(\Delta_n)& T(0)& T(\Delta_n)& \cdots & T((n-2)\Delta_n)\\
& \vdots& & \ddots & \vdots\\
T((n-1)\Delta_n)& T((n-2)\Delta_n)&T((n-3)\Delta_n) & \cdots & T(0)\\
\end{bmatrix}$$


Since, it is a Hermit Toeplitz matrix, we can applying Szeg ̈os theorem, See next reference.

> Gray, Robert M. "Toeplitz and circulant matrices: A review." Foundations and Trends® in Communications and Information Theory 2.3 (2006): 155-239.

Thus, for generating function, $f_g(\lambda)$ defined as,

$$f_g(\lambda) = \sum_{k =-\infty}^\infty t_k e^{ik\lambda}; \lambda \in [0 , 2 \pi] \\
= \sum_{k =-\infty}^\infty T(k \cdot \Delta_n) e^{ik\lambda}$$

It is a discrete Fourier transform of radiation function.


The Fourier Transform of radiation functions are

$$\mathcal{F}[R_g] = \frac{1}{H}\sqrt{\frac{s}{2 \pi}} \exp(- \frac{2 H^2 \pi^2}{s} t^2)\\
\mathcal{F}[R_l] = \frac{\sqrt{H^{s-1}}}{\sqrt{2}^s\Gamma[1+\frac{s}{2}]} |t|^{s/2+1} K_{\frac{s+1}{2}}(H |t|)$$

and those are positive function for all $t\in \mathbb{R}$.

Thus any sampled Toeplitz matrix $\mathbf{T}_n$ its minumum value of eigen values are greater or equal to zero.

However, with positive definite function Theorem. we can get more strict result. See next reference.

> Gregory E Fasshauer, Meshfree Approximation Methods With Matlab, 6th SEries of Interdisciplinary MAthematical Sciences, World Scientific Publishing Company, 2007, ISBN: 9813101571.

**Definition**: Positive definite function

A complex-valued continuouse function $\Phi: \mathbb{R}^s \rightarrow \mathbb{C}$ is called *positive definite* on $\mathbb{R}^s$ if 

$$\sum_{j=1}^{N} \sum_{k=1}^{N} c_j \overline{c_k} \Phi(x_j - x_k) \geq 0$$

for any $N$ pairwise different points $x_1, \dots , x_N \in \mathbb{R}^s$ and $c = [c1, \dots, c_N]^T \in \mathbb{C}^N$. If the function $\Phi$ is zero only for $\mathbb{c} =0$ on the above equation, then called *strictly positive definite* on $\mathbb{R}^s$.

**Properties**: $\Phi_i$ are positive definite functions.

1. Non-negative finite linear combinations(coefficients are non-negative) of positive definite functions are positive definite.
2. From non-negative finite linear combinations of positive definite functions, if at least one of function is strictly positive definite and corresponding coefficient is non-zero then combination also be a strictly positive definite.
3. $\Phi_i(0) \geq 0$
4. If $\Phi_i(0) =0$ then, $\Phi_i(x) = 0, \forall x$
5. $\Phi_i(-x) = \overline{\Phi_i(x)}$
6. Any $\Phi_i$ is bounded, preciously, $|\Phi(x)|\leq\Phi_i(0)$
7. Products of (strictly) positive definite functions is (strictly) positive definite function.


Thus, $\min(\{\lambda_{T_n}\})$ is monotonic decreasing function by interger $n$ increasing with lower bound $0$. 

With Cauchy Interace theorem, for $n$ dim Toeplitz matrix $T_n$ and its principle submatrix of dimension $n-1$, $T_{n-1}$ their ordered eigenvalues $\{\lambda_i \}, \{\mu_i \}$. By the above result, the lower bound is $0$ so that



$$\text{Tr}(T_n) = na_0 = \sum_{i=1}^n \lambda_i\\
\text{Tr}(T_{n-1}) = (n-1) a_0 = \sum_{i=1}^{n-1} \mu_i
$$

$$a_0 = \sum_{i=1}^n \lambda_i - \sum_{i=1}^{n-1} \mu_i\\
= \lambda_1  + \sum_{i=1}^{n-1}(\lambda_{i+1} - \mu_i)$$

Schur Complement: 

$$\det(T_n) = \det(T_{n-1}) (a_0 - l_{n-1}^t T_{n-1}^{-1}l_{n-1})$$

$$T_{n-1} y = l_{n-1}$$

<table style="border-radius:8px;width:100%;">
<th style="text-align:center;background-color:rgb(0, 0, 0); color:white; border-top-left-radius: 10px;width:20%;">
Thm</th>
<th style="text-align:left; width:80%;">
Cauchy interlace </th>
<tr style="text-align:center;">
<td colspan="2">


Let $T_n$ be a Toeplitz matrix of sequence $\{a_k\}_{k=-(n-1)}^{n-1}$, if $T_{n-1}$ is non-singular and all egien value 


</td>
</tr>

</table>

Now $T_{n+1}$

