# LED design for uniform illumination area of rectangular shape

## Basic terminologies and functions.

Lambertian model:

$$I_s(r, \theta) = \frac{I_0}{r^2}\cos^s(\theta)$$

System parameter(Lambertian)

* $s$: Optical model parameter: Lambertian coefficient.
* $H$: Distance between measure plane and source plane.
* $W$ : Width of permitted area(Linear). 
* $W_1, W_2$: Widths of permitted area(Rectangular).

## Criterion Methods

### ESC method

<table style="border-radius:8px;width:100%;">
<th style="text-align:center;background-color:rgb(0, 0, 0); color:white; border-top-left-radius: 10px;width:20%;">
Def</th>
<th style="text-align:left;">
Sparrow Criterion </th>
<tr style="text-align:center;">
<td colspan="2">
For two LEDs and superposition of those intensity function,

$$I(t) = I(-d/2,t) + I(d/2,t)$$

The resolution limit $d$ that combined intensity curve remains single source is yielded from $\frac{\partial^2 I}{\partial t^2}|_{t=0} = 0$.


</td>
</tr>
<tr><td> Def </td> <td > <b>Expended Sparrow Criterion</b></td></tr>
<tr style="text-align:center;">
    <td colspan="2">

For $N$ number of LEDs which are uniformly distributed with distance $d$ each other, the flat(or uniform) condition of intensity is 

$$\frac{\partial^2 I_N}{\partial t^2}|_{t=0} = 0$$

where, 

$$I_N (t) = \sum_{i=1}^N I( (- \frac{N}{2} + (i-1)) d, t)$$


</td>
</tr>
</table>

LEDs are uniformly distributed, their uniform distance $r$ is calculated from coefficient $d$ and distance between target plane and LED plane $h$.

$$r = d \cdot h$$

Finding root of next equations for linear and rectangular array yield flat coefficient $d$.

linear

$$f(d) = \sum_{n=1}^N \frac{[1-(m+3)(N+1-2n)^2 (\frac{d}{2})^2]}{[(N+1-2n)^2 (\frac{d}{2})^2]^{(m+6)/2}}$$

rectangular

$$f(d) = \sum_{i=1}^N \sum_{j=1}^M  \frac{\{1-[(m+3)(N+1-2i)^2 -(M+1-2j)^2](\frac{d}{2})^2\})}{(\{ [N+1-2i)^2 + (M+1-2j)^2](\frac{d}{2})^2+1\}^{(m+6)/2}}$$

Special cases

1. 2 source linear array: $d = \sqrt{\frac{4}{s+3}}$.
2. 4 source square array: $d =\sqrt{\frac{4}{s+2}} $
3. Approximation, $s > 30$, $N \times M > 4 \times 4$ square: $d \approx \sqrt{\frac{1.2125}{s-3.349}}$
4. Approximation, $s> 30$, $N>4$: $d \approx  \sqrt{\frac{3.2773}{s+4.2539}}$

References:
> R. Barakat, “Application of apodization to increase two-point resolution by the sparrow criterion. i. coherent248 illumination,” J. Opt. Soc. Am. 52, 276–283 (1962).

> Moreno, I., Avendaño-Alejo, M., & Tzonchev, R. I. (2006). Designing light-emitting diode arrays for uniform near-field irradiance. Applied optics, 45(10), 2265-2272.


---

### Boundary Center Matching

Boundary Center Matching expansion method.

$$I(\theta, r) = I(H, s;x,t) = I_0 \frac{H^m}{(H^2 + (x-t)^2)^(\frac{s}{2}+1)}$$

The difference function $Di(x)$ of two LEDs, located in center symmetry, intensity contribution between two points, center and boundary is defined as

$$Di(x) : = (I(x, \frac{W}{2} + I(-x, \frac{W}{2})) - I(x,0) + I(-x,0))$$

It can be nomailized to $D(\alpha, d)$ with next equation.

$$D(\frac{W}{H}; \frac{x}{H}) = \frac{H^2}{I_0}Di(x)$$

$x_e$ is a point that makes boundary and center contribution of LEDs same.

$x_m$ is a point for center corresponding to $x=W/2$ for boundary. 

$$x_e = \text{root}(D(d))*H$$
$$x_m = \text{root}(D(d)+D(\alpha/2))*H$$

$\mathbf{P} := [0,x_m), \mathbf{Q} := [x_m,x_e], \mathbf{R} := (x_e, W/2]$

For $x_q \in \mathbf{Q}$ and $x_r \in \mathbf{R}$, there exist $x'_r \in \mathbf{R}$ and $x'_q \in \mathbf{Q}$ such that,

$$I_c(x_q) + I_c(x'_r) = I_b(x_q) + I_c(x'_r)$$
$$I_c(x_r) + I_c(x'_q) = I_b(x_r) + I_c(x'_q)$$

Then, for Esc array which consist of LEDs in $\mathbf{Q}$ region in given area $W$, the LEDs of which the center and boundary intensity contributions are same can be found.

$$\{\overbrace{x_1, x_2, ... , x_{m-1} }^\mathbf{Q}, \underbrace{x_m ,...,x_{n-1} x_n}_\mathbf{R}\}$$


$$\sum_{i=1}^n I_c (x_i) = \sum_{i=1}^n I_b (x_i) $$

Now considering $\mathbf{P}$ region sources, 

$$\{x_j\}_{j=1}^{n_p}, x_j \in \mathbf{P}$$

KDE + $D$ function interpolation can be adopted to get $\mathbf{R}$ sources corresponding sources for $\mathbf{P}$.


1. $x_{q, i}, x_{r, i} \rightarrow Di(x_{q, i}) + Di(x_{r, i}) = 0$
2. For $k$ such that $x_k \in \mathbf{Q}, x_{k+1} \in \mathbf{R}$, $|x_{k+1} -x_k| \approx d = avg(d_q)$.
3. Search corresponding points on $R$ region: 
    1. Calculate radiation, $I_{pqr}$ of $p, q, r_q$ sources.
    2. $\Epsilon: = (I_{pqr})_{max} - I_{pqr}$.
    3. Continuouse gauss elimination on $[x_e, W/2]$ region.


Practical approximation are next equations.

$$x_m = (\sqrt{2^{\frac{2}{s+2}}-1}) \cdot H$$
$$x_e = (\frac{1}{6}d_m + \frac{1}{4}\alpha) \cdot H$$

These are assumed as initial guessing points in finding algorithm. 

---

## Integral Equation

With radiation pattern, $R$, on $xy$ measure plane with source point $(t, l)$ on $tl$ plane with distance $h$.

$$R_h(x, y, t, l)$$ 

For example, Lambertian radiation with inverse law is,

$$R_h(x, y, t, l) := \frac{h^s}{(h^2 + (x-t)^2 + (y-l)^2)^(s/2+1)}$$

The illumination values on the measure plane point, $(x_i, y_i)$ can be represented as integral equation with source distribution, $\rho(t, l)$.

$$I(x_i, y_i) = \int_{-W_y/2}^{W_y/2} \int_{-W_x/2}^{W_x/2} R_h(x_i, t, y_i, l) \rho(t, l) dt dl$$

Because, $R_h(x, t, y, l) = R_h(x-t, y-l)$ and we want constant illumination on plane,

$$c = \int_{-W_l/2}^{W_l/2} \int_{-W_t/2}^{W_t/2} R_h(x_i - t, y_i - l) \rho(t, l) dt dl$$ 

Utilizing dirac-delta distribution yields next approximation,

$$\rho(t, l) = \sum_{i=1}^N \delta(t- t_i, l-l_i) =\sum_{i=1}^N \delta(t- t_i) \delta( l-l_i)$$

### 1 dim case

Suppose that the distribution $\rho(x,y)$ is separable. 
That is,

$$\rho(x,y) = \rho(x) \cdot \rho(y)$$

For separable input the higher dimension convolution is decomposed into products of 1-dim convolutions.

$$c  = \left(\int_{-W_l
/2}^{W_l/2} R_h(y_j - l) \rho(l) dl \right) \cdot \left(\int_{-W_t/2}^{W_t/2} R_h(x_i - t) \rho(t) dt \right)  $$

This integral equation is a Fredholem integral equation of first kind, since it has a difference kernel, $R_h$, by the choice of radiation pattern, Fourier transform method can be used to get solution. Unfortunately, centeral dominant radiation pattern like Lambertian or Gaussian, the solution function on frequncy domain diverges as $|\omega| \rightarrow \infty$ which is invaild Dirichlet's conditions. Therefore, we cannot get exact distribution formula from the equation. 

In this case, to get an approximation of system, we can matrize above system. 1-dim convolution equation could be directly transformed to simple linear equation with discretization.

Let, for $i = \N^+ \backslash \{0 \}$,
$$\Delta_n = \frac{W}{n}$$
$$x_i = t_i = -\frac{W}{2} + \frac{1}{2}\Delta_n + (i-1) \Delta_n$$

$$\mathbf{R}_{ij} = R_h(x_i - x_j) \\
\vec{t}_i = \text{power weight on } t_i \text{ point} \\
\vec{i}_i = \text{illumination on } x_i \text{ point}
$$

$$\mathbf{R}\vec{t} = \vec{i} = c \vec{1}$$
 
where $\vec{t} ,\vec{1} \in \mathbb{R}^n$, $c \in \mathbb{R}^+$, $\mathbf{R} \in M_{n \times n}(\mathbb{R^+})$ for $n \in \mathbb{N}^+$.

$c = 1$ below equations for convinence of computation.

$$\mathbf{R}\vec{t} = \vec{i} = c \vec{1}$$

The ideal solution for system is zero-one solution  which is $\vec{t}_i = \{0, 1\} \forall i \in[1, n]$. However, this is one example of NP-hard problem in mathematics programming which means no standard algorithm exists. Focusing on weight of each points in space, then the positivness, where $\vec{t}_i \geq 0,   \forall i \in[1, n]$, is a given resitrction and it allow us various attemptions. Most well-known optimization method is Non-negative least square method(NNLS). With NNLS method the positive approximation is easily achieved but, the solvability and positiveness of system must be discussed. 

#### **Positiveness of system**

The association function $R_h(x,t): \mathbb{R} \times \mathbb{R} \rightarrow \mathbb{R}$ with Lambertian radiation pattern is 

$$R_h(x,t) = \frac{1}{H} (1+ (\frac{x-t}{H})^2 )^{-(s/2 +1)}$$

This is an *inverse multiquadratics* and strictly positive definite function, when $\frac{s+2}{2} >1$ and by the definition of Lambertian, $s \leq 1$, this condition always achieved.

See details of positive definite functions from below reference.

>  G. Fasshauer, Meshfree Approximation Methods with MATLAB, Interdisciplinary mathematical sciences (World244 Scientific, 2007).

By the definitions, the system is solvable and there exist unique solution. There are some other examples of positive definite system. ...

Now the remained property is positiveness. Kaykobad studied conditions of positive solution of positive system. 

> M. Kaykobad, Positive solutions of positive linear systems, Linear Algebra and its Applications, Volume 64, January 1985, pp 133-140, doi:10.1016/0024-3795(85)90271-X



<table style="border-radius:8px;width:100%;">
<th style="text-align:center;background-color:rgb(0, 0, 0); color:white; border-top-left-radius: 10px;width:20%;">
Thm</th>
<th style="text-align:left;">
Positive solution of Positive Linear system </th>
<tr style="text-align:center;">
<td colspan="2">

For linear system, $\mathbf{K} \vec{\sigma} = \vec{\beta}$, 
where 
$$\mathbf{K} \in M_{n\times n}(\R^+), 
\text{ diag}(\mathbf{K}) \in (\R^{+}\backslash \{0\})^{n} \\
\vec{\sigma} \in \R^n, \\
\vec{\beta} \in (\R^{+}\backslash \{0\})^{n}$$

if next condition is satisfied, $\forall i, j = 1, 2, \dots, n$,

$$\beta_i > \sum_{j=1, j \neq i}^n \mathbf{K}_{ij} \frac{\beta_j}{\mathbf{K}_{jj}}$$


then, the system is invertible and  $\vec{\sigma} = \mathbf{K}^{-1} \vec{\beta} \in (\R^{+}\backslash \{0\})^{n}$.

</td>
</tr>

</table>


---

The system we considering is a positive system whose all elements are positive, and because of the association function, 

$$\text{diag}(\mathbf{R}) = max(\mathbf{R}) \cdot \vec{1}^t$$

For symmetric positive system with $\vec{\beta} = \vec{1}$, using Kaykobad's theory next condition guarantees the positive solution of the system,

$$1 \geq   max(\sum_{i=1, i\neq j}^n \frac{\mathbf{R}_{ij}}{\mathbf{R}_{ii}}, j=[1,n]) $$

In addition, with rectangular approximation of integration we can get next,

$$2 \cdot \int_{W/2n}^{W/2} \frac{1}{max(R_h(x))}R_h(x)  dx \geq   max(\sum_{i=1, i\neq j}^n \frac{\mathbf{R}_{ij}}{\mathbf{R}_{ii}}, j=[1,n]) \cdot \Delta_n$$

Combining above two inequalities yields next inequality for the dimension $n$ which provides the positive solution of system.

$$\Delta_n \geq  2 \cdot \int_{W/2n}^{W/2} \frac{1}{max(R_h(x))} R_h(x)  dx$$

Thus,

$$\frac{W}{2 \cdot \int_{W/(2n)}^{W/(2)} R_h(x)  dx} \geq n$$

For Lambertian case, we can get approximated dimension that guarantee that positive solution of system,

$$\int_{-t}^t \frac{1}{(1+x^2)^{s/2+1}} dx =(2 t) {}_2F_1 \left(\frac{1}{2};\frac{s+2}{2} ;  \frac{3}{2}; -t^2 \right)$$

${}_2 F_1$ is a Gausse hypergeometric function. 

Evaluating inequality, we get an inequality of the dimension of the system for positive solution,

$$f(s, \alpha, n) := \frac{1/2 + {}_2F_1 \left(\frac{1}{2};\frac{s+2}{2} ;  \frac{3}{2}; -\frac{\alpha^2}{4 n^2}\right)}{{}_2F_1 \left(\frac{1}{2};\frac{s+2}{2} ;  \frac{3}{2}; -\frac{\alpha^2}{4} \right)}$$

$$ f(s, \alpha, n) \geq n$$

The LHS of inequality, $f(s, \alpha, n)$, is monotonic increasing function with upper bound as $n$ increasing, (See `nmax.md`). Therefore, upper bound value $1$ yields next inequality.

$$n_{max} \leq \sup{f(s, \alpha, n) } = \frac{1.5}{{}_2F_1 \left(\frac{1}{2}, \frac{s+2}{2}, \frac{3}{2}; -\frac{W}{2H}^{2} \right)}$$

---

With active set method and Non-negative least square optimization, we can get positive optimizated solution for system of arbitary dimension, but as far from the dimension, $n_{max}$, nosiy signal will occured between dominant positions of source in solution. Practically, we want to find minimum number of sources generate uniform illumination on rectangular plane. 

The solution can be treated as space-weight signal.
Even we get a high nosied solution for higher dimension of system, passing thorough low pass filter with blocking value $[-f_{max}, f_{max}]$, where $f_{max} = \frac{2\pi}{W}$, provide us some constant solution for given optimization parameter no matter how dimension are changed.

### 2 dim case

The general case is 2 dimension convolution and it is the original problem we want to solve. Luckily, the convolution is one of linear transform. Therefore, we can construct some linear system corresponding to each convolution system,

$$\mathbf{A} * \mathbf{X} = \lambda \mathbf{1} \Leftrightarrow \mathbf{A}'  \vec{x} = \lambda \vec{1}$$


where, $\mathbf{A} \in M_{k \times k }(\R^+)$, $\mathbf{X} \in M_{n \times n}(\R^+)$, $\mathbf{1} \in (\R^+)^n$, and
$\mathbf{A'} \in M_{m \times m}(\R^+)$, $\vec{x} \in (\R^+)^m$ , $\vec{\mathbf{1}} \in (\R^+)^m$ for $\lambda >0, k = 2n-1, m = n^2$.


Since, positive kernel is implicated to convolution system the transformed matrix system is always solvable. However, in this document, the matrix representation allow NNLS optimization for convolution equation. See detiails in supplment document. 

We can get power weight map for the given region.

If we can allocating power to the dimension, it is enough to apply to design optimization algorithm.
For whom want to get a result for source distribution, we can use local quatization of power weight map for each non zero points.

1 dimension:



2 dimension:

We have 4 kinds of grid and two types for each kinds.

Triangular, Rectangular, Petagonal, and Hexagonal
In addition, $x, y$ axis scale parameters could be considered.


## Converting to source location

This section describes sampling methods especially *inverse method*. 
The power weight on the given region $W= W_x \times W_y$ could be viewed as a
probability distribution of random variable $X, Y \in [-W_x/2 W_x/ 2] \times [-W_y/2,W_y/2 ]$. With sampling method from the distribution we can obtain some desired 
source locations. We want some unique sampling from the distribution, however most of  sampling method depend on some uniform distribution data during the procedure.
The inverse method also uses uniform distribution sample while the sampling process, but distribute them uniformly on $[0, 1]$ we can get some exact and location fitting for $N$ number of location sample. 

General acknowledgements

> Luc Devroye, Non-Uniform Random Variate Generation, Springer-Verlag, 1986, DOI:10.1007/978-1-4613-8643-8, ISBN:1461386454 

Chebyshev approximation and 2dim inverse algorithm 
> Olver, Sheehan, and Alex Townsend. "Fast inverse transform sampling in one and two dimensions." arXiv preprint arXiv:1307.1223 (2013).


Starting from next equation,

$$\mathbf{R}\vec{t} = \vec{i} = c \vec{1}$$

with NNLS method or ordinary matrix system solve with Kaykobad's condition we get some power weight map from the given geometric and optical system.

The constant $c$ is an arbitary positive real number. Whatever it is does not affect to the system solver. This means that we can restrict the constant, $c$, without loss of generality of the solve process.

We treat the output result as a probability density function defined on geometric area $W = W_x \times W_y$ and sampling from the probability distribution determined the solution. There are many sampling techniques, in this paper use **inverse transform sampling**. 

The definition of probability density function is 

<table style="border-radius:8px;width:100%;">
<th style="text-align:center;background-color:rgb(0, 0, 0); color:white; border-top-left-radius: 10px;width:20%;">
Def</th>
<th style="text-align:left;">
Probability density function</th>
<tr style="text-align:center;">
<td colspan="2">

For continuous random variable $X$ and function $f$, $f$ is a *probability density* for $X$ that possesses the following properties:


$$f(x) \geq 0$$

$$\int_{-\infty}^\infty f(x) dx = 1$$

$$\int_{a}^b f(x) dx = P(a < X <b) $$

where $a,b$ are any two values of $x$ satisfying $a<b$.

</td>
</tr>

</table>



> Paul G. Hoel, Introduction to Mahtemtical statistics, 5th edition, Wiley, 1984, ISBN:978-0471890454, p. 41


By the property of our previous solution the positiveness is achieved. Let $\vec{t}_\lambda$ is a solution when $c = \lambda$ on equation (marked). We can calcualte the exact value of $c$ which makes the ouput to satisfy the density function.

$$c = \frac{1}{\Delta_n ||\vec{t_1}||_{L^1}} $$

$$ = \frac{n}{W ||\vec{t_1}||_{L^1}} $$

for 2dim solution

$$c = \frac{n \cdot m}{Wx \cdot Wy \cdot ||\mathbf{X}||_{L^1}}$$


Cumulative distribution function of the discrete density


$$F_x(x) = P(X \leq x) = \sum_{x_i \leq x} P(X= x_i) = \sum_{x_i \leq x} p(x_i) $$

$$F_{xy}(x, y) = P(X \leq x, Y \leq y) = \sum_{x_i \leq x} \sum_{y_i \leq y } p(x_i, y_i)$$

The 2-dim inverse transform method is well described the preprint written by Olver, and Alex.

> Olver, Sheehan, and Alex Townsend. "Fast inverse transform sampling in one and two dimensions." arXiv preprint arXiv:1307.1223 (2013).


---

>Optical design of a compact multispectral LED light source with high irradiance and uniformity

>Modeling the irradiation pattern of LEDs at short distances 

>https://sci-hub.ru/10.1002/9781119076230.ch5

>LED irradiance pattern at short distances 

>http://www.chenq.ecei.tohoku.ac.jp/common/item/pdf/ronbun/2022_5_K-Mochiki.pdf

>Optical design of a compact multispectral LED light source with high irradiance and uniformity 



From the non-parametric estimation method, Weighted kernel density estimation. 

$$\hat{f_k} = \frac{1}{n} \sum_{i} w_i K(X-x_i, h)$$

The binary solution, $\{x_j'\}_{j=1}^{n'}$, is 

$$\frac{1}{n} \sum_{i} w_i K(X- x_i, h) = \lambda_{n'}\sum_{j}^{n'}K'(X - x_i', h')$$

We get weighted samples from NNLS solution of 2D convolution. The inverse sampling method with uniform $Y$, on $0-1$ range.
The 2D inverse sampling method is investigated by *Chebfun2* and its related project, the python implmentation is presented in *Supplement Materials*. 