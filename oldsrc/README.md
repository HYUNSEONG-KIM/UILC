<script type="text/x-mathjax-config">
MathJax.Hub.Config({
    displayAlign: "center"
});
</script>


# UILC (**U**niform **I**llumination **L**ED arrangement **C**alculatiuon)
This repository is for the calculation of the uniform illumincation LED arrangement for target plane

This project is written with C and offer as C library

Core calculation is written with GSL(GNU SCIENTIFIC LIBRARY)

By the GNU General Public License which is the license of the gsl package, this library is distributed with GNU General Public License.
You can use, modify, redistribute to other without no financial barrier only if you also follow the same license; GPL.
 

---
# Uniform irradiation LED arrangement optimization using approximation of integral equation

Kim, Hyeaonsung

*School of physics and photonics, Gwangji Insititute of Scoence and Technology, 123 Chamdan St. Gwangju, Korea*

\* *qwqwhsnote@gm.gist.ac.kr*

**Abstract**

© 2021 Optical Society of America under the terms of the [OSA Open Access Publishing Agreement](https://www.osapublishing.org/library/license_v1.cfm)

## 1. Introduction

## 2. Mathematical Construction

Common LED intensity model is a Lambertian intensity distribution model[^Lamber] with an inverse square law.

$$I'(r, \theta) = \frac{I_0}{r^2}\cos^m(\theta)$$

Where $I_0$ is an intensity that LED emits perpendicular to its surface, $\theta$ is an angle between perpendicular vector and direction vector to point from center of the LED, and $m$ is a number, determing the optical property. 

![geometrc](./Fig/Fig0.PNG)
**Fig. 1**  Geometrical Schematic diagram of target plane and LED array plane.$(x_i,y_i,0)$ is a i-th LED location coordinate on below plane. $(x_j,y_j,0)$ is a sample point on target plane above LED plane.

### 2-1. 1-dimension array

Consider 1-dimensional array of LED and line with a distance $H$ from the array surface. For, i-th LED and j-th sample point, the distance $r$ and view-angle $\theta$ can be calculated as next with a location coordinate value $x,t$ of each line whose origins are located on same perpendicular line.

$$r_{ij} = \sqrt{H^2 + (x_i -t)^2}$$

$$\theta_{i} = \arctan(\frac{|x_i - t|}{H})$$

With a Lambertian distribution model, the intensity of the j-th sample point $I$ is

$$I'(t) := \sum_{i=1}^N I'(x_i,t)= I_0 \sum_{i=1}^N \frac{H^m}{(H^2 + (x_i - t)^2)^(\frac{m}{2} +1)}$$

for $N$ LEDs. To reduce the complexity of the calculation we consider all LEDs have same optical properties in radiation power and distribution.

The $I_0$ is differ as the power consumption condition of the LED and does not matter in this case so, we wikll define the $I(t)$ as next. 

$$I(t) =\sum_{i=1}^N I'(x_i,t)/I_0=  \sum_{i=1}^N \frac{H^m}{(H^2 + (x_i - t)^2)^(\frac{m}{2} +1)}$$


Now, assume that the LEDs are located in the closed region which length is $W$. We will only consider the intensity in this region on the above line. In other word, it means $x, t \in [-\frac{W}{2}, \frac{W}{2}]$ region. Using Kronecker delta function, we can represent the specific arrangement $\{ x_i \}_{i=1}^N$ of the LEDs as array function $\sigma$ and define the function $f(x,t)$ as next.

$$\sigma(x, \{ x_i \}_{i=1}^N) := \sum_{i=1}^N \delta(x-x_i)$$
$$f(x,t) := \frac{H^m}{(H^2 + (x - t)^2)^{(\frac{m}{2} +1)}}$$

With these definition, we can rewrite the intensity at $t$ point.

$$I(t) = f(x,t) \sigma(x, \{ x_i \}_{i=1}^N)$$

Our goal is finding an optimized arrangement $\{x_i \}_{i=1}^N$ which determines the array function $\sigma(x, \{ x_i \}_{i=1}^N)$.  However, the number of the LED $N$ is not important. We will concentrate on tendency of the $\sigma(x, \{ x_i \}_{i=1}^N)$ by $x$ values for general modification. Thus, we assume the continous function $\sigma_C(x)$ defined on $[-\frac{W}{2}, \frac{W}{2}]$. We can get $\sigma(x, \{ x_i \}_{i=1}^N)$ with discretizing $\sigma_C(x)$. 

$$D_N(\sigma_C(x)) = \sigma(x, \{ x_i \}_{i=1}^N)$$

Intuitively, we can infer some properties of the function $\sigma_C(x)$. For example, it will show a narrow distance between two consecutive elements $x_n$, $x_{n+1}$ near boundary than they are located near center point, because, when they have same distance, those array showed smaller intensity value near the boundary than central value. 

$$\sigma_C(x) \geq 0 , \forall x \in [-\frac{W}{2}, \frac{W}{2}]$$
$$\sigma_C(x) = \sigma_C(-x)$$
$$\sigma_C(x_1) \leq \sigma_C(x_2), |x_1| < |x_2|, x_1,x_2 \in [-\frac{W}{2}, \frac{W}{2}]$$


The initial step is determine the function $\sigma_C (x)$. We cannot construct it directly. Therefore, we can approximate using interpolation method with sample point set $\{x_j \}_{j=1}^n , \{ t_k \}_{k=1}^n$. These point set is discretizing the region $[-\frac{W}{2}, \frac{W}{2}]$ with $n$ number of interval points. 

$$x_j = (j-\frac{1}{2})\frac{W}{n}-\frac{W}{2}, t_k = (k-\frac{1}{2})\frac{W}{n}-\frac{W}{2}$$

For sample point $t_k$, the intsnsity $I_k = I(t_k)$ is

$$I(t_k) = \sum_{j=1}^n f(x_j, t_k) \sigma_C(x_j)$$

Therefore, we can construct linear system

$${\bf{F}} {\sigma_C} = {I}$$

,where $\bf{F} \in M_{n\times n}(\mathbb{R})$ and ${\sigma_C}, {I} \in \mathbb{R}^n$. Each elements are defined as

$${\bf{F}}_{jk} := f(x_j, t_k)$$
$${\sigma_C}_j := \sigma_C(x_j)$$
$${I}_k := I(t_k)$$

By the definition of the $f(x_j, t_k) = f(|j-k|\frac{W}{n})$, the matrix $\bf{F}$ is bisymmetric matrix and it is invertible **Cite_proof_invertible**. However, we want get positive solution which all elements of $\sigma_C$ are positive by the first property of $\sigma_C$ with constant vector ${I} = (I_1, I_2, \dots , I_n)^T, I_{i+1} = I_i \forall i$. 

M. Kaykobad studied the condition of the positive solution for specific linear system, whose diagnal elements of ${\bf{F}}$ and elemenst of ${I}$ are larger than zero and ${\bf{F}}$ is non-negative matrix. If for all $j,k = 1, 2, \dots n$,

$${I}_j > \sum_{k=1, k \neq j}^n \frac{{\bf{F}}_{jk}}{{\bf{F}}_{kk}} {I}_k $$

then, matrix ${\bf{F}}$ is invertible and solution ${\sigma_C} = {\bf{F}}^{-1}{I}$ is postive solution whose elements are all positive. 

Our system satisfies all pre-condition of the Kaykobad's Theorem. In addition, our matrix ${\bf{F}}$ has unit diagonal entry and all ${I}$ elements are same. With that, the condtion state the diagonal dominant property of ${\bf{F}}$. Define the number $n_{max}$ as the maximum number of sample points which have positive solution. The solution ${\sigma_C}>>0$ for $n\leq n_{max}$, and $\exist {\sigma_C}_k \leq 0$ for $n>n_{max}$.

The equation (diagonal dominant) can be rewritten as

$$ I_0 > \sum_{k=1, k \neq j}^n \frac{f(|j-k|\frac{W}{n})}{f(0)} I_0$$

Let's define normalized funtion $f_n(x)$.

$$f_n(x) = \frac{f(x)}{f(0)}$$

$max(f_n) = f_n(0) = 1$ and with a relationship of Reimann sum and integral we can show next,

$$\int_{a}^{a+W} f_n(x) \, dx \geq \sum_{k=1, k \neq j}^n f_n(|j-k|\frac{W}{n}) (\frac{W}{n})$$

$a \in [-\frac{W}{2}, \frac{W}{2}]$. 

**[Picture?]**



For specific $m$, the maximum case of the summation is the row that the diagonal element exist on middle of it. If matrix satiesfies the equation in maximum case, it holds other rows automatically. Such case is $a \approx - \frac{W}{2}$. Therefore, we can determine the $n_{max}$ as next.

$$ n_{max} = \lfloor {W}{/\int_{-\frac{W}{2}}^{\frac{W}{2}} f_n(x) \, dx} \rfloor = \lfloor 1/ _2F_1(\frac{1}{2}, \frac{m+2}{2}; \frac{3}{2}; - (\frac{W}{2H})^2) \rfloor$$

where $_2F+1$ is a Gaussian Hypergeometic function, or just directly test using middle row of the matrix with $p=n/2$ for even $n$ and $p = (n+1)/2$ for odd $n$. The $n_{max}$ is a maximum $n$ which satisfies next

$$1 > \sum_{k=1, k \neq p}^{n} f_n(|p-k|\frac{W}{n}) $$

When Lambertian model number $m$ increases, the distribution function $f(x,t)$ become more similar to Kronecter delta function. If $m$ decreases, the non-diagonal term become larger than before, so that the matrix easily loose diagonal dominant charateristic at smaller $n$.

**[Fig: m= 80, n= 20, n =50 comparsion]**

The problem arise, when $m$ is small. Such case, we cannot find enough sample points for $\sigma_C(x)$. Even though, $n$ overcome the $n_{max}$, we can obtain meaningful samples, if $n - n_{max}$ is small. However, the greater value of gap $n - n_{max}$, the solution becomes less meaningful.

**[m=3, $n$= $n_{limt}$, $n$>$n_{limt}$, $n$>>$n_{limt}$,] comparsion**

Therefore, when $m$ is too small that we cannot get meaningful number of samples, we have to use another way for solving the system. It is equivalent with Non-Negative Least Square(NNLS) problem that finding a positive solution of the linear system. *Active set method* is commonly used to solve NNLS problem. It is not only possible to find the exact solution with positive constraint, but also providing an approximation solution with cutting some precison by limited iteration numbers.

**[Comparsion with n < $n_{max}$ case, n > $n_{max}$ with inverse matrix and Active set method comparsion, solution and the $F$ multiplication]**

After we get the $\sigma_C$ solution, we need to transform it to discrete LED arrangement $\sigma$. We can assume the value of the $\sigma_C$ as the density of the LED. To determine the precison value of the density, we need some standard case. In this research, we will use Morena's analytic solution $\sigma_M$.
$$\sigma_M (x) = \sum_{i=1}^N \delta(x-x_i)$$

$$x_{i+1} - x_i = d_M(m)*H$$

and define the system as

$${\bf{F}} ({\sigma}_{MC} + \sigma_{add})= {I}$$

Where $I_i = I_M(0)$, and $\sigma_C = \sigma_{MC} + \sigma_{add}$. The $\sigma_{MC} = (\sigma_{MC}(x_1), \sigma_{MC}(x_2), \dots , \sigma_{MC}(x_n))^T$

$${\bf{F}} {\sigma}_{MC}= {I_M}$$

$I_M = (I_{M1},I_{M2}, \dots , I_{Mn})^T$, $I_{Mk} = I_M(t_k)$

With previous steps, we can calculate the ${\sigma}_{MC}$ and with a $\sigma_{M}(x)$. We can construct accurate discretization transform $D$ such that

$$D(\sigma_{MC}) = \sigma_M$$


Finally, applying the transform $D$ to $\sigma_{MC}$ and we get a optimized solution $\sigma(x)$
$$D(\sigma_{MC}) = \sigma(x)$$ 


### 2-2. 2-dimension Rectangular array

### 2-3. 2-dimension Cicular array
## 3. Experiment 

## 4. Result and discussion


## 5. Conclusion


## References

[^Lamber]: D. Wood, Optoelectronic Semiconductor Devices, Prentice-Hall international series in optoelectronics (Prentice Hall,
1994).



