# Convolution of Multi-dimension

## Definition

**Milti-dimension convolution**

For $n$-dim function, $n \in \mathbf{N}+\backslash \{0\}$ $f,g: \mathbf{F}^n \rightarrow \mathbf{F}$, it convolution $h$ is defined as

$$h(\vec{v}) = (f * g)(\vec{v}) := \int_{\mathbf{F}^n} f(t)g(x-t) \,dt$$

### Discrete version


$$h(\vec{v}) = (f *^{(n)} g)(\vec{v})\\
= \sum_{t_n=-\infty}^{\infty} \sum_{t_{n-1}=-\infty}^{\infty} \cdots \sum_{t_1=-\infty}^{\infty} f(t_1, t_2, \dots, t_n) g(v_1 - t_1, v_2 - t_2, \dots, v_n - t_n)
$$

## Properties

1. Communtivity: $f*g = g*f$
2. Associativity: $f*(g *h) = (f * g) *h$
3. Distributivity: $f*(g+h) = f*g + f*h$
4. Associativity with scalar multiplication: $(\alpha f) *g = \alpha (f*g)$

### Linearity
For arbitary function $f$ on function space $\mathbf{F}$ defined on field $\mathbb{F}$, the operator $(f *)$ is linear operator on $\mathbf{F}$,

1. $\forall g, h \in \mathbf{F}$, $(f*)(g+h) = (f*)(g) + (f*)(h)$
2. $\forall g, \in \mathbf{F}, \alpha \in \mathbb{F}$, $(f*)(\alpha g) = \alpha(f*)(g)$

## Edge handling

Edge handling in convolution is a handling edge for calculating convolution near edge of the input, because near the edge of the given data there exist kernel parts requiring outside datas of the input. In mathematical definition and signal analysis, it is not mentioned before you treat discrete calculation. Commonly, this technique is significantly treated in image processing field:filtering.

### List of edge handling method

For kernel, $K$ and input $A$, where $K\in M_{5 \times 5}, A \in M_{n \times m}, n, m>5$, and $\epsilon_{ij}$ is a required value for calculating convolution outside of the input.

$$\begin{array}{c|c}
 {
    \begin{matrix}
    \epsilon_{11}, \epsilon_{12}, \epsilon_{13}\\ 
    \epsilon_{21}, \epsilon_{22}, \epsilon_{23}\\
    \epsilon_{31}, \epsilon_{32}, \epsilon_{33}
    \end{matrix}
 }&{
   \begin{matrix}
    \epsilon_{14}, \epsilon_{15}\\ 
    \epsilon_{24}, \epsilon_{25}\\
    \epsilon_{34}, \epsilon_{35}
    \end{matrix}
 } \\
 \hline
 {
   \begin{matrix}
    \epsilon_{41}, \epsilon_{42}, \epsilon_{43}\\ 
    \epsilon_{51}, \epsilon_{52}, \epsilon_{53}\\
    \end{matrix}
 }&{
   \begin{matrix}
    a_{11},  a_{12}\\ 
    a_{21},  a_{22}\\
   \end{matrix}
 }
\end{array}$$

* Extend:
$\begin{array}{c|c}
 {
    \begin{matrix}
    a_{11}, a_{11}, a_{11}\\ 
    a_{11}, a_{11}, a_{11}\\
    a_{11}, a_{11}, a_{11}
    \end{matrix}
 }&{
   \begin{matrix}
    a_{11}, a_{12}\\ 
    a_{11}, a_{12}\\
    a_{11}, a_{12}
    \end{matrix}
 } \\
 \hline
 {
   \begin{matrix}
    a_{11}, a_{11}, a_{11}\\ 
    a_{21}, a_{21}, a_{21}\\
    \end{matrix}
 }&{
   \begin{matrix}
    a_{11},  a_{12}\\ 
    a_{21},  a_{22}\\
   \end{matrix}
 }
\end{array}$

* Wrap: 
$\begin{array}{c|c}
 {
    \begin{matrix}
    a_{(n-2) (m -2)}, a_{(n-2) ( m-1)}, a_{(n-2) (m)}\\ 
    a_{(n-1) (m -2)}, a_{(n-1) ( m-1)}, a_{(n-1) (m)}\\
    a_{(n )(m -2)}, a_{  (n) ( m-1)}, a_{(n)(m)}
    \end{matrix}
 }&{
   \begin{matrix}
    a_{n-2 1}, a_{n-2 2}\\ 
    a_{n-1 1}, a_{n-1 2}\\
    a_{n 1}, a_{n2}
    \end{matrix}
 } \\
 \hline
 {
   \begin{matrix}
    a_{1 m-2}, a_{1 m-1}, a_{1m}\\ 
    a_{2 m-2}, a_{2 m-1}, a_{2m}\\
    \end{matrix}
 }&{
   \begin{matrix}
    a_{11},  a_{12}\\ 
    a_{21},  a_{22}\\
   \end{matrix}
 }
\end{array}$

* Mirror:
$\begin{array}{c|c}
 {
    \begin{matrix}
    a_{33}, a_{32}, a_{31}\\ 
    a_{23}, a_{22}, a_{21}\\
    a_{13}, a_{12}, a_{11}
    \end{matrix}
 }&{
   \begin{matrix}
    a_{31}, a_{32}\\ 
    a_{21}, a_{22}\\
    a_{11}, a_{12}
    \end{matrix}
 } \\
 \hline
 {
   \begin{matrix}
    a_{13}, a_{12}, a_{11}\\ 
    a_{23}, a_{22}, a_{12}\\
    \end{matrix}
 }&{
   \begin{matrix}
    a_{11},  a_{12}\\ 
    a_{21},  a_{22}\\
   \end{matrix}
 }
\end{array}$

* Constant:
For constant $c$, 
$\begin{array}{c|c}
 {
    \begin{matrix}
    c, c, c\\ 
    c, c, c\\
    c, c, c
    \end{matrix}
 }&{
   \begin{matrix}
    c, c\\ 
    c, c\\
    c, c
    \end{matrix}
 } \\
 \hline
 {
   \begin{matrix}
    c, c, c\\ 
    c, c, c\\
    \end{matrix}
 }&{
   \begin{matrix}
    a_{11},  a_{12}\\ 
    a_{21},  a_{22}\\
   \end{matrix}
 }
\end{array}$


* **Cropping**: Usually, cropping input data or kernel also discussed in edge handling topic. However, cropping method can be combined with the above methods or the kernel cropping can be implemented with constant method, therefore it will be treated in further section as preprocessig option in convolution. 


## Matrix representation of convolution

By its linerity, all discrete convolution can be represented with matrix mulitplication, 
even if dimension of the given system is higher than 2 by allowing larger rows and column number. 
Inspite of that, using tensor equation is recommanded for convenience. By the way, common matrix convolution is represented as next.

$$\mathbf{A} * \mathbf{X} = \mathbf{X} * \mathbf{A} = \mathbf{B}$$

where, $\mathbf{A} \in M_{n \times m}(\mathbf{F}), \mathbf{X} \in M_{l \times k}(\mathbb{F}), \mathbf{B} \in M_{p \times q}(\mathbb{F})$. 

However, the above is just represent 2-dim functions and their sample data as matrix form. **It is not a matrix multiplication representation**.
Matrix representation system corresponding to the above convolution is,

$$\mathbf{A_c} \cdot \mathbf{X_c} = \mathbf{B_c}$$

Allowing some transformation on $\mathbf{B_c}$, we can choose various matrix equation, matrix-matrix or matrix-vector forms can be presented by the representation. 
No matter the shape of system, $\mathbf{A_c}$ is always Toeplitz matrix or at least blocked Toeplitz matrix. 
About $\mathbf{X_c}$ and $\mathbf{B_c}$, they may need reshaping functions respectively by the representiation. 
For example, Michal and Krystian suggested $\mathbf{B_c} = \mathbf{B}$ transform in 2-dim convolution calculation. 

> Michal Gnacik, Krystian Lapa, Using Toeplitz matrices to obtain 2D convolution, posted: October 27th, 2022, doi:https://doi.org/10.21203/rs.3.rs-2195496/v1

In this document, two matrix-vector form transformations are presented and solvability condition of convolution system in matrix equation also discussed. 


### Dimension of convolution and cropping

The convolution of matrix in this document is defined as,

$$\begin{pmatrix}
x_{11} & x_{12} & \cdots &x_{1m}\\
x_{21} & x_{22} & \cdots &x_{2m}\\
\vdots & \vdots & \ddots& \vdots\\
x_{n1} & x_{n2} & \cdots & x_{nm} 
\end{pmatrix} 
* 
\begin{pmatrix}
h_{11} & h_{12} & \cdots &h_{1k}\\
h_{21} & h_{22} & \cdots &h_{2k}\\
\vdots & \vdots & \ddots& \vdots\\
h_{l1} & h_{l2} & \cdots & h_{lk} 
\end{pmatrix} = \mathbf{Y}$$

$$\mathbf{Y}_{ij} = 
\sum_{p=1}^{l}\sum_{q=1}^{k} x_{(i+p+1)(j+q+1)} h_{pq}$$

where, $1 \leq p \leq n-l +1, 1\leq q \leq m-k+1$.

The dimension of convolution is 

$$(n,m) * (l, k) = (|n-l|+1, |m-k|+1)$$

This method may be seems to be special case of convolution that does not allowing outside overlapping. 
However, if the input matrix is expanded to considering dimension and edge handling type, all convolutions can be transformed to above definition.
Before the expansion description, details of cropping are needed.

#### Cropping

Cropping is defining a minimum overlapping dimension between kernel and input data.
$(c_r, c_c)$ tuple indicates row and column cropping level respectively. In addition, $1 \leq c_r \leq l, 1\leq c_c \leq k$ for dimension of kernel $H \in M_{l \times k}$.
Default convolutions defined on signal and mathematics regard $(c_r, c_c) = (1, 1)$.
Next example shows overlapping of kernel and data which calculate $[1,1]$ element of convolution for $(1,1)$ and $(2, 3)$ cropping condition respectively,

$$\begin{array}{cccc|c}
  c& c& c& c& c\\
  c& c& c& c& c\\
  c& c& c& c& c\\
  c& c& c& c& c\\
  \hline
  c& c& c& c& c*a_{11}
\end{array}, 
\begin{array}{cc|ccc}
 c& c& c& c& c\\
 c& c& c& c& c\\
 c& c& c& c& c\\
 \hline
 c& c& c*a_{11}& c*a_{12}& c*a_{13}\\
 c& c& c*a_{21}& c*a_{22}& c*a_{23}\\
\end{array}$$

The outside of the input data usually omitted, $0$ element unless some edge handling method is defined.
Considering cropping coefficients $c_r, c_c$, overall dimension of convolution equation is

$$(n, m) * (l, k) = (n+l -2c_r +1, m+k -2 c_c +1) \\
=( (n +2l-2c_r ) -l +1, (m +2k-2c_c ) -k +1)\\
= ( n' -l +1, m' -k +1)$$

That is,

$$X * H, (c_r, c_c) = X' * H, (1,1)$$

$$X' := \left[
\begin{array}{c c c}
{} & R_1 &{}\\
\hline
C_1 &  X & C_2\\
\hline
{}& R_2 & {}
\end{array}\right]$$

where, $\dim(X) = (n ,m), \dim(H) = (l, k), \dim(X') = \left(n +2(l-c_r) , m +2(k-c_c)\right)$. 
The rectangle matrix $R_1 R_2, C_1, C_2$ are symmetric for each row and column filp operation and determined by edge handling method.

For example, see details of next convolution with cropping $(2, 2)$ and kernel cropping: $0$ constant handling.

$$C = \left(\begin{bmatrix}
a & b & c & d\\
e & f & g & h\\
i & j & k & l
\end{bmatrix} * \begin{bmatrix}
1 & 2 & 3\\
4 & 5 & 6\\
7 & 8 & 9
\end{bmatrix}\right) \in M_{(3) \times (3)}$$

The dimension of convolution is $3 = 3+3 -2*2+1, 4 = 4+3 -2*2 +1$ and some elements are weight combination of kernel covered area of input data.

* $C_{1,1} = 5*a + 6*b+ 8*e+ 9*f$
* $C_{1,2} = 4*a + 5*b+ 6*c+ 7*e + 8*f + 9*g$
* $C_{2,1} = 2*a + 3*b+ 5*e+ 6*f + 8*i + 9*j$
* $C_{3,3} = 1*g +2*h + 4*k + 5*l$

And this is same with next convolution with $(c_r, c_c) = (l, k) = (3, 3)$

$$C = \left(\begin{bmatrix}
0& 0 & 0 & 0 & 0& 0\\
0& a & b & c & d& 0\\
0& e & f & g & h& 0\\
0& i & j & k & l& 0\\
0& 0 & 0 & 0 & 0& 0
\end{bmatrix} * \begin{bmatrix}
1 & 2 & 3\\
4 & 5 & 6\\
7 & 8 & 9
\end{bmatrix}\right) \in M_{(3) \times (3)}$$


**List of Methods**

* $X * H \rightarrow \mathbf{H} \cdot \vec{x}$
* $X * H \rightarrow \mathbf{X} \cdot \vec{h}$
* $X * H \rightarrow \mathbf{H} \cdot \mathbf{X}$


## Method 1: $X * H \rightarrow \mathbf{H} \cdot \vec{x}$ 


$\mathbf{H} \in M_{(n-l+1)(m-k+1) \times nm} (\mathbb{F}), \vec{x} \in \mathbb{F}^{nm}$


$$X = \left[\begin{array}{}
\vec{x_1} \\
\hline \vec{x_2} \\
\hline \vdots \\
\hline \vec{x_{n-1}}\\
\hline \vec{x_{n}}
\end{array}\right], \vec{x} = \left[\begin{array}{}
\vec{x_1}^t \\
\hline \vec{x_2}^t \\
\hline \vdots \\
\hline \vec{x_{n-1}}^t\\
\hline \vec{x_{n}}^t
\end{array} \right]$$

$$ \vec{x_i} = [ x_{i1}, x_{i2}, \dots, x_{i (m-1) }, x_{i m}]$$

$$\mathbf{H} =\left[ \begin{array}{c|c|c|c}
H_{11}     & H_{12} & \cdots & H_{1,n} \\
\hline
\mathbf{0} & H_{11} & \cdots & H_{1,n-1}\\
\hline
\vdots & \vdots& \ddots & \vdots\\
\hline
\mathbf{0} & \mathbf{0} & \cdots &  H_{11}
\end{array}\right]$$

$$(H_{1i})[p, q] = \begin{cases} 
 h_{i, q-p +1} &  h_{i, q-p +1} \in H \\
 0 & q-p <0 | q-p\geq k
\end{cases} $$

$H_{1i} \in M_{(m-k+1) \times m}$ thus, $\mathbf{H} \in M_{((|n-l|+1)(|m-k|+1))\times(nm)}$

For example, $H_{1i}$ is

$$\begin{bmatrix}
h_{i1} & h_{i2} & h_{i3}&\cdots & h_{ik}   & 0         & \cdots & 0 &\cdots & 0\\
0      & h_{i1} & h_{i2}&\cdots & h_{ik-1} & h_{ik}    & \cdots & 0 &\cdots & 0\\
0      & 0      & h_{i1}&\cdots & h_{ik-2} & h_{ik-1}  & \cdots & 0 &\cdots & 0\\
\vdots & \vdots & \vdots& \ddots& \vdots   & \vdots    & \ddots & {}&\vdots & {}\\
0      & 0      & 0     & \cdots& h_{i1}   & h_{i2}    & \cdots & {h_{ik-1}}& \cdots &0 \\
0      & 0      & 0     & \cdots& 0        & h_{i1}    & \cdots & {h_{ik-2}}& \cdots &0 \\
\vdots & \vdots & \vdots& \ddots& \vdots   & \vdots    & \ddots & {} & {\vdots} & {} \\
0      & 0      & 0     & \cdots& 0        & 0         & \cdots & h_{i1} &\cdots & h_{ik}
\end{bmatrix}$$


#### Condition for square matrix

<table style="border-radius:8px;width:100%;">
<th style="text-align:center;background-color:rgb(0, 0, 0); color:white; border-top-left-radius: 10px;width:20%;">
Thm</th>
<th style="text-align:left;">
Square condition of transformed matrix of convolution </th>
<tr style="text-align:center;">
<td colspan="2">

For given convolution,

$$ X * K$$

where $X \in M_{n\times m},  K \in M_{l \times k}$ with cropping coefficient $(c_r, c_c)$,

corresponding matrix equation with () tranformation, 

$$\mathbf{K} \cdot \vec{\mathbf{x}}.$$ 

If $l = 2 c_r -1$ and $k = 2 c_c -1$ and kernel cropping then, the system is a square system of which $\mathbf{K} \in M_{nm \times nm}$

</td>
</tr>

</table>

**Proof**

Start from extend matrix, $X_{ext}$ by the kernel and cropping dimensions,

Note: It is a roundabout way but very convinence in flow.

$$X_{ext} := \left[
\begin{array}{c c c}
{} & R_1 &{}\\
\hline
C_1 &  X & C_2\\
\hline
{}& R_2 & {}
\end{array}\right]$$


$$X_{ext} \in M_{n_{ext}\times m_{ext}} \\n_{ext} = n + 2 e_r, m_{ext} = m + 2 e_c$$

where, $e_r = l- c_r, e_c = k -c_c$.

The corresponding matrix system is 

$$\mathbf{K} \in M_{a\times b}\\
a = (n_{ext}-l+1)(m_{ext}-k+1)\\
b = n_{ext} m_{ext}$$

First with kernel cropping, the additional $e_r = e_c =0$, thus $n_{ext} = n , m_{ext} = m$.

and if we extend $a$ with $n, m, l, k, c_r, c_c$ then,

$$a = nm \\
+ \left( (n+m) + (nk + ml) - lk +1 \right)\\
- \left( 2(n c_c + m_cr) -(l c_c + kc_r) + (c_r + c_c)\right)$$

$$\because l = 2c_r -1, k = 2c_c -1 \\
(n+m) + (nk + ml) - lk +1 \\
 = (n+m) +2(nc_c + mc_r ) - (n+m) + (l c_c + kc_r) - (c_r + c_c) \\

 \therefore a = nm, a=b
$$

---

if $l = 2n-1 \& k = 2m-1$ then, all submatrices $H_{1,i}$ become upper triangle Toeplitz matrix, thus we get square matrix system.

Determinant: $\Pi_{i=1}^{n} h_{i1}^{(m)}$




## 2d convolution with numpy

Data 의 각 행을 filter의 행으로 numpy.convolve(mode="vaild")로

Data[0, n-l+1 +0] * filter 0  번째 행 = (n-l+1)(m-k+1)
Data[1, n-l+1 +1] * filter 1  번째 행 = (n-l+1)(m-k+1)
Data[2, n-l+1 +2] * filter 2  번째 행 = (n-l+1)(m-k+1)
Data[3, n-l+1 +3] * filter 3  번째 행 = (n-l+1)(m-k+1)

\+= 2d convolution, mode="vaild" 와 같음

```.{py}
def convolve2d(  
    data:np.ndarray, 
    filter:np.ndarray,
    crop:Tuple[int, int]=(1, 1),
    edge_mode:Literal["extend", "wrap", "mirror","constant"]="constant",
    edge_params = [0]
    ):
    
    l, k = filter.shape
    er, ec = get_dim_ext((l, k), crop)

    data_ext = _expand_matrix(data, (er, ec), [edge_mode]+edge_params)
    n,m = data_ext.shape

    result = np.zeros(shape = (n-l+1, m-k+1))
    for i, h_row in enumerate(filter):
        data_i = data_ext[i:n-l+i+1]
        result += np.stack([np.convolve(r, h_row, mode="vaild") for r in data_i])

    return result
```