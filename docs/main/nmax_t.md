

## Hypergeometric function ${}_2 F_1$

Analytic continuation, integrated introduction

> https://core.ac.uk/download/pdf/82108003.pdf



Definition: 

$${}_2F_1(a, b, c; z) := \sum_{n=0}^\infty \frac{(a)_n(b)_n}{(c)_n} \frac{z^n}{n!}$$
converges on disk $|z|<1$ and $c \neq 0, -1, -2, ...$.

### Properties

1. ${}_2F_1(a, b, c; z) = {}_2F_1(b, a, c; z)$: by the definition
2. ${}_2F_1(a, b, b; z) = (1-z)^{-a}$

### Special cases

$|z|<1$

1. ${}_2F_1\left(\frac{1}{2}, 1, \frac{3}{2}; z^2\right) = \frac{1}{2z} \ln(\frac{1+z}{1-z})$
2. ${}_2F_1\left(\frac{1}{2}, 1, \frac{3}{2}; -z^2\right) =\frac{1}{z} \arctan(z)$
3. ${}_2F_1\left(\frac{1}{2}, \frac{1}{2}, \frac{3}{2}; z^2\right) = \frac{1}{z} \arcsin(z)$
4. ${}_2F_1\left(\frac{1}{2}, \frac{1}{2}, \frac{3}{2}; -z^2\right) = \frac{1}{z} \ln(z + \sqrt{1+z^2})$

### Derivation

$$\frac{d}{dz} {}_2F_1(a, b, c ;z) = \frac{a b}{c} {}_2F_1(a+1, b+1, c+1; z)\\
= \frac{a}{z} \left({}_2F_1(a+1, b, c ;z) - {}_2F_1(a, b, c ;z)  \right)
$$

See "Contiguous relations for ${}_2F_1$ hypergeometric series", for the general version of above formula of the derivation.


### Transformation

$${}_2F_1(a, b, c ;z) \\
= (1-z)^{-a} {}_2F_1(a, c-b, c ;\frac{z}{z-1})\\
= (1-z)^{-b} {}_2F_1(c-a, b, c ;\frac{z}{z-1})\\
= (1-z)^{c-a-b} {}_2F_1(c-a, c-b, c ;z)$$

Initial two: Pfaff transformation


## Lambertian radiation on plane

$$R_{Lamber}(r, \theta) = \frac{1}{r^2} \cos^s(\theta)$$

Ignoring $\phi$ variation then for source location $x$ and target point $t$, with distance $h$ between of them,

$$R_{Lamber}(r, \theta) = R_{Lamber}(x,t) = \frac{h^s}{(h^2 + (x-t)^2)^{(s/2+1)}}\\
= \frac{1}{h^2} \frac{1}{(1 + (\frac{x-t}{h})^2)^{(s/2+1)}}$$

Let $k = |x-t| $ then,

$$ R_{Lamber}(x,t) = R_{Lamber}(|x-t|) \\
= R_{Lamber}(k) = \frac{1}{h^2} (1+k^2)^{-(s/2+1)}$$


---

See similar and more general proves in

> https://www.sciencedirect.com/science/article/pii/S0021904508002116#sec3

Differentiation of 2F1 function

> https://www.sciencedirect.com/science/article/pii/S1110256X1200020X

For, $t \in \R^+\backslash \{ 0 \}, b \in \R^+, n \in \N^+ \backslash \{ 0 \}$,

$$\int_{-t}^t (1+x^{2n})^{-b} dx = 2 \int_0^t (1+x^{2n})^{-b} dx \\
= 2 t \cdot {}_2F_1 \left( \frac{1}{2n}, b, \frac{1}{2n}+1, -t^{2n} \right)$$

Proof:

$x= t \cdot s \rightarrow dx = t ds$

$$2 \int_0^t (1+x^{2n})^{-b} dx = 2t \int_0^1 (1-(-t^{2n}) s^{2n})^{-b} dx$$

$s^{2n} = l \rightarrow ds = \frac{l^{(1/2n) -1}}{2n} dl$

$$2t \int_0^1 (1-(-t^{2n}) s^{2n})^{-b} dx = \frac{2t}{2n} \int_0^1 \frac{l^{1/{2n} -1}}{(1-(-t^{2n}) l)^{b}} dl$$

Let, $c = \frac{1}{2n}+1$ and $a = \frac{1}{2n}$

$${}_2F_1(a, b, c; z) = \frac{\Gamma(c)}{\Gamma(b)\Gamma(c-b)} \int_0^1 \frac{t^{b-1} (1-t)^{c-b-1}}{(1-tz)^a} dt $$ 
where, $\mathbf{R}(c) > \mathbf{R}(b) >0$

and (See Gamma recurrence relationship of DMLF)

$$\Gamma(z+1) = z \Gamma(z)$$

Thus,

$${}_2F_1(a, \frac{1}{2n}, \frac{1}{2n} +1 ; z) = \frac{1}{2n} \int_0^1 \frac{t^{1/2n -1}}{(1-zt)^a} dt$$

and

$$\int_{-t}^t (1+x^{2n})^{-b} dx = 2 t \cdot {}_2F_1 \left( \frac{1}{2n}, b, \frac{1}{2n}+1, -t^{2n} \right) $$

---

$$\because  (1+x^{2n})^{-b_1} \geq (1+x^{2n})^{-b_2} > 0 \\ \forall x \in \R, \forall b_1, b_2 \text{ s.t } 0< b_1 < b_2$$

$$2\int_0^t (1+x^{2n})^{-b_1} - (1+x^{2n})^{-b_2} dx \geq 0$$

$$\therefore  {}_2F_1 \left( \frac{1}{2n}, b_1, \frac{1}{2n}+1, -t^{2n} \right) \geq {}_2F_1 \left( \frac{1}{2n}, b_2, \frac{1}{2n}+1, -t^{2n} \right)$$

$$1 \geq (1+t^{2n})^{-\frac{1}{2n}}\\ \geq {}_2F_1 \left(\frac{1}{2n}, b, \frac{1}{2n}+1; -t^{2n}\right)\geq 0 $$


* ${}_2F_1 \left(\frac{1}{2n}, b, \frac{1}{2n}+1; 0 \right) = 1$
* $\lim_{t->\infty} {}_2F_1 \left(\frac{1}{2n}, b, \frac{1}{2n}+1; -t^{2n} \right) = 0$


$$\frac{d}{dt}{}_2F_1 \left(\frac{1}{2n}, b, \frac{1}{2n}+1; -t^{2n}  \right) = -2n t \cdot {}_2F_1^{(1)} \left(\frac{1}{2n}, b, \frac{1}{2n}+1; -t^{2n}  \right)\\
= \frac{1}{t} \left(\left(1+t^{2 n}\right)^{-b} - _2F_1\left(b,\frac{1}{2 n};1+\frac{1}{2 n};-t^{2 n}\right)\right) \\
$$

$$\because   {}_2F_1(a, b, c ; z) = (1-z)^{-b} {}_2F_1(c-a, b, c ;\frac{z}{z-1}) $$

$$\frac{d}{dt}{}_2F_1 \left(\frac{1}{2n}, b, \frac{1}{2n}+1; -t^{2n}  \right) = \frac{1} {t(1+t^{2n})^{b}} \left( 1 - {}_2F_1 \left(1, b, \frac{1}{2n}+1; \frac{t^{2n}}{1+ t^{2n}}  \right) \right) $$

Since ${}_2F_1 \left(1, b, \frac{1}{2n}+1; \frac{t^{2n}}{1+ t^{2n}}  \right) > 1,  \forall t >0$, $f(t) = {}_2F_1 \left(\frac{1}{2n}, b, \frac{1}{2n}+1; -t^{2n}  \right)$ is monotonic decreasing $t >0$ and bounded $0 \leq {}_2F_1 \leq 1$.

Proof of ${}_2F_1 \left(1, b, \frac{1}{2n}+1; \frac{t^{2n}}{1+ t^{2n}}  \right) >1$, let $z = \frac{t^{2n}}{1+ t^{2n}}$ then $0<z<1, \forall t>0$

$${}_2F_1 \left(1, b, \frac{1}{2n}+1; z  \right) = \frac{\Gamma(\frac{1}{2n}+1)}{\Gamma(1) \Gamma(\frac{1}{2n})} \int_0^1 \frac{(1-s)^{\frac{1}{2n}-1}}{(1-zs)^b} ds$$

$$$$
$$\because \frac{\Gamma(\frac{1}{2n}+1)}{\Gamma(1) \Gamma(\frac{1}{2n})} = \frac{1}{2n}, (1-zs)^b < 1$$

That is, 

$$f(n) := \frac{0.5 + {}_2F_1 \left(\frac{1}{2}, \frac{s+2}{2}, \frac{3}{2}; -\frac{W}{2Hn}^{2}  \right)}{{}_2F_1 \left(\frac{1}{2}, \frac{s+2}{2}, \frac{3}{2}; -\frac{W}{2H}^{2} \right)} \leq \frac{1.5}{{}_2F_1 \left(\frac{1}{2}, \frac{s+2}{2}, \frac{3}{2}; -\frac{W}{2H}^{2} \right)}$$

$$n_{max} \approx \left\lfloor{\frac{1.5}{{}_2F_1 \left(\frac{1}{2}, \frac{s+2}{2}, \frac{3}{2}; -\frac{W}{2H}^{2} \right)}}\right\rfloor$$