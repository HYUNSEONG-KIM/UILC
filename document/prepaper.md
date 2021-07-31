

# LED array design for uniform irradiance on rectangular area with boundary matching

Author : HYEON-SUNG KIM $^1$ <br>

$^1$ Department of Physics and Photon Science, GIST (Gwangju Institute of Science and Technology), Gwangju, Korea <br>
$^1$ qwqwhsnote@gm.gist.ac.kr <br>
$^1$ https://sites.google.com/view/hf3dlab/

## ABSTRACT


## INTRODUCTION


## FORMULAS AND EQUATIONS

### Assumption 


In this calculation, it is assumed that permitted area for LED location is same as target area. It is only consider center point of LEDs, so the total LED device can overcome boundary as much as half of dimension of iteself. 

### Optical model of LED

 The imperpect Lambertian or generalied Lambertian model, this model is commonly used for hemisphere lenz LED device and far field irradation. With inverse-sqaure law, this model can be written as next form

$$I(r,\theta) = \frac{I_0}{r^2} \cos^m(\theta)$$

where $\theta$ is a viewing angle the angle between target vector and LED orientation, $r$ is a distance between LED emitting point and target point, $I_0$ is irradiation flux value at $r = 1$m with $0^\circ$ viewing angle, and value $m$ is determine the shape of distribution of imperpect Lambertian model which affected by curvature of lenz.

For given distance $H$, between two parallel plane; LED plane and target plane, Eq() can be rewritten as next form

$$I(x,t) = \frac{I_0 H^m}{(H^2 + (t-x)^2)^{\frac{m}{2}+1}}$$

where $x$ is a location of LED and $t$ is location of target plane in same vertical plane as Figure ().




### Seperation of Area

The illuminous distribtuion of center and boundary point $t= 0, W/2$ induced by one pair of LEDs with distance $x$ from origin can be expressed as function $I_c(x), I_b(x)$. It is trivial that every LEDs are symmetric to z-axis. Therefore, it is enough to deal one pair of LEDs as one unit in calculation.

$$I_c (x) := \\ I_b(x) := $$

The difference function of such two distributions is defined as $$\text{Di}(x) := (I_b - I_c )(x)$$

By defining two point $x_e, x_m$ with $\text{Di}(x)$ function, the given area can be separated with three region *P*, *Q*, *R*. The $x_m$ is a root of $\text{Di}(x)$, and $x_m$ is a point such that $|\text{Di}(x_m)| = |\text{Di}(\frac{W}{2})|$. Since, the function $\text{Di}(x)$ is monotonically increasing on $[0, W/2]$ and $\text{Di}(0) < 0, \text{Di}(W/2) >0$, $x_e$ and $x_m$ exist and uniquely determined on region $[0, W/2]$.

For large $m$ or small $H$ value, it is hard to find these two points, because the function $\text{Di}(x)$ show almost horizontal graph except near $x=0$ and $x = W/2$ points. Such cases we can approximate $x_m$ as point which $I(0,0) = 2 I(x_m,0)$ which is 

$$x_m \approx H\sqrt{2^{\frac{2}{m+2}}-1} $$

and $x_e$ is half point of between $x_m$ and $W/2$ point.

This approximation is vaild if $20/H^2 > (1+ (W/H)^2)^{(m/2 +1)}$ which indicates condition that $I_b(W/2) \approx I(W/2, W/2)$. 

If width of LED device $w$ is given, the domain of function $\text{Di}(x)$ must be set as $[w/2 , W/2]$, because for $x < w/2$, the two LEDs are overlaped each other and it is physically impossible. However, practically 

This propagation method can be combined with various method not only analytic methods but also stochastic methods such as previous simulated annealing, genetic algorithms, local search algorithm etc with finding flat condition for center with constraint LEDs are only located in *Q* region. With expanded Sparrow's criterion method [Morena](morena), it becomes fully analytical soltuion. 



### Linear 


$$dH \geq 2x_m $$ 
where N is even
$$dH \geq x_m$$ 
where N is odd

#### Even LEDs

### Linear - Odd LEDs

Since, the total LEDs are located symmetric for z-axis at origin, if the number of LEDs is odd, then the one LED consequnetly exists at $x=0$ point. 

### A region LEDs

### Rectangular 
The rectangular array can be acheived vertical overlapping of two linear array. Even though, their irradiation distribution of 1 linear LED array by vertical line for LED array is different for single LED. The $m$ value increases as forming . As $m$ value increases, their vertical irradiation distribution will be remain same shape with single LED

## SIMULATION & EXPERIMENTAL EXAMPLES

### Region A example 

### Without Region A example

### Compare with previous papers

* Rectangular 
* Linear
* 

## CONCLUSION
