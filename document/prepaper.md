

# LED array design for uniform irradiance on rectangular area with boundary matching

Author : HYEON-SUNG KIM $^1$ <br>

$^1$ Department of Physics and Photon Science, GIST (Gwangju Institute of Science and Technology), Gwangju, Korea <br>
$^1$ qwqwhsnote@gm.gist.ac.kr <br>
$^1$ https://sites.google.com/view/hf3dlab/

## ABSTRACT

This paper studied design method for uniform irradiance intensity on linear and rectangular area in analytical approach. This method aims to enhance boundary intensity which rapidly decrease near boundary in previous papers. The complete process is combination of two analytical method; central uniformity and boundary reinforce propagation. Central design can be achieved with various way not only analytical but also stochastic methods.

## INTRODUCTION

Uniform irradiance light source has been widely studied and used for many industries and research fields. Since, Light Emitting Diode(LED) is developed, most optical devices have been replaced with LED device because of its high energy efficiency, fast reaction time, compact size and robustness, and long usage time \cite{khatib2009sustainability, kalyani2017principles}. However, because normal LED sources show spatial irradiation intensity distribution for viewing angle, it is hard to construct uniform areal intensity light device with multiple LED components. 

In such limitation, various attemptions have been studied by researchers. The target plane intensity flatness depends on three parameters of each LEDs in array; optical power, irradiation distribution, and location. Generally, restraining two of those parameters and optimize rest value. Practically, optimizing location of light source requires less components to construct real device from calculation, since, the other methods require lens design for change irradiation distribution \cite{Whang2009} or additional control unit to modulate power of each LEDs\cite{PraneethVarma:17}. For LED location, it is trivial that infinite uniform LED distribution shows flat irradiance intensity for any target plane. In addition, central area intensity induced by finite LED array will approximate flat intensity well if they are symmetric around center and uniformly located. However, minimizing number of LED for light device is critical consideration. As their number increasing, maintains and power efficiency are get worse, especially, for high power LED, heat flow is added to them. Morena et al \cite{Moreno:06}, suggested the analytic solution using generalized Lambertian LED model for linear, rectangular, circular array cases with expansion of Sparrow resolution criterion \cite{Barakat:62}. This study restricts LED power and irradiation distribution(;generalized Lambertian) and optimize LED array to find maximum LED distance between two LED in array.  
After, this analytical solution is introduced, various stochastic methods are studied for optimizing LED location. Zhoupig et al, first studied with stochastic method; stimulated annealing algorithm for different optical model LEDs for finite area\cite{Su:12}. Sourav Pal used genetic algorithm and evolutionary programming \cite{Pal:15}. Lei et al, tried basic local search method \cite{Lei:0}. Yu et al, suggested particle swarm optimization\cite{yu2018research}.

There are several differences between analytical solutions and stochastic solutions. For analytic solutions, all LEDs are assumed share same optical properties(intensity and Lambertian coefficient) and plane area is not determined at initial calculation. The total area is an outcome of calculation. Meanwhile in stochastic methods, plane area is given value for calculation and have less constraints for optical models of LEDs. They even can use various LED models for one calculations such as from different Lambertian value to general polynomial model\cite{}. These freedom of stocastic methods are useful for many case, however there are several disadvantages.  
Practically, there are many case which LED locations are only permitted for area equal or smaller than target plane. Even those case, expanded Sparrow's criterion method still useful. However, they shows rapid decreasing distribution near boundary. 



Since, results of previous papers show rapid decreasing near boundary of LEDs array. For example, ---- show xx\% area show under 60\% intensity for center point. 

In this paper, the analytical algorithm is introduced for computing LED array with given rectangular target area. This algorithm use center-boundary irradiance intensity matching in order to prevent rapid decreasing near boundary. 

## FORMULAS AND EQUATIONS

### Assumption 


In this calculation, it is assumed that permitted area for LED location is same as target area. It is only consider center point of LEDs, so the total LED device can overcome boundary as much as half of dimension of iteself. 

### Optical model of LED

<figure class = "image">
    <img src="./img/Total_img.svg" alt="{{hello}}">
    <figcaption>
    </figcaption>
</figure>

 The imperpect Lambertian or generalied Lambertian model, this model is commonly used for hemisphere lenz LED device and far field irradation. With inverse-sqaure law, this model can be written as next form

$$I(r,\theta) = \frac{I_0}{r^2} \cos^m(\theta)$$

where $\theta$ is a viewing angle the angle between target vector and LED orientation, $r$ is a distance between LED emitting point and target point, $I_0$ is irradiation flux value at $r = 1$m with $0^\circ$ viewing angle, and value $m$ is determine the shape of distribution of imperpect Lambertian model which affected by curvature of lenz.

For given distance $H$, between two parallel plane; LED plane and target plane, Eq() can be rewritten as next form

$$I(x,t) = \frac{I_0 H^m}{(H^2 + (t-x)^2)^{\frac{m}{2}+1}}$$

where $x$ is a location of LED and $t$ is location of target plane in same vertical plane as Figure ().

#### $\phi$ rotated LED




If LED orientation is not perpedicular to surface with $\phi$ angle for plane vector, then intennsity distribution Eq() is changed to Eq()

$$I(x,t,\phi)
= \frac{I_0}{(H^2 + (t-x)^2)^{\frac{m}{2}+1}} (H \cos (\phi) + (t-x)\sin(\phi))^m.$$ 

In this case, the maximum point location of distribution $t_{max}$ can be expressed as next form.

$$t_{max} = x + \frac{2mH \sin(\phi)}{\sqrt{(m+2)^2-(m-2)^2 \sin^2(\phi)} +(m+2)\cos(\phi)}$$




### Area Separation

The illuminous distribtuion of center and boundary point $t= 0, W/2$ induced by one pair of LEDs with distance $x$ from origin can be expressed as function $I_c(x), I_b(x)$. It is trivial that every LEDs are symmetric to z-axis. Therefore, it is enough to deal one pair of LEDs as one unit in calculation.

$$I_c(x) := I(x,0) + I(-x,0)  \\  \text{}\\  I_b(x) := I(x,\frac{W}{2}) + I(-x,\frac{W}{2})$$

The difference function of such two distributions is defined as 

$$\text{Di}(x):=(I_b - I_c)(x)$$

By defining two point $x_e, x_m$ with $\text{Di}(x)$ function, the given area can be separated with three region $P, Q, R$ as in Fig(). The $x_m$ is a root of $\text{Di}(x)$, and $x_m$ is a point such that $|\text{Di}(x_m)| = |\text{Di}(\frac{W}{2})|$. Since, the function $\text{Di}(x)$ is monotonically increasing on $[0, W/2]$ and $\text{Di}(0) < 0, \text{Di}(W/2) >0$, $x_e$ and $x_m$ exist and they are uniquely determined on region $[0, W/2]$. Additionally, for any point $x_q \in Q$ region, there exists point $x_r \in R$ region such that 

$$I_c(x_q) + I_c(x_r) = I_b(x_q) + I_b(x_r)$$

That means if LEDs of array only exist in $Q,R$ region with $\text{Di}$ corresponding relation. The intensity distribution will show

$$\sum_{i=1}^{n_{LED}} I(x_i,0) = \sum_{i=1}^{n_{LED}} I(x_i, W/2)$$



For large $m$ or large $W/H$ ratio value, function $\text{Di}(x)$ show almost horizontal line except near $x=0$ and $x = W/2$ points. In such situations, $x_m$ can be approximated by point such that $I(0,0) = 2 I(x,0)$ which is 

$$x_m \approx H\sqrt{2^{\frac{2}{m+2}}-1} $$

and $x_e$ is approximated half point of between $x_m$ and $W/2$ points.

This approximation is vaild if $20/H^2 > (1+ (W/H)^2)^{(m/2 +1)}$ which indicates condition that $I_b(W/2) \approx I(W/2, W/2)$. 

Pratically, if width of LED device $w$ is given, the domain of function $\text{Di}(x)$ must be set as $[w/2 , W/2]$, because for $x < w/2$, the two LEDs are overlaped each other and it is physically impossible. In the case which $x_m < w$, the given value $m, H ,W$ must be modified.

This propagation method can be combined with various method not only analytic methods but also stochastic methods such as previous simulated annealing, genetic algorithms, local search algorithm etc with finding flat condition for center with constraint LEDs are only located in *Q* region. With expanded Sparrow's criterion method [Morena](morena), it becomes fully analytical soltuion. 



### Linear 

 For given numerical value, $W, H, m$ the vaild condition can be acheived by Eq(). The one tendency of expanded Sparrow's criterion in Morena et al, LED distance decreases as number of LEDs increasing. Therefore, with $n = 2, 3$ cases 

$$dH \geq 2x_m $$ 
where N is even
$$dH \geq x_m$$ 
where N is odd

N =2 

$$d= \sqrt{\frac{4}{m+3}}$$

N = 3

$$d= \sqrt{\frac{3}{m+3}}$$

#### Even LEDs

### Linear - Odd LEDs

Since, the total LEDs are located symmetric for z-axis at origin, if the number of LEDs is odd, then the one LED consequnetly exists at $x=0$ point. If there are no LEDs at $|x| = x_m$ points, it does not matter. Additional two LEDs at points $|x| = W/2$ can make up for intensity gap between center and boundary induced by center LED. This boundary intension can not make gap zero but at least within small percentage $\epsilon$. For fixed $m$, we can modify $W/H$ ratio to achieve proper $\epsilon$ value using Eq().

$$\epsilon^{- \frac{2}{m+2}} -1 < (\frac{W}{H})^2$$

If LEDs exist at $|x| = x_m$ and origin, such origin LED must be treated as *P* region LED. If double intensity LED with same optical property is possible then, same.

Problem arises if one LED pair located near $|x| = x_m$ and there exists LED on origin point $x=0$. The $x_m$ represents the two LEDs at $|x| = x_m$ contribute almost same irradiance to center point with one LED below center.

### $P$ region LEDs 

The one of the problem of $P$ region LEDs except one which located on origin point $x=0$ is it is impossible to achieve almost same intensity on $x=W/2$ point, even if, the more LEDs are added to region $R$ with vertical direction. Although, the LED plane is expanded over $x=W/2$ point,  it is much harder as $m$ value increases. Therefore, LED orientation must be modified. If the orient of LED is rotated $\phi$ to outward direction as in Fig(), the intensity distribution can be expressed Eq(3) and maximum intensity point is far from LED location as $t_{max}$ in Eq(4)

The maximum value of $t_{max}$ is $t_{1/2}$ when $\phi = \pi/2$. 
Then, with constraint $x_{LED}+t_{max} = W/2$, added LEDs on region $R$  can raise boundary irradiation intensity minimizing loss of uniformity of intensity on region $R$ . The expansion area $W/2 <x < W/2+ \sqrt{m/2} H$ LEDs so do too. Since the concave of distribution function is positive for $|t_{max} + x_{LED}| < |x|$ and negative for $|x_{LED}| < |x| < |t_{max} + x_{LED}|$, the LEDs added to region $R$ are more.



### Rectangular 
The rectangular array can be acheived vertical overlapping of two linear array. Even though, their irradiation distribution of linear LED array by vertical line for LED array is different for single LED. The $m$ value increases as forming . As $m$ value increases, their vertical irradiation distribution will be remain same shape with single LED

## SIMULATION & EXPERIMENTAL Set Up

### Region A example 

### Including Region P example

### Compare with previous papers

* Rectangular 
* Linear
* 

## Results

## CONCLUSION
