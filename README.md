# CapacitorCalculationsExtensionModule
 Extension module for python to calculate interdigital capacitor related things

# Install
In a terminal type:

`pip install [path to git repos]\CapacitorCalculations\CapacitorCalculations`

If there is a permission error add the `--user` tag to the end of the command.

# Documentation
This package includes functions which make calculating the capacitance for interdigital capacitors possible.

#### References

(1) Huey-Daw Wu, Zhihang Zhang, F. Barnes, C. M. Jackson, A. Kain and J. D. Cuchiaro, "Voltage tunable capacitors using high temperature superconductors and ferroelectrics," in IEEE Transactions on Applied Superconductivity, vol. 4, no. 3, pp. 156-160, Sept. 1994, doi: 10.1109/77.317831.

(2) S. S. Bedair and I. Wolff, "Fast, accurate and simple approximate analytic formulas for calculating the parameters of supported coplanar waveguides for (M)MIC's," in IEEE Transactions on Microwave Theory and Techniques, vol. 40, no. 1, pp. 41-48, Jan. 1992, doi: 10.1109/22.108321.

### Functions
`ellint_ratio(k: float) -> float`

`k`: modulus $k$ for complete elliptic integral of the first kind.

This is used to calculate $\frac{K(k)}{K'(k)}$ where $K$ is the complete elliptic integral of the first kind and $K'$ is its complement such that $K'(k)=K(k')$ and $k'=\sqrt{1-k^2}$. This uses the `comp_ellint1` from the cmath standard library. $K(0)=\frac{\pi}{2}$, but $\lim\limits_{k\rightarrow 1}K(k)\rightarrow \tilde{\infty}$, so for small $k$, $K'(k)$ becomes problematic. This starts to happen for $k<10^{-7}$. For these cases of small $k$, you should use `ellint_ratio_aprrox`.

---

`ellint_ratio_approx(k: float) -> float`

`k`: modulus $k$ for complete elliptic integral of the first kind.

This is used to calculate $\frac{K(k)}{K'(k)}$ for small $k$ ($k<10^{-6}). This is based on the approximation im ref. (2).

$\pi/ln[2\frac{1+(1-k^2)^{1/4}}{1-(1-k^2)^{1/4}}]$

---

`k_thin(g: float, h: float, u: float) -> float`

This is used to calculate the modulus $k$ for the complete elliptic integral of the first kind $K(k)$ 
all three values `g`, `h`, and `u` must be in the same length units.

`g`: the "gap distance" i.e. the spacing between fingers

`h`: the height i.e. the thickness 
