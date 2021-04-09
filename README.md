# TIME DIFFERENCE OF ARRIVAL (TDOA) LOCALISATION

Time difference of arrival (TDOA) is a method for target localisation. Let the 2D position of the target be $\mathbf{x} = [x \; y]^{\mathsf{T}}$, and let $\{\mathbf{x}_{i}\}_{i=1}^{M}$ be the position of the $M$ sensors used to locate the target. The TDOA system measures difference in ranges using the time-delay of a test signal and the measurements are modelled as:

$$
    r_{i} = cT_{0} + d_{i} + w_{i}, \; i=1,\cdots,M,
$$

where $T_{0}$ is the time at which the target emits the beacon, $c$ is the speed of propagation of the wave, $d_{i} = \Vert \mathbf{x}-\mathbf{x}_{i} \Vert_{2}$ is the distance between the sensors and the target and $w_{i}$ is additive white Gaussian noise with uniform variance $\sigma^{2}$. In the case of this example, the number of sensors $M=4$. In TDOA, to eliminate $T_{0}$, the difference in the ranges is used, i.e.,

$$
    r_{ij} = d_{i} - d_{j} + n_{ij}, \; i,j=1,2,3,4,
$$

where $n_{ij} = w_{i}-w_{j}$. The differences in ranges are used for estimation of the unknown location of the target $\mathbf{x}$.

To estimate the 2D position of the target, three range measurements are sufficient. This implies three range-difference measurements are sufficient. We consider the following measurements for estimation:

$$
    \underbrace{\begin{bmatrix}
		r_{12} \\ r_{23} \\ r_{34}
	\end{bmatrix}}_{\mathbf{r}} = 
	\underbrace{\begin{bmatrix}
		d_{1}-d_{2} \\ d_{2}-d_{3} \\ d_{3}-d_{4}
	\end{bmatrix}}_{\mathbf{d}(\mathbf{x})} +
	\underbrace{\begin{bmatrix}
		n_{12} \\ n_{23} \\ n_{34}
	\end{bmatrix}}_{\mathbf{n}}.

$$

The noise $\mathbf{n}$ is a linear transformation of the i.i.d Gaussian random vector with zero mean and variance $\sigma^{2}$. The transformation is given by:

$$
    \mathbf{n} = \begin{bmatrix}
		w_{1} - w_{2} \\ w_{2} - w_{3} \\ w_{3} - w_{4}
	\end{bmatrix} = 
	\underbrace{\begin{bmatrix}
		1 & -1 & 0 & 0 \\
		0 & 1 & -1 & 0 \\
		0 & 0 & 1 & -1 \\
	\end{bmatrix}}_{\mathbf{A}}
	\begin{bmatrix}
		w_{1} \\ w_{2} \\ w_{3} \\ w_{4}
	\end{bmatrix}.

$$

Therefore, distribution of $\mathbf{n}$ is Gaussian with mean $\mathbb{E} \left[ \mathbf{n} \right] = \mathbf A \mathbb E \left[ \mathbf w \right] = \boldsymbol{0}$, and covariance matrix $\boldsymbol \Sigma = \mathbb E \left[ \mathbf n \mathbf{n}^{\mathsf T} \right] = \mathbf A \mathbb E \left[ \mathbf w \mathbf w^{\mathsf T} \right] \mathbf A = \sigma^{2} \mathbf A \mathbf A^{\mathsf T}$. Therefore, the joint data-distribution of the measurements $\mathbf r$ is also Gaussian, i.e., $\mathbf r \sim \mathcal{N} \left( \mathbf d(\mathbf x), \boldsymbol \Sigma \right)$.

## Cramer-Rao Lower Bound

Let $\hat{\mathbf x}$ be any unbiased estimator for $\mathbf x$. Then, using Cramer-Rao lower bound (CRLB) theorem:

$$
    \mathrm{cov}(\hat{\mathbf x}) \geq \mathbf I (\mathbf x)^{-1},
$$

where $\mathbf I (\mathbf x) = \left( \frac{\partial \mathbf d(\mathbf x)}{\partial \mathbf x} \right)^{\mathsf T} \boldsymbol \Sigma^{-1} \left( \frac{\partial \mathbf d(\mathbf x)}{\partial \mathbf x} \right)$ and therefore the variance of the estimates for the coordinates, $\mathrm{var}(\hat{x}) \geq [\mathbf I (\mathbf x)]^{-1}_{11}$ and $\mathrm{var}(\hat{y}) \geq [\mathbf I (\mathbf x)]^{-1}_{22}$.

## Maximum-Likelihood Estimation

The maximum-likelihood estimator (MLE) for the parameter $\mathbf x$ is given by:

$$
\begin{align}
    \hat{\mathbf x}_{MLE} = \mathrm{arg }\max_{\mathbf x \in \mathbb R^{2}} p(\mathbf r ; \mathbf x), \\
    = \mathrm{arg }\max_{\mathbf x \in \mathbb R^{2}} \ln p(\mathbf r ; \mathbf x), \\
    = \mathrm{arg }\min_{\mathbf x \in \mathbb R^{2}} J(\mathbf x) \doteq \frac{1}{2} (\mathbf r - \mathbf d(\mathbf x))^{\mathsf T} \boldsymbol \Sigma^{-1} (\mathbf r - \mathbf d(\mathbf x)).
\end{align}
$$