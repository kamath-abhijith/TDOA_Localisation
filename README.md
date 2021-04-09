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

The maximum-likelihood estimator (MLE) for the parameter $\mathbf x$ is given by $\hat{\mathbf x}_{MLE} = \mathrm{arg }\max_{\mathbf x \in \mathbb R^{2}} p(\mathbf r ; \mathbf x)$:

$$
    \hat \mathbf{x}_{MLE} = \mathrm{arg }\min_{\mathbf x \in \mathbb R^{2}} J(\mathbf x) \doteq \frac{1}{2} (\mathbf r - \mathbf d(\mathbf x))^{\mathsf T} \boldsymbol \Sigma^{-1} (\mathbf r - \mathbf d(\mathbf x)).
$$

The MLE is given by the solution to the unconstrained optimisation programme, and we choose gradient descent with fixed step-size to solve the optimisation programme. The gradient descent updates with constant step-size $\alpha > 0$ are given by:

$$
    \hat{\mathbf x}^{(k+1)} = \hat{\mathbf x}^{(k)} + \alpha \left( \frac{\partial \mathbf d(\mathbf x^{(k)})}{\partial \mathbf x} \right)^{\mathsf T} \boldsymbol \Sigma^{-1} \left( \mathbf r - \mathbf d(\mathbf x^{(k)}) \right),
$$

The iterations are run until the error between successive updates is below tolerance or the maximum number of iterations is reached. The estimates of the coordinates of the target position from the MLE iterations is given directly by the entries of the vector.

## Best-Linear-Unbiased Estimation

The linearisation of the range measurements gives:

$$
	r_{ij}^{2} + d_{j}^{2} + 2r_{ij}d_{j} = d_{i}^{2} + e_{ij},
$$

where the noise is approximates as $e_{ij} = n_{ij}^{2} = (w_{i} - w_{j})^{2}$. Using $d_{i}^{2} = \Vert \mathbf x - \mathbf x_{i} \Vert_{2}^{2}$, we get:

$$
	r_{ij}^{2} - \Vert \mathbf x_{i} \Vert_{2}^{2} + \Vert \mathbf x_{j} \Vert_{2}^{2} = -2(\mathbf x_{i} - \mathbf x_{j})^{\mathsf T} \mathbf x - 2r_{ij}d_{j} + e_{ij},
$$

which is linear in the parameter $\mathbf x$. In matrix-vector form, and taking the parameter vector $\boldsymbol \theta = [x \; y \; d_{2} \; d_{3} \; d_{4}]^{\mathsf T}$, we have:

$$
	\underbrace{\begin{bmatrix}
		r^{2}_{12} - \Vert \mathbf x_{1} \Vert_{2}^{2} + \Vert \mathbf x_{2} \Vert_{2}^{2} \\
		r^{2}_{13} - \Vert \mathbf x_{1} \Vert_{2}^{2} + \Vert \mathbf x_{3} \Vert_{2}^{2} \\
		r^{2}_{14} - \Vert \mathbf x_{1} \Vert_{2}^{2} + \Vert \mathbf x_{4} \Vert_{2}^{2} \\
		r^{2}_{23} - \Vert \mathbf x_{2} \Vert_{2}^{2} + \Vert \mathbf x_{3} \Vert_{2}^{2} \\
		r^{2}_{24} - \Vert \mathbf x_{2} \Vert_{2}^{2} + \Vert \mathbf x_{4} \Vert_{2}^{2} \\
		r^{2}_{34} - \Vert \mathbf x_{3} \Vert_{2}^{2} + \Vert \mathbf x_{4} \Vert_{2}^{2} \\
	\end{bmatrix}}_{\boldsymbol \gamma} = 
	\underbrace{\begin{bmatrix}
		-2(x_{1} - x_{2}) & -2(y_{1} - y_{2}) & -2r_{12} & 0 & 0 \\
		-2(x_{1} - x_{3}) & -2(y_{1} - y_{3}) & 0 & -2r_{13} & 0 \\
		-2(x_{1} - x_{4}) & -2(y_{1} - y_{4}) & 0 & 0 & -2r_{14} \\
		-2(x_{2} - x_{3}) & -2(y_{2} - y_{3}) & 0 & -2r_{23} & 0 \\
		-2(x_{2} - x_{4}) & -2(y_{2} - y_{4}) & 0 & 0 & -2r_{24} \\
		-2(x_{3} - x_{4}) & -2(y_{3} - y_{4}) & 0 & 0 & -2r_{34} \\
	\end{bmatrix}}_{\mathbf H}
	\underbrace{\begin{bmatrix}
		x \\ y \\ d_{2} \\ d_{3} \\ d_{4}
	\end{bmatrix}}_{\boldsymbol \theta} +
	\underbrace{\begin{bmatrix}
		e_{12} \\ e_{13} \\ e_{14} \\ d_{23} \\ e_{24} \\ e_{34}
	\end{bmatrix}}_{\mathbf e}.
$$

The complete distribution of the random vector $\mathbf e$ is not important for estimation, and only mean and covariance of the distribution are sufficient. The mean of each entry in $\mathbf e$ is $2\sigma^2$ and the covariance matrix of $\mathbf e$ is given by:

$$
	\mathbf C = 2\sigma^{4} \begin{bmatrix}
		4 & 1 & 1 & 1 & 1 & 0 \\
		1 & 4 & 1 & 1 & 0 & 1 \\
		1 & 1 & 4 & 0 & 1 & 1 \\
		1 & 1 & 0 & 4 & 1 & 1 \\
		1 & 0 & 1 & 1 & 4 & 1 \\
		0 & 1 & 1 & 1 & 1 & 4 \\
	\end{bmatrix}.
$$

Then, the BLUE estimate for the parameter vector $\boldsymbol \theta$ is given by $\displaystyle \hat{\boldsymbol \theta} = \left( \mathbf H^{\mathsf T} \mathbf C^{-1}\mathbf H \right)^{-1} \mathbf H^{\mathsf T} \mathbf C^{-1} \left( \boldsymbol \gamma - 2\sigma^{2} \mathbb{1} \right)$. The estimates of the coordinates of the target position from the BLUE estimate for $\boldsymbol \theta$ is then $\hat{x}_{BLUE} = [\hat{\boldsymbol \theta}]_{1}$ and $\hat{y}_{BLUE} = [\hat{\boldsymbol \theta}]_{2}$.
