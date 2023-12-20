# MRI ReconLab ShanghaiTech

## [Iterative SENSE](IterativeSENSE.ipynb)

Iterative SENSE solves linear problem $E\hat{m}=b$ as a least square minization $\min_m\|Em-b\|_2^2$. Applying the gradient descent method, the iterative SENSE algorithm is given by solving $\min_m\|(E^HE)m-E^Hb\|_2^2$. Where $E=UFC$

**Gradient Descent:**

Initialize:
$$\begin{aligned}
m_0&=E^Hb\\
r_0&=m_0-E^HEm_0\\
\end{aligned}$$
Iteration:
$$\begin{aligned}
\alpha_k&=\cfrac{r_k^Tr_k}{r_K^TE^HEr_k}\\
m_{k+1}&=m_k+\alpha_k r_k\\
r_{k+1}&=r_k-\alpha_kE^HEr_k\\
\end{aligned}$$

## [Compressed Sensing](CompressedSensing.ipynb)

SALSA (split augmented Lagrangian shrinkage algorithm) algorithm is a fast and robust algorithm for solving formulations of a unconstrained optimization problem with a non-smooth regularization term. SALSA is based on variable splitting and augmented lagrangian methods.
 
The CS reconstruction problem is formulated as follows:

$$\begin{aligned}
&\min\limits_{m,\theta} \|Em-b\|_2^2 + \lambda\Phi(\theta)\\
&\mathrm{s.t.}\quad m = \theta\\
\end{aligned} $$

Where $E=UF$ is the forward operator for image to k-space and undersampling, $b$ is the undersampled k-space data, $m$ is the image to be reconstructed, $\theta$ is the auxiliary variable introduced by the variable splitting method, $\Phi$ is the regularization term (Total Variation), and $\lambda$ is the regularization parameter.

The algorithm is given by the following steps:

**Initialization:**

$$\textit{choose}\quad \mu>0, m_0, \theta_0, d_0$$

$$\bar{b} = E^Hb$$

**Iteration**

$$\begin{aligned}
m'_k&=\theta_k+d_k\\
r_k&=\bar{b}+\mu m'_k\\
m_{k+1}&=\frac{1}{\mu}(I-E^HE)r_k\\
\theta'_k&=m_{k+1}-d_k\\
\theta_{k+1}&=\Psi_{r\Phi/\mu}(\theta'_k)\\
d_{k+1}&=d_k-\beta_{k+1}+\theta_{k+1}\\
\end{aligned}$$

Where $\Psi_{r\Phi/\mu}$ denotes shrinkage/thresholding function for denoiseing the image.