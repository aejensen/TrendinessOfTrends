## Data-driven non-stationarity (better title needed)
By construction, the distribution of the derivative at time $t$ conditional on the value of the function at the same time is a linear function, and the latent dynamics is therefore governed by linear first-order differential equation with time-varying coefficients that can be written as
\begin{align*}
df(t) \mid f(t), \beta, \theta = a(t) + b(t)f(t) + \epsilon(t)
\end{align*}
where $b(t) = \Cov[df(t), f(t) \mid \beta, \theta]\Cov[f(t), f(t) \mid \beta, \theta]^{-1}$.

[\textbf{TODO:} Fix from here] The generating model described hitherto depends on a parametric covariance function $C_\theta(\cdot, \cdot)$ governing the dynamics of the latent function. Usually in Gaussian Process Regression a stationary covariance function, such as the SE function in Equation (\ref{eq:seKernel}) is chosen out of convenience [@rasmussen2003gaussian]. Choosing a stationary covariance function, however, implies an important model assumption. It can be shown that for all stationary covariance functions $C_{\theta}(s,t) = C_\theta(|s-t|)$, it follows that $\partial_2 C_\theta(t,t) = 0$ and therefore $f(t) \!\perp\!\!\!\perp df(t)$ [@cramer1967stationary].
Using a model with a stationary covariance function implies that $b(t) = 0$, and the model must therefore require non-stationarity in order to have non-trivial latent dynamics with $b(t) \ne 0$.

It can be argued that a priori selecting a parametric non-stationary covariance function that fits data is a very difficult problem. We now present an extension to the previous model that facilitates a data-driven adaptation to non-stationarity in $f$. We do this by introducing a latent time-varying categorical variable $G(t)$ taking values in the discrete set $\left\{1, \ldots, K\right\}$ at every time $t$. We then assume that $f$ conditionally on $G(t) = k$ is a Gaussian Process with a stationary covariance function indexed by a parameter $\theta_k$. This corresponds to the following functional latent class model
\begin{align}
w(t) \mid \phi &\sim \mathcal{S}_{K-1}^\infty(\phi)\label{eq:latentClassModel}\\
G(t) \mid w(t), \phi &\sim \Multinomial(w(t))\nonumber\\
f(t) \mid G(t) = k, \theta_k, w, \phi &\sim \mathcal{GP}(m, C_{\theta_k}(\cdot,\cdot))\nonumber
\end{align}
where $w$ is a $K$-dimensional probability function and $\mathcal{S}_{K-1}^\infty$ is a distribution on the interior of the functional $K-1$ dimensional simplex
\begin{align*}
\left\{w_k(t), k = 1,\ldots, K, t \in \mathcal{T}, w_k(t) > 0, \sum_{k=1}^K w_k(t) = 1\right\}
\end{align*}
Marginalizing out the distribution of $G(t)$ in $f(t) \mid G(t) = k, \theta_k, w, \phi$ we may write the joint distribution of $f$ and $df$ conditional on the probability weights in the samme form as in Equation (\ref{eq:latentJoint}) but with the following covariance and cross-covariance functions
\begin{align}
C(s, t) &= \sum_{k=1}^K w_k(s)C_{\theta_k}(s,t)w_k(t)\label{eq:latentClassCov}\\
\partial_1 C(t, s) &= \sum_{k=1}^K dw_k(s)C_{\theta_k}(s,t)w_k(t) + \sum_{k=1}^K w_k(s) \partial_1 C_{\theta_k}(s,t)w_k(t)\nonumber\\  
\partial_2 C(s, t) &= \sum_{k=1}^K w_k(s)C_{\theta_k}(s,t)dw_k(t) + \sum_{k=1}^K w_k(s) \partial_2 C_{\theta_k}(s,t)w_k(t)\nonumber\\
\partial_1\partial_2 C(s, t) &= \sum_{k=1}^K dw_k(s)C_{\theta_k}(s,t)dw_k(t) + \sum_{k=1}^K w_k(s)\partial_1 C_{\theta_k}(s,t)dw_k(t)\nonumber\\
&+  \sum_{k=1}^K dw_k(s)\partial_2 C_{\theta_k}(s,t)w_k(t) + \sum_{k=1}^K w_k(s)\partial_1\partial_2 C_{\theta_k}(s,t)w_k(t)\nonumber
\end{align}
The construction with a latent class implies that the covariance function of $f$, $C$, becomes a time-varying convex sum of stationary covariance functions with separate parameters. The expressions for the covariance function of $df$ and its cross-covariance with $f$ follows according to Equation (\ref{eq:latentJoint}) and the chain rule. From the expression of cross-covariance it can be seen that $\partial_2 C(t,t)$ is generally non-zero due to the existence of the terms in the first sum. The terms in the second sum will, however, still equal to zero on the diagonal when each class uses a stationary covariance function. The joint posterior distribution for this extended model is then given similarly to the expressions in Equation (\ref{eq:jointPost}) with the appropriate covariance substitutions.

To model the functional class probabilities, $w(t)$, in Equation (\ref{eq:latentClassModel}) we suggest mapping $K-1$ linear basis expansions onto the simplex. The basis could for example be polynomial splines where we require that the degree is chosen such they are at least once continuously differentiable on $\mathcal{T}$. Let $\left\{B_p(t) : t \in \mathcal{T}, p = 1,\ldots,P\right\}$ be a set of basis functions on $\mathcal{T}$ with $P$ degrees of freedom and $(\phi_0^{(k)}, \ldots, \phi_P^{(k)})$ a vector of coefficients for each $k = 1, \ldots, K-1$ expansion. We then construct the $K-1$ functions as
\begin{align}
\eta_k(t) = \phi_{0}^{(k)} + \sum_{p=1}^P \phi_{p}^{(k)} B_p(t), \quad k = 1,\ldots, K - 1\label{eq:splineExpansion}
\end{align}
which we map into the class probabilities by the inverse log-additive ratio transform
\begin{align*}
w_k(t) = \frac{\exp(\eta_k(t))}{1 + \sum_{j=1}^{K-1} \exp(\eta_{j}(t))}, \quad k = 1,\ldots, K-1, \quad w_{K}(t) = \frac{1}{1 + \sum_{j=1}^{K-1} \exp(\eta_j(t))}
\end{align*}
Note that the derivatives of $w_k$ are required in the calculation of the covariance matrices in Equation (\ref{eq:latentClassCov}). These can be derived analytically using the chain rule and pre-computed derivatives of the basis functions.

The parameters of this extension to the model are therefore $\Theta = \left(m, (\phi_0^{(k)}, \ldots, \phi_P^{(k)})_{k=1}^{K-1}, (\theta_k)_{k=1}^K, \sigma^2\right)$ and the Trend Direction Index in Proposition \ref{prop:TDIposterior} is as before defined conditionally on these variables.


####

To estimate the extended model with a latent class we further need to impose prior distributions on the spline coefficients for the probability functions $w_k$. A well-known problem when using spline expansions is selecting the number of knots or equivalently the degrees of freedom of the basis. To mitigate this problem we suggest using a large number of knots and then regularize the splines through the prior distributions on their coefficients. It can be shown that if all spline coefficients are equal, then the resulting expansion is a constant functions. This motivates imposing a random walk prior on the coefficients of each $\eta_k$ in Equation (\ref{eq:splineExpansion}) in order to penalize local variability. For each $k = 1,\ldots, K - 1$ spline expansion we therefore use the following prior distribution
\begin{align*}
\phi_{0}^{(k)} \sim N(0, 1), \quad \phi_{1}^{(k)} \sim N(0, 1), \quad \{\phi_{p}^{(k)} \sim N(\phi_{p-1}^{(k)}, \tau^{(k)})\}_{p=2}^P, \quad \tau^{(k)} \sim \text{Half-Normal}(0, 1)
\end{align*}



#####

We also fitted a model with $K = 2$ latent classes to the data. We used a cubic spline basis with $15$ degrees of freedom for the functional class probabilities. The left panel of Figure \ref{fig:smoking2-1} shows the posterior median of the probability that the latent class is equal to class $1$ at any given time. The right panel shows kernel density estimates of the posterior distributions of the length-scale parameter for each latent class. The plot suggests that some degree of non-stationarity is present in the latent dynamics. The latent function seems to develop at two length-scales, one with a posterior mode of $2.78$ and another with a posterior mode of $3.79$. The smaller length-scale has a higher influence during the first ten years with a bump around 2004 corresponding to the close to constant region seen in the upper left panel of Figure \ref{fig:smoking1}. Figure \ref{fig:smoking2-2} shows results similar to Figure \ref{fig:smoking1} but for the latent class model. The conclusions from this model fit are similar to the ones from before.


\begin{figure}[htb]
%\center\includegraphics{finalFitStan2}
\caption{}
\label{fig:smoking2-2}
\end{figure}




# Brokkassen
\begin{proposition}
Let the data generating model be defined as in Equation (\ref{eq:generatingProcess}) and $\mathbf{t}^\ast$ any finite vector of time points. Then the posterior distribution of $(f, df)$ at $\mathbf{t}^\ast$ is given by
\begin{align*}
\begin{bmatrix}f(\mathbf{t}^\ast)\\ df(\mathbf{t}^\ast) \end{bmatrix} \mid \mathbf{Y}, \mathbf{t}, \Theta &\sim
N\left(\begin{bmatrix}\mu_f(\mathbf{t}^\ast)\\ \mu_{df}(\mathbf{t}^\ast)\end{bmatrix}, \begin{bmatrix}\Sigma_{f,f}(\mathbf{t}^\ast, \mathbf{t}^\ast) & \Sigma_{f,df}(\mathbf{t}^\ast, \mathbf{t}^\ast)\\ \Sigma_{f,df}(\mathbf{t}^\ast, \mathbf{t}^\ast)^T & \Sigma_{df,df}(\mathbf{t}^\ast, \mathbf{t}^\ast)\end{bmatrix}\right)\\
\mu_f(\mathbf{t}^\ast) &= m + C_\theta(\mathbf{t}^\ast, \mathbf{t})\left[C_\theta(\mathbf{t}, \mathbf{t}) + \sigma^2 I\right]^{-1}\left(\mathbf{Y} - m\right)\nonumber\\
\mu_{df}(\mathbf{t}^\ast) &= \partial_1 C_\theta(\mathbf{t}^\ast, \mathbf{t})\left[C_\theta(\mathbf{t}, \mathbf{t}) + \sigma^2 I\right]^{-1}\mathbf{Y}\nonumber\\
\Sigma_{f,f}(\mathbf{t}^\ast, \mathbf{t}^\ast) &= C_\theta(\mathbf{t}^\ast, \mathbf{t}^\ast) - C_\theta(\mathbf{t}^\ast, \mathbf{t})\left[C_\theta(\mathbf{t}, \mathbf{t}) + \sigma^2 I\right]^{-1} C_\theta(\mathbf{t}^\ast, \mathbf{t})^T\nonumber\\
\Sigma_{df,df}(\mathbf{t}^\ast, \mathbf{t}^\ast) &= \partial_1 \partial_2 C_\theta(\mathbf{t}^\ast, \mathbf{t}^\ast) - \partial_1 C_\theta(\mathbf{t}^\ast, \mathbf{t})\left[C_\theta(\mathbf{t}, \mathbf{t}) + \sigma^2 I\right]^{-1} \partial_1 C_\theta(\mathbf{t}^\ast, \mathbf{t})^T\nonumber\\
\Sigma_{f,df}(\mathbf{t}^\ast, \mathbf{t}^\ast) &= \partial_1 C_\theta(\mathbf{t}^\ast, \mathbf{t}^\ast) - C_\theta(\mathbf{t}^\ast, \mathbf{t})\left[C_\theta(\mathbf{t}, \mathbf{t}) + \sigma^2 I\right]^{-1}\partial_1 C_\theta(\mathbf{t}^\ast, \mathbf{t})^T\nonumber
\end{align*}
\end{proposition}


[\textbf{TODO}: Since the outcomes are proportions calculated from known sample size we know the marginal variance of the outcomes by the formula for the variance of a sample proportion. In the stationary case with a standardized covariance function such that $C_\theta(t, t) = \alpha^2$ we obtain according to Equation (\ref{eq:generatingProcess}) the following equality constraint between the variance of the latent function and the residual variance
  \begin{align*}
  \Var[Y_i, \mid t_i, \Theta] = \alpha^2 + \sigma^2 = \frac{Y_i(1-Y_i)}{n_i}
  \end{align*}
  where $n_i$ is the number of responders at time $t_i$. How do we best enforce this constraint? \textbf{TODO:} Maybe all of this is an unnecceary complication that we don't need to mention.]




#####
\begin{align*}
  \bm{\mu}(\mathbf{t}^\ast) &= \begin{bmatrix}\mu_\beta(\mathbf{t}^\ast)\\ d\mu_\beta(\mathbf{t}^\ast)\\ d^2\!\mu_\beta(\mathbf{t}^\ast)\end{bmatrix} + \begin{bmatrix}C_\theta(\mathbf{t}^\ast, \mathbf{t})\\ \partial_1 C_\theta(\mathbf{t}^\ast, \mathbf{t})\\ \partial_1^2 C_\theta(\mathbf{t}^\ast, \mathbf{t})\end{bmatrix}K_{\theta,\sigma}(\mathbf{t}, \mathbf{t})^{-1} \left(\mathbf{Y} - \mu_\beta(\mathbf{t})\right)\\
  \bm{\Sigma}(\mathbf{t}^\ast,\mathbf{t}^\ast) &= \begin{bmatrix}C_\theta(\mathbf{t}^\ast, \mathbf{t}^\ast) & \partial_2 C_\theta(\mathbf{t}^\ast, \mathbf{t}^\ast) & \partial_2^2 C_\theta(\mathbf{t}^\ast, \mathbf{t}^\ast)\\ \partial_1 C_\theta(\mathbf{t}^\ast, \mathbf{t}^\ast) & \partial_1 \partial_2 C_\theta(\mathbf{t}^\ast,\mathbf{t}^\ast) & \partial_1 \partial_2^2 C_\theta(\mathbf{t}^\ast, \mathbf{t}^\ast)\\ \partial_1^2 C_\theta(\mathbf{t}^\ast, \mathbf{t}^\ast) & \partial_1^2 \partial_2 C_\theta(\mathbf{t}^\ast, \mathbf{t}^\ast) & \partial_1^2 \partial_2^2 C_\theta(\mathbf{t}^\ast,\mathbf{t}^\ast)\end{bmatrix} - \begin{bmatrix}C_\theta(\mathbf{t}^\ast, \mathbf{t})\\ \partial_1 C_\theta(\mathbf{t}^\ast, \mathbf{t})\\ \partial_1^2 C_\theta(\mathbf{t}^\ast, \mathbf{t})\end{bmatrix}K_{\theta,\sigma}(\mathbf{t}, \mathbf{t})^{-1}\begin{bmatrix}C_\theta(\mathbf{t}, \mathbf{t}^\ast)\\ \partial_2 C_\theta(\mathbf{t}, \mathbf{t}^\ast)\\ \partial_2^2 C_\theta(\mathbf{t}, \mathbf{t}^\ast)\end{bmatrix}
\end{align*}
and $K_{\theta,\sigma}(\mathbf{t}, \mathbf{t}) = C_\theta(\mathbf{t}, \mathbf{t}) + \sigma^2 I$.
