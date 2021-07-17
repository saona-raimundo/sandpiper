# mixture

## Objective

Compute expected polymorphisms.

## Background

### Polymorphisms
With random allele frequency x, it is given by $\mathbb{E}[2 x (1 - x)]$ where $f(x | s, h)$ is proportional to $\exp(2Ns(x^2 + 2hx(1-x))) * x^{(4NU-1)} (1-x)^{(4NU-1)}$. 

### Modelling choices

- s follows a SkewNormal distribution (mu, sigma, alpha)
- h equals $1 / (1 + e^{-\beta s})$  

### Parameters

- alpha
- beta
- mu
- sigma

## Algorithm

Consider beta fixed. Then, given the values $(\mathbb{E}[2 x (1 - x) \mid s])_{s \in \mathbb{R}}$, we have that
$$
\mathbb{E}[2 x (1 - x)] = \int_{\mathbb{R}} \mathbb{E}[2 x (1 - x) \mid S] f_S(s) ds \,.
$$
Therefore, one can make a bottom-up construction.

### Example

Consider $\beta = 0, \alpha = 0$. Then, an approximation is to apply the heat equation in the plane with boundary conditions $(\mathbb{E}[2 x (1 - x) \mid s])_{s \in \mathbb{R}}$ on the line.

## Interpretation

- $\beta$ and the cutting of max frequency
  - Determines initial conditions
    - $\beta$ interpolates as in the poster
    - Cutting of max frequency elevates positive selection
- $\alpha = 0$
  - Makes heat equation
- $\sigma$
  - Time
- $\mu$
  - Space  