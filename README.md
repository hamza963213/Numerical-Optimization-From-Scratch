# Benchmark of numerical optimization algorithms
This repository provides a comparative analysis of various iterative descent optimization methods (Gradient Descent, Conjugate Gradient, Newton, Quasi-Newton, Gauss-Newton, Levenberg-Marquardt). These methods are implemented from scratch. 

> **Academic Context:** This work was done during Lab 1 (TP1) of the **CSOPT** (Scientific Computing and Optimization) course at **Centrale Nantes** (2024-2025).

*Note: The MATLAB code is commented in French. A Python translation with English comments will be added.*


### General Optimization Principle
All methods implemented in this project rely on the principle of iterative update:

$$x_{k+1} = x_k + \alpha_k d_k$$

Where:
* **$x_k$** is the current position at iteration $k$.
* **$d_k$** is the descent direction (determined by the gradient, the Hessian, etc., depending on the chosen method).
* **$\alpha_k$** is the step size (fixed, or variable using a line search technique such as the Armijo condition).



These algorithm performance are evaluated thgrough the minimiatiom of the Rosenbrock function (with parameters $n=2$ and $b=2$):
$$f(x_1, x_2) = 2(x_2 - x_1^2)^2 + (1 - x_1)^2$$

* **Global mimimum :** $(x_1, x_2) = (1, 1)$ avec $f(1, 1) = 0$.

> **Visualization**
> ![Plot of f0](chemin/vers/fonction_rosenbrock.png) 
> ![Level curves](chemin/vers/lignes_niveaux_rosenbrock.png)

---

The convergence of the algorithms depend heavily on the initial point $x_0$. We consider two cases.


#### Case 1 : Initial point $x_0 = (-4, 10)$

For this initial point, **second-order methods successfully converge** to the global minimum (1, 1), while **first-order methods fail** to reach it within the allowed iterations (the gradient stagnates at $\nabla f \neq 0$). While all second-order methods converge, their computational efficiency varies.

**Why?** The Rosenbrock function represents a flat, curved valley.
* **First-order methods** only use the slope, causing them to zig-zag and stagnate in the valley.
* **Second-order methods** use curvature (Hessian matrix) to adapt to the valley's geometry.

| Method | Iterations | Computation Time (s) | Final Point | Termination |
| :--- | :---: | :---: | :---: | :---: |
| Gradient Descent | 1000 | 1.3337 | (-3.50, 12.57) | Maxiter |
| Conjugate Gradient | 1000 | 1.7750 | (-2.48, 6.55) | Maxiter |
| Newton (unit step) | 6 | 1.4611 | (1.00, 1.00) | TolG |
| Newton (variable step) | 286 | 1.3117 | (1.00, 1.00) | TolG |
| Quasi-Newton (BFGS) | 310 | 1.6219 | (1.00, 1.00) | TolG |

**Specific Methods (Non-Linear Least Squares):**

| Method | Iterations | Computation Time (s) | Final Point | Termination |
| :--- | :---: | :---: | :---: | :---: |
| Gauss-Newton | 103 | 1.0817 | (1.00, 1.00) | TolF |
| Levenberg-Marquardt | 999 | 3.4641 | (1.00, 1.00) | TolF |

> **Evolution for $x_0 = (-4, 10)$**
> ![Gradient](chemin/vers/g_p2_-410.png) ![Criterion](chemin/vers/f_p2_-410.png)

#### Case 2 : Initial point $x_0 = (15, 3)$

For this initial point, only the Newton method with a unit step converge to the global minimum. The others either diverge or terminate prematurely.


| Method | Iterations | Computation Time (s) | Final Point | Termination |
| :--- | :---: | :---: | :---: | :---: |
| Gradient Descent | 1000 | 1.7149 | (-44.21, 94.14) | Maxiter |
| Conjugate Gradient | 1000 | 1.5176 | (426.69, -10.76) | Maxiter |
| Newton (unit step) | 6 | 0.9864 | (1.00, 1.00) | TolG |
| Newton (variable step) | 1000 | 1.4584 | (-1.17e3, 1.38e6) | Maxiter |
| Quasi-Newton | 136 | 1.7434 | (-2.39e3, 5.75e6) | TolF |

> **Evolution for $x_0 = (15, 3)$**
> ![Gradient](chemin/vers/g_p2_315.png) ![Criterion](chemin/vers/f_p2_315.png)
