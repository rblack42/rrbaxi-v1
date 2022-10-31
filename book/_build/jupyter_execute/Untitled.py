#!/usr/bin/env python
# coding: utf-8

# # Axisymmetric Navier Stokes Equations
# 
# The general form of the axisymmetric Navier-Stokes equations is given as follows:
# 
# \begin{equation}
# \frac{\partial U}{\partial t} +
# \frac{\partial (E_i - E_v)}{\partial x} +
# \frac{\partial (F_i - F_v)}{\partial r} + (H_i - H_v) = 0
# \end{equation}

# The vectors in the above equation are:
#     
# \begin{equation}
# U = \begin{bmatrix}
# \rho \\
# \rho u \\
# \rho v \\
# \rho e
# \end{bmatrix}
# \end{equation}
# 
# \begin{equation}
# E_i = \begin{bmatrix}
# \rho u \\
# \rho u^2 + p \\
# \rho uv \\
# (\rho e + p)u
# \end{bmatrix}
# \end{equation}
# 
# \begin{equation}
# G_i = \begin{bmatrix}
# \rho v \\
# \rho uv \\
# \rho v^2 + p \\
# (\rho e + p)v
# \end{bmatrix}
# \end{equation}
# 
# \begin{equation}
# H_i = \frac{1}{r}\begin{bmatrix}
# \rho v \\
# \rho uv \\
# \rho v^2 \\
# (\rho e + p)v
# \end{bmatrix}
# \end{equation}

# \begin{equation}
# E_v = \begin{bmatrix}
# 0 \\
# \lambda(\nabla\cdot q)u + 2\mu u(\xi_x u_\xi + \eta_x u_\eta) 
# \mu v(\xi_x v_\xi + \eta_x v_\eta + \xi_r u_\xi + \eta_x u_\eta) \\
# \lambda(\nabla\cdot q) + 2\mu(\xi_x u_\xi + \eta_x u_\eta) \\
# \mu(\xi_x v_\xi + \eta_x v_\eta + \xi_r u_\xi + \eta_r u_\eta) \\
# \lambda(\nabla\cdot q) + 2\mu u(\xi_x u_\xi + \eta_x u_\eta) +
# \mu v(\xi_x v_\xi + \eta_x v_\eta + \xi_r u_\xi + \eta_r u_\eta) +
# k(\xi_x T_\xi + \eta_x T_\eta)
# \end{bmatrix}
# \end{equation}

# \begin{equation}
# F_v = \begin{bmatrix}
# 0 \\
# \mu(\xi_r u_\xi + \eta_r u_\eta +\xi_x v_\xi + \eta_x v_\eta) 
# \lambda(\nabla\cdot q) + 2\mu(\xi_r v_\xi + \eta_r v_\eta) \\
# \mu(\xi_x v_\xi + \eta_x v_\eta + \xi_r u_\xi + \eta_r u_\eta) \\
# \lambda(\nabla\cdot q)v + \mu u(\xi_r u_\xi + \eta_r u_\eta + \xi_xv_\xi + \eta_x v_\eta) +
# 2\mu v(\xi_r v_\xi + \eta_r v_\eta) +
# k(\xi_r T_\xi + \eta_r T_\eta)
# \end{bmatrix}
# \end{equation}

# \begin{equation}
# h_v = \frac{1}{r}\begin{bmatrix}
# 0 \\
# \mu(\xi_r u_\xi + \eta_r u_\eta + \xi_x v_\xi + \eta_x v_\eta) \\
# 2\mu(\xi_r v_\xi + \eta_r v_\eta) = 2\mu\frac{v}{r} \\
# \lambda(\nabla\cdot q)v + \mu u(\xi_r u_\xi + \eta_r u_\eta + \xi_x v_\xi + \eta_x v_\eta) +
# 2\mu v(\xi_r v_\xi + \eta_r v_\eta) + k(\xi_r T_\xi + \eta_r T_\eta)
# \end{bmatrix}
# \end{equation}

# \begin{equation}
# \nabla\cdot q = \xi_x u_\xi + \eta_x u_\eta + \xi_r v_\xi + \eta_r v_\eta + \frac{v}{r}
# \end{equation}

# For this study, $\xi$ = $x$ and $\eta = \frac{r - r_b}{r_s - r_b}$
# 
# Therefore $\xi_x = 1$ and $\xi_r = 0$ and $\eta_x = 0$ and $\eta_r = 1/(r_s - r_b)$

# In[ ]:




