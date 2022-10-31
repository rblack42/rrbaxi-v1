#!/usr/bin/env python
# coding: utf-8

# # Axisymmetric PNS for Ogive-Cylinder

# The governing Navier-Stokes equations in cylindrical coordinates are given below:
#     
# \begin{equation}
# {\bf E}_x + {\bf F}_r + \frac{1}{r}{\bf G}_\theta = \frac{1}{r}{\bf R}
# \end{equation}

# Where:

# \begin{equation}
# {\bf E} =
# \begin{Bmatrix}
# \rho u \\
# \rho u^2 + P - \frac{1}{Re}\sigma_{xx} \\
# \rho u v - \frac{1}{Re}\tau_{xr} \\
# \rho u w - \frac{1}{Re}\tau_{x\theta} \\
# \rho u H + \frac{1}{PrRe}\mu T_x -
#     \frac{1}{Re}(u\sigma_{xx} +v\tau_{xr} + w\tau_{x\theta}
# \end{Bmatrix}
# \end{equation}

# \begin{equation}
# {\bf F} =
# \begin{Bmatrix}
# \rho v \\
# \rho u v - \frac{1}{Re}\tau_{xr} \\
# \rho v^2 + P - \frac{1}{Re}\tau_{rr} \\
# \rho v w - \frac{1}{Re}\tau_{r\theta} \\
# \rho v H + \frac{1}{PrRe}\mu T_r -
#     \frac{1}{Re}(u\tau_{xr} +v\sigma_{rr} + w\tau_{r\theta})
# \end{Bmatrix}
# \end{equation}

# \begin{equation}
# {\bf G} =
# \begin{Bmatrix}
# \rho  w \\
# \rho u w - \frac{1}{Re}\tau_{x\theta} \\
# \rho v w - \frac{1}{Re}\tau_{r\theta}) \\
# \rho w^2 + P - \frac{1}{Re}\sigma_{\theta\theta} \\
# \rho w H + \frac{1}{PrRe}\mu T_\theta -
#     \frac{1}{Re}(u\tau_{x\theta} +v\tau_{r\theta} + w\sigma_{\theta\theta})
# \end{Bmatrix}
# \end{equation}

# \begin{equation}
# {\bf R} =
# \begin{Bmatrix}
# 0 \\
# 0 \\
# \rho w^2 + P - \frac{1}{Re}\sigma_{\theta\theta} \\
# \frac{1}{Re}\tau_{r\theta} - \rho v w \\
# 0
# \end{Bmatrix}
# \end{equation}

# The shear stress terms age given by:
# 
# \begin{equation}
# \sigma_{xx} = -P +
# \lambda\Bigl\{\frac{\partial u}{\partial x} +
# \frac{\partial v}{\partial r} +
# \frac{1}{r}\frac{\partial w}{\partial \theta} + \frac{v}{r}
# \Bigr\}  + 2\mu\frac{\partial u}{\partial x}
# \end{equation}

# \begin{equation}
# \sigma_{rr} = -P +
# \lambda \Big\{
# \frac{\partial u}{\partial x} +
# \frac{\partial v}{\partial r} +
# \frac{1}{r}
# \frac{\partial w}{\partial \theta} + \frac{v}{r} 
# \Bigr\} + 2\mu\frac{\partial v}{\partial r}
# \end{equation}

# \begin{equation}
# \sigma_{\theta\theta} = -P +
# \lambda \Big\{
# \frac{\partial u}{\partial x} +
# \frac{\partial v}{\partial r} +
# \frac{1}{r}
# \frac{\partial w}{\partial \theta} + \frac{v}{r} 
# \Bigr\} + \frac{2\mu}{r}\frac{\partial w}{\partial \theta}
# \end{equation}

# \begin{equation}
# \tau_{xr} = \mu\Bigl\{
# \frac{\partial u}{\partial x} + 
# \frac{\partial v}{\partial r}
# \Bigr\}
# \end{equation}

# \begin{equation}
# \tau_{xr} = \mu\Bigl\{
# \frac{\partial u}{\partial x} + 
# \frac{\partial v}{\partial r}
# \Bigr\}
# \end{equation}

# \begin{equation}
# \tau_{r\theta} = \mu\Bigl\{
# \frac{\partial w}{\partial r} + 
# \frac{1}{r}
# \frac{\partial v}{\partial \theta} -
# \frac{w}{r}
# \Bigr\}
# \end{equation}

# 
# ## Generalized Coordinate Transformation

# For the computational grid surrounding the ogive-cylinder body, we will create a coordinate system that extends from the body to a point outside of the shock cone:
# 
# \begin{equation}
# \{x,r,\theta\} \Leftrightarrow \{\xi,\eta,\zeta\}
# \end{equation}

# \begin{align}
# \xi &= \xi(x) \\
# \eta &= \eta(x,r,\theta) \\
# \zeta &= \zeta(x,r,\theta)
# \end{align}

# The partial derivatives become:
#     
# \begin{equation}
# \frac{\partial}{\partial x} = \xi_x\frac{\partial}{\partial\xi} +
# \eta_x\frac{\partial}{\partial \eta} +
# \zeta_x\frac{\partial}{\partial \zeta}
# \end{equation}

# \begin{equation}
# \frac{\partial}{\partial r} = \eta_r\frac{\partial}{\partial\eta} +
# \zeta_r\frac{\partial}{\partial \zeta}
# \end{equation}

# \begin{equation}
# \frac{\partial}{\partial \theta} = \eta_\theta\frac{\partial}{\partial\eta} +
# \zeta_\theta\frac{\partial}{\partial \zeta}
# \end{equation}

# \begin{equation}
# J = x_\xi\{ r_\eta \theta_\zeta - r_\zeta \theta_\eta \}
# \end{equation}

# \begin{equation}
# \xi_x{\bf E_p^*}_\xi +
# \eta_x{\bf E}_\eta^* +
# \zeta_x{\bf E}_\zeta^* +
# \eta_r{\bf F}_\eta^* +
# \zeta_r{\bf F}_\zeta^* +
# \eta_\theta{\bf G}_\eta^* +
# \zeta_\theta{\bf G}_\zeta^*
# = {\bf R}
# \end{equation}

# \begin{equation}
# {\bf E_p^*} =
# \begin{Bmatrix}
# (\rho u) r\\
# (\rho u^2 + P)r  \\
# (\rho u v)r \\
# (\rho u w)r  \\
# \rho u H)r 
# \end{Bmatrix}
# \end{equation}

# \begin{equation}
# {\bf E^*} =
# \begin{Bmatrix}
# (\rho u) r\\
# (\rho u^2 + P - \frac{1}{Re}\sigma^*_{xx})r \\
# (\rho u v)r - \frac{1}{Re}\tau^*_{xr})r \\
# (\rho u w)r - \frac{1}{Re}\tau^*_{xr})r \\
# (\rho u H - \frac{1}{PrRe}\mu T^*_x -
#     \frac{1}{Re}(u\sigma^*_{xx} +v\tau^*_{xr} + w\tau^*_{x\theta})r
# \end{Bmatrix}
# \end{equation}

# \begin{equation}
# {\bf F^*} =
# \begin{Bmatrix}
# (\rho v) r\\
# (\rho u v - \frac{1}{Re}\tau^*_{xr})r \\
# (\rho v^2 + P - \frac{1}{Re}\sigma^*_{rr})r \\
# (\rho v w)r - \frac{1}{Re}\tau^*_{r\theta})r \\
# (\rho v H - \frac{1}{PrRe}\mu T^*_r -
#     \frac{1}{Re}(u\tau^*_{xr} +v\sigma^*_{rr} + w\tau^*_{r\theta})r
# \end{Bmatrix}
# \end{equation}

# \begin{equation}
# {\bf G^*} =
# \begin{Bmatrix}
# (\rho w) r\\
# (\rho u w - \frac{1}{Re}\tau^*_{x\theta})r \\
# (\rho v w - \frac{1}{Re}\tau^*_{r\theta})r \\
# (\rho w^2 - \frac{1}{Re}\sigma^*_{\theta\theta})r \\
# (\rho w H - \frac{1}{PrRe}\mu T^*_\theta -
#     \frac{1}{Re}(u\tau^*_{x\theta} +v\tau^*_{r\theta} + w\sigma^*_{\theta\theta})r
# \end{Bmatrix}
# \end{equation}

# \begin{equation}
# {\bf R^*} =
# \begin{Bmatrix}
# 0 \\
# 0 \\
# \rho w^2 + P -\frac{1}{Re}\sigma_{\theta\theta}
# \frac{1}{Re}\tau_{r\theta} - \rho v w \\
# 0
# \end{Bmatrix}
# \end{equation}

# The transformed shear stress terms now become:

# \begin{equation}
# \sigma^*_{xx} = 
# \lambda (\overline{\nabla} \cdot \overline{\bf V})^*  +
# 2\mu\{\eta_x u_\eta + \zeta_x u_\zeta\}
# \end{equation}

# \begin{equation}
# \sigma^*_{rr} = 
# \lambda (\overline{\nabla} \cdot \overline{\bf V})^*  +
# 2\mu\{\eta_r v_\eta + \zeta_r v_\zeta\}
# \end{equation}

# \begin{equation}
# \sigma^*_{\theta\theta} = 
# \lambda (\overline{\nabla} \cdot \overline{\bf V})^*  +
# \frac{2\mu}{r}\{\eta_\theta w_\eta + \zeta_\theta w_\zeta + v\}
# \end{equation}

# \begin{equation}
# \tau^*_{xr} = 
# \mu\{\eta_r u_\eta + \zeta_r u_\zeta + \eta_x v_\eta + \zeta_x v_\zeta\}
# \end{equation}

# \begin{equation}
# \tau^*_{x\theta} = 
# \mu\{\eta_x w+\eta + \zeta_x w_\zeta + \frac{1}{r}(\eta_\theta u_\eta + \zeta_\theta u_\eta \}
# \end{equation}

# \begin{equation}
# \tau^*_{r\theta} = 
# \mu\{\eta_r w_\eta + \zeta_r w_\zeta + \frac{1}{r}(\eta_\theta v_\eta + \zeta_\theta v_\zeta - w)\}
# \end{equation}

# \begin{equation}
# \left(\overline{\nabla} \cdot \overline{\bf V}\right)^* =
# \{
# \eta_x u_\eta + 
# \zeta_x u_\zeta +
# \eta_r v_\eta + 
# \zeta_r v_\zeta +
# \frac{1}{r}(
# \eta_\theta w_\eta + 
# \zeta_\theta w_\zeta + v
# )
# \}
# \end{equation}

# ## Axisymmetric Equations

# For initial testing, we will conside the axisymmetric case for the ogive-cylinder at zero angle of attack. This reduces the problem from a 3D one to a 2D 1. In this case all derivatives in the $\theta$ ($\zeta$) direction are zero and the corresponding **w** velocity is zero as well. We This reduces the equation set to this:
# 
# \begin{equation}
# \xi_x{\bf E}_\xi + \eta_x{\bf E}_\eta +
# \eta_r{\bf F}_\eta
# = {\bf R}
# \end{equation}

# \begin{equation}
# {\bf E} =
# \begin{Bmatrix}
# \rho u \\
# \rho u^2 + P - \frac{1}{Re}\sigma_{xx} \\
# \rho u v - \frac{1}{Re}\tau_{xr} \\
# \rho u H + \frac{1}{PrRe}\mu T_x -
#     \frac{1}{Re}(u\sigma_{xx} +v\tau_{xr}
# \end{Bmatrix}
# \end{equation}

# \begin{equation}
# {\bf F} =
# \begin{Bmatrix}
# \rho v \\
# \rho u v - \frac{1}{Re}\tau_{xr} \\
# \rho v^2 + P - \frac{1}{Re}\tau_{rr} \\
# \rho v H + \frac{1}{PrRe}\mu T_r -
#     \frac{1}{Re}(u\tau_{xr} +v\sigma_{rr})
# \end{Bmatrix}
# \end{equation}

# \begin{equation}
# {\bf R} =
# \begin{Bmatrix}
# 0 \\
# 0 \\
# - \rho v w \\
# 0
# \end{Bmatrix}
# \end{equation}

# ## Solving For Primative Variables

# \begin{align}
# \rho u &=  A \\
# \rho u^2 + P &= B \\
# \rho u v &= C \\
# \rho u w &= D \\
# \rho u H &= E
# \end{align}
# 
# \begin{equation}
# M_x = \frac{u}{c} = \frac{u}{\sqrt{(\gamma - 1) T}}
# \end{equation}
# or
# \begin{equation}
# u^2 = (\gamma - 1) T M_x^2
# \end{equation}

# \begin{equation}
# \rho = \frac{\gamma P}{(\gamma - 1)T}
# \end{equation}

# \begin{equation}
# B = P + \rho u^2 = P +\rho\{(\gamma - 1)T M^2_x \} = P\{1 + \gamma M_x^2\}
# \end{equation}

# \begin{equation}
# T = H -\frac{1}{2}\bigl\{(\gamma-1)T
# \end{equation}

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




