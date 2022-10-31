#!/usr/bin/env python
# coding: utf-8

# # The Navier-Stokes Equations in Cylinderical Coordinates

# Define a cylindrical coordinate system $(x,r,\phi)$. In this system, a velocity vector is defined as:
#    
# \begin{equation}
# \overrightarrow{V} = 
# u\hat{e}_x + v\hat{e}_r + w\hat{e}_\phi
# \end{equation}
# 
# where $\hat{e}_x$, $\hat{e}_r$, and $\hat{e}_\phi$ are unit vectors.
# 
# The *del* operator $\overrightarrow{\nabla}$ is given by:
# 
# \begin{equation}
# \overrightarrow{\nabla} = 
# \hat{e}_x\frac{\partial}{\partial x} +
# \hat{e}_r\frac{\partial}{\partial r} +
# \frac{\hat{e}_\phi}{r}\frac{\partial}{\partial \phi}
# \end{equation}

# ## Derivation of Unit Vector Derivatives
# 
# A position vector $\overrightarrow{r}$ in cartesian coordinates is given by
# 
# \begin{equation}
# \overrightarrow{r} = x\hat{i} + y\hat{j} + z\hat{k}
# \end{equation}
# 
# The transformation from cartesian to cylindrical coordinates is given by:
# 
# \begin{align}
# x &= x \\
# y &= r\cos\phi \\
# z &= r\sin\phi
# \end{align}
# 
# Which gives this:
# 
# \begin{equation}
# \overrightarrow{r} = x\hat{i} + r\cos\phi\hat{j} + r\sin\phi\hat{k}
# \end{equation}
# 
# We can now form the derivatives, noting that $r$ is not dependent on either $x$ or $\phi$:

# \begin{equation}
# \hat{e}_x = \frac{\partial \overrightarrow{r}}{\partial x} = \hat{i}
# \end{equation}

# \begin{equation}
# \hat{e}_r = \frac{\partial \overrightarrow{r}}{\partial r} =
# \cos\phi\hat{j} + \sin\phi\hat{k}
# \end{equation}

# \begin{equation}
# \hat{e}_\phi = \frac{\partial \overrightarrow{r}}{\partial \phi} =
# -\sin\phi\hat{j} + \cos\phi\hat{k}
# \end{equation}

# From these equations, the derivatives of the unit vectors become:
# 
# \begin{equation}
# \frac{\partial\hat{e}_x}{\partial x} = \frac{\partial\hat{e}_x}{\partial r} = \frac{\partial\hat{e}_x}{\partial\phi} = 0
# \end{equation}

# \begin{equation}
# \frac{\partial\hat{e}_r}{\partial x} = 0
# \end{equation}

# \begin{equation}
# \frac{\partial\hat{e}_r}{\partial \phi} = - \sin\phi\hat{j} + \cos\phi\hat{k} = \hat{e}_\phi
# \end{equation}

# \begin{equation}
# \frac{\partial \hat{e}_\phi}{\partial x} = 
# \frac{\partial \hat{e}_\phi}{\partial r} = 0
# \end{equation}

# \begin{equation}
# \frac{\partial\hat{e}_\phi}{\partial\phi} =
# -\cos\phi\hat{j} - sin\phi\hat{k} = -\hat{e}_r
# \end{equation}

# We can ow write the derivatives of a general vector function $\overrightarrow{F}$ as follows:
# 
# \begin{equation}
# \overrightarrow{F} = F_1\hat{e}_x + F_2\hat{e}_r + F_3\hat{e}_\phi
# \end{equation}
# 
# \begin{equation}
# \frac{\partial\overrightarrow F}{\partial x} =
# \frac{\partial F_1}{\partial x}\hat{e}_x +
# \frac{\partial F_2}{\partial x}\hat{e}_r +
# \frac{\partial F_3}{\partial x}\hat{e}_\phi
# \end{equation}
# 
# \begin{equation}
# \frac{\partial\overrightarrow F}{\partial r} =
# \frac{\partial F_1}{\partial r}\hat{e}_x +
# \frac{\partial F_2}{\partial r}\hat{e}_r +
# \frac{\partial F_3}{\partial r}\hat{e}_\phi
# \end{equation}
# 
# \begin{equation}
# \frac{\partial\overrightarrow F}{\partial \phi} =
# \frac{\partial F_1}{\partial \phi}\hat{e}_x +
# \left(\frac{\partial F_2}{\partial \phi} - F_3\right)\hat{e}_r +
# \left(\frac{\partial F_3}{\partial \phi} - F_2\right)\hat{e}_\phi
# \end{equation}

# Using these equations we can come up with our governing equations:

# ## Continuity Equation
# 
# \begin{equation}
# \frac{\partial\rho}{\partial t} +
# \frac{\partial\rho u}{\partial x} +
# \frac{1}{r}\frac{\partial \rho v r}{\partial r} +
# \frac{1}{r}\frac{\partial\rho w}{\partial\phi} = 0
# \end{equation}
# 
# In writing the rest of the equations, we use the continuity equation to create the *conservative* form of the final equations.

# ## X Momentum Equation
# 
# \begin{equation}
# \frac{\partial \rho u}{\partial t} +
# \frac{\partial \rho u^2}{\partial x} +
# \frac{\partial \rho uv}{\partial r} +
# \frac{\rho w}{r}\frac{\partial u}{\partial\phi} = 
# -\frac{\partial p}{\partial x} +
# \frac{1}{Re}\left\{
# \frac{\partial\sigma_{xx}}{\partial x} +
# \frac{\partial\tau_{xr}}{\partial r} +
# \right\}
# \end{equation}

# ## r Momentum Equation
# 
# \begin{equation}
# \rho\frac{\partial v}{\partial t} +
# \rho u\frac{\partial v}{\partial x} +
# \rho v\frac{\partial v}{\partial r} +
# \frac{\rho w}{r}\left(
# \frac{\partial v}{\partial\phi} - w
# \right) =
# -\frac{1}{r}\frac{\partial p}{\partial r} +
# \frac{1}{Re}\left\{
# \frac{\partial\tau_{xr}}{\partial x} +
# \frac{\partial\sigma_{rr}}{\partial r} +
# \frac{1}{r}\left(\frac{\partial\tau_{r\phi}}{\partial \phi} - \sigma_{\phi\phi} + \sigma_{rr}\right)
# \right\}
# \end{equation}

# ## $\phi$ Momentum E
# quation
# 
# \begin{equation}
# \rho\frac{\partial w}{\partial t} +
# \rho u\frac{\partial w}{\partial x} +
# \rho v\frac{\partial w}{\partial r} +
# \frac{\rho w}{r}
# \left(
# \frac{\partial w}{\partial\phi} + v
# \right) =
# -\frac{1}{r}\frac{\partial p}{\partial \phi} +
# \frac{1}{Re}\left\{
# \frac{\partial\tau_{x\phi}}{\partial x} +
# \frac{\partial\tau_{r\phi}}{\partial r} +
# \frac{1}{r}\left(\frac{\partial\sigma_{\phi\phi}}{\partial \phi} + 2\tau_{r\phi}\right)
# \right\} =
# \end{equation}

# ## Energy Equation
# 
# \begin{equation}
# \rho\frac{\partial}{\partial t}
# \left(H - \frac{p}{\rho}\right) +
# \rho u\frac{\partial}{\partial x}
# \left(H - \frac{p}{\rho}\right) +
# \rho v\frac{\partial}{\partial r}
# \left(H - \frac{p}{\rho}\right) +
# \frac{\rho w}{r}\frac{\partial}{\partial \phi}
# \left(H - \frac{p}{\rho}\right) =
# \frac{1}{PrRe}
# \left\{
# \frac{\partial}{\partial x}
# \left(
# \mu \frac{\partial T}{\partial x}
# \right) +
# \frac{1}{r}\frac{\partial}{\partial r}
# \frac{\partial}{\partial x}
# \left(
# \mu\frac{\partial T}{\partial r}
# \right) +
# \frac{1}{r^2}\frac{\partial}{\partial \phi}
# \left(
# \mu\frac{\partial T}{\partial \phi}
# \right)
# \right\} -
# \frac{\partial}{\partial x}(up) -
# \frac{1}{r}\frac{\partial rvp}{\partial r} -
# \frac{1}{r}\frac{\partial wp}{\partial\phi} +
# \frac{1}{Re}
# \left\{
# \frac{\partial}{\partial x}\left(
# u\sigma_{xx} + v\tau_{xr} + w\tau_{x\phi}
# \right) +
# \frac{1}{r}
# \frac{\partial}{\partial r}\left(
# ru\tau_{xr} + rv\sigma_{rr} + rw\tau_{r\phi}
# \right) +
# \frac{1}{r}
# \frac{\partial}{\partial \phi}\left(
# u\tau_{x\phi} + v\tau_{r\phi} + w\sigma_{\phi\phi}
# \right)
# \right\}
# \end{equation}

# ## State Equation
# 
# \begin{equation}
# p = \frac{\gamma - 1}{\gamma}\rho t
# \end{equation}

# ## Vector Form
# 
# \begin{equation}
# {\bf E}_x + {\bf F}_r + {\bf G}_\phi = {\bf R}
# \end{equation}
# 
# 
# where:

# \begin{equation}
# \frac{\partial{(\bf E_I - E_V})}{\partial x} +
# \frac{\partial{(\bf F_I - F_V})}{\partial r} +
# \frac{1}{r}\frac{\partial{(\bf G_I - G_V})}{\partial y} + \frac{1}{r}R
# = 0
# \end{equation}

# \begin{equation}
# {\bf E_I} =
# \begin{bmatrix}
# \rho u r\\
# r \left(\rho u^2 + p\right) \\
# \rho u v r \\
# \rho u w r \\
# \rho u H r
# \end{bmatrix}
# \end{equation}

# \begin{equation}
# {\bf E_V} =
# \begin{bmatrix}
# 0 \\
# \sigma_{xx} \\
# \tau_{xr} \\
# \tau_{\theta z} \\
# u\sigma_{xx} + v\tau_{xr} + w\tau_{\theta x} + - q_x
# \end{bmatrix}
# \end{equation}

# \begin{equation}
# {\bf F_I} =
# \begin{bmatrix}
# \rho v \\
# \rho u v \\
# \rho v^2 + p \\
# \rho v w \\
# \rho v H
# \end{bmatrix}
# \end{equation}

# \begin{equation}
# {\bf F_V} =
# \begin{bmatrix}
# 0 \\
# \tau_{xr} \\
# \sigma_{rr} \\
# \tau_{\theta r} \\
# u\tau_{xr} + v\sigma_{rr} + w\tau_{\theta r} - q_r
# \end{bmatrix}
# \end{equation}

# \begin{equation}
# {\bf G_I} =
# \begin{bmatrix}
# \rho w \\
# \rho u w \\
# \rho v w \\
# \rho w^2 + p \\
# \rho w H
# \end{bmatrix}
# \end{equation}

# \begin{equation}
# {\bf G_V} =
# \begin{bmatrix}
# 0 \\
# \tau_{x\theta} \\
# \tau_{r\theta} \\
# \sigma_{\theta\theta} \\
# u\tau_{x\theta} + v\tau_{r\theta} + w\sigma_{\theta\theta} - q_z
# \end{bmatrix}
# \end{equation}

# \begin{equation}
# {\bf R} =
# \begin{bmatrix}
# \rho v \\
# \rho u v - \tau_{xr} \\
# \rho w^2 - \sigma_{rr} - \sigma_{\theta\theta} \\
# 2\rho v w - 2 \tau_{r\theta} \\
# \rho v H + q_r - u\tau_{xr} -v\sigma_{rr} - w\tau_{r\theta}
# \end{bmatrix}
# \end{equation}

# \begin{equation}
# E_t = 
# \frac{T}{\gamma(\gamma - 1)M^2} +
# \frac{1}{2}( u^2 + v^2 + w^2)
# \end{equation}

# \begin{equation}
# H = E + \frac{p}{\rho}
# \end{equation}

# \begin{equation}
# \sigma_{xx} =
# \frac{2\mu}{3Re} \left[
# 2\frac{\partial u}{\partial x} -
# \frac{\partial v}{\partial r} -
# \frac{1}{r}(
# \frac{\partial w}{\partial \theta} + v
# \right)
# \end{equation}

# \begin{equation}
# \sigma_{yy} =
# \frac{2\mu}{3Re_L} \left(
# 2\frac{\partial v}{\partial x} -
# \frac{\partial u}{\partial y} -
# \frac{\partial w}{\partial z}
# \right)
# \end{equation}

# \begin{equation}
# \sigma_{zz} =
# \frac{2\mu}{3Re_L} \left(
# 2\frac{\partial w}{\partial z} -
# \frac{\partial u}{\partial x} -
# \frac{\partial v}{\partial y}
# \right)
# \end{equation}

# \begin{equation}
# \tau_{xy} =
# \frac{\mu}{Re_L} \left(
# \frac{\partial u}{\partial y} +
# \frac{\partial v}{\partial x}
# \right)
# \end{equation}

# \begin{equation}
# \tau_{xz} =
# \frac{\mu}{Re_L} \left(
# \frac{\partial u}{\partial z} +
# \frac{\partial w}{\partial x}
# \right)
# \end{equation}

# \begin{equation}
# \tau_{yz} =
# \frac{\mu}{Re_L} \left(
# \frac{\partial v}{\partial z} +
# \frac{\partial w}{\partial y}
# \right)
# \end{equation}

# \begin{equation}
# q_x = \frac{\mu}{(\gamma-1)M_\infty Re_L Pr}
# \frac{\partial T}{\partial x}
# \end{equation}

# \begin{equation}
# q_y = \frac{\mu}{(\gamma-1)M_\infty Re_L Pr}
# \frac{\partial T}{\partial y}
# \end{equation}

# \begin{equation}
# q_z = \frac{\mu}{(\gamma-1)M_\infty Re_L Pr}
# \frac{\partial T}{\partial z}
# \end{equation}

# In[ ]:




