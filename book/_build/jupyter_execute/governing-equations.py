#!/usr/bin/env python
# coding: utf-8

# # Governing Equations

# ## General Form

# \begin{equation}
# \frac{\partial U}{\partial t} +
#     \frac{\partial A}{\partial z} +
#     \frac{\partial B}{\partial z} +
#     \frac{1}{r}\frac{\partial C}{\partial\theta} +
#     \frac{1}{r}D = 0
# \end{equation}

# \begin{equation}
# A = \begin{Bmatrix}
# \rho u\\
# \rho u u + p - \tau_{zz} \\
# \rho u v -\tau_{rz} \\
# \rho u w - \tau_{\theta z} \\
# \rho u H + q_z -u\tau_{zz} - v\tau_{rz} - w\tau_{\theta z}
# \end{Bmatrix}
# \end{equation}

# \begin{equation}
# B = \begin{Bmatrix}
# \rho v \\
# \rho u v - \tau_{rz} \\
# \rho v v + p - \tau_{rr} \\
# \rho v w - \tau_{\theta r}\\
# \rho v H + q_r - u\tau_{rz} -v\tau_{rr} - w\tau_{\theta r} 
# \end{Bmatrix}
# \end{equation}

# \begin{equation}
# C = \begin{Bmatrix}
# \rho w \\
# \rho u w  - \tau_{\theta z} \\
# \rho v w - \tau_{\theta r}  \\
# \rho w w + p - \tau_{\theta\theta} \\
# \rho w H + q_\theta - u\tau_{\theta z} - v\theta_{\theta r} - w\theta_{\theta\theta} 
# \end{Bmatrix}
# \end{equation}

# \begin{equation}
# D = \begin{Bmatrix}
# \rho v \\
# \rho u v - \tau_{rz} \\
# \rho v v - \rho w w - \tau_{rr} + \tau_{\theta\theta} \\
# 2\rho v w - 2\tau_{\theta r} \\
# \rho v H + q_r - u\tau_{rz} - v\tau_{rr} - w\tau_{\theta r}
# \end{Bmatrix}
# \end{equation}

# where:
# 
# \begin{align}
# E &= \frac{T}{\gamma(\gamma - 1)M^2} + \frac{1}{2}(u^2 + v^2 + w^2) \\
# H &= E + p/\rho \\
# \gamma &= 1.4
# \end{align}

# ## Shear Stress Terms
# 
# \begin{align}
# \tau_{zz} &= \frac{2\mu}{3Re}\left[
#     2\frac{\partial u}{\partial z} 
#     - \frac{\partial v}{\partial r} 
#     - \frac{1}{r}\Big(\frac{\partial w}{\partial\theta} + v \Big)
# \right] \\
# \tau_{rr} &= \frac{2\mu}{3Re}\left[
#     -\frac{\partial u}{\partial z} 
#     + 2\frac{\partial v}{\partial r} 
#     - \frac{1}{r}\Big(\frac{\partial w}{\partial\theta} + v \Big)   
# \right] \\
# \tau_{\theta\theta} &= \frac{2\mu}{3Re}\left[
#     -\frac{\partial u}{\partial z} 
#     - \frac{\partial v}{\partial r} 
#     + 2\frac{1}{r}\Big(\frac{\partial w}{\partial\theta} + v \Big)
# \right] \\
# \tau_{rz} &= \frac{\mu}{Re}\left[
#     \frac{\partial u}{\partial r} 
#     + \frac{\partial v}{\partial z} 
# \right] \\
# \tau_{\theta z} &= \frac{\mu}{Re}\left[
# \frac{\partial w}{\partial z} + \frac{1}{r}\frac{\partial u}{\partial\theta}
# \right] \\
# \tau_{\theta r} &= \frac{\mu}{Re}\left[
# \frac{1}{r}\Big(\frac{\partial v}{\partial\theta} - w\Big) + \frac{\partial w}{\partial r}
# \right] \\
# \end{align}

# \begin{align}
# q_z &= \frac{-\mu}{Pr(\gamma - 1)M^2Re}\frac{\partial T}{\partial z} \\
# q_z &= \frac{-\mu}{Pr(\gamma - 1)M^2Re}\frac{\partial T}{\partial r} \\
# q_z &= \frac{-\mu}{Pr(\gamma - 1)M^2Re}\frac{\partial T}{\partial\theta}
# \end{align}
# 
# where:
# \begin{equation}
# Pr = 0.72
# \end{equation}

# \begin{equation}
# p = \frac{\rho T}{\gamma M^2}
# \end{equation}

# The unknowns in this equation set are:
#     
# \begin{align}
# \rho &- density
# u &- axial velocity
# v &- radial velocity
# w &- azmuthal velocity
# p &- static pressure
# T &- static temperature
# H &- total enthalpy
# E &- total energy
# \mu &- viscosity
# M &- Mach Number
# \end{align}
# 
# That is ten equations. We have shown only seven equations. We need two more:
# 
# ### Sutherland's Law
# 
# \begin{equation}
# \mu = \mu(T)
# \end{equation}
# 
# ### Mach Number
# 
# \begin{align}
# c = \sqrt{\frac{\gamma p}{\rho}}
# \end{align}
# 
# \begin{equation}
# V = \sqrt{u^2 + v^2 + w^2}
# \end{equation}
# 
# \begin{equation}
# M = \frac{V}{c}
# \end{equation}
# 

# In[ ]:




