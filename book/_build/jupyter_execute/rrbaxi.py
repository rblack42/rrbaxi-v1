#!/usr/bin/env python
# coding: utf-8

# # RRBAXI Decoded

# ## Coordinate Transformation

# \begin{align}
# (x,r,\theta) &= (\xi, \eta, \phi) \\
# \xi &= x \\
# \eta &= \frac{r(x) - r_b(x)}{r_s(x) - r_b(x)} \\
# \phi &= \theta
# \end{align}

# ## Partial Derivatives
# 
# The partial derivatives of the coordinate transformations are shown next:
# 
# \begin{equation}
# \xi_x = 1 \\
# \phi_\theta = 1
# \end{equation}
# 
# The equations for $\eta$ are a bit more complicated:

# \begin{equation}
# \eta = \frac{f(x)}{g(x)}
# \end{equation}
# 
# \begin{equation}
# \eta_x = \frac{g f_x - f g_x}{g^2}
# \end{equation}
# 
# \begin{equation}
# \eta_x = \frac{1}{g}(f_x - \eta g_x)
# \end{equation}
# 
# \begin{align}
# f &= r(x) - r_b(x) \\
# f_x &= -\frac{\partial r_b}{\partial x} \\
# g &= r_s - r_b \\
# g_x &= \frac{\partial r_s}{\partial x} - \frac{\partial r_b}{\partial x}
# \end{align}
# 
# therefore

# \begin{equation}
# \eta_x = \frac{1}{g}\biggl\{
# -\frac{\partial r_b}{\partial x} - 
# \eta \left(\frac{\partial r_s}{\partial x} - \frac{\partial r_b}{\partial x}\biggr\}
# \right)
# \end{equation}
# 
# \begin{equation}
# \eta_x = \frac{1}{g}\biggl\{
# (\eta - 1)\frac{\partial r_b}{\partial x} - 
# \eta \left(\frac{\partial r_s}{\partial x}\biggr\}
# \right)
# \end{equation}

# Let's try this with **sympy**:

# In[1]:


import sympy
x = sympy.symbols('x')
r = sympy.Function('r')(x)
r_s = sympy.Function('r_s')(x)
r_b = sympy.Function('r_b')(x)


# In[2]:


eta = (r - r_b)/(r_s - r_b)
eq2 = sympy.diff(eta,x)
eq2.simplify()


# \begin{equation}
# \eta_x = 
# \frac{\eta}{r_s - r_b}
# \left(\frac{\partial r_b}{\partial x} - \frac{\partial r_s}{\partial x}\right) -
# \frac{\frac{\partial r_b}{\partial x}}{r_s - r_b}
# \end{equation}

# \begin{equation}
# \eta_x = 
# \frac{1}{r_s - r_b}
# \left((\eta - 1)\frac{\partial r_b}{\partial x} - \eta\frac{\partial r_s}{\partial x}
# \right)
# \end{equation}
# 
# \begin{equation}
# \eta_r = \frac{1}{r_s - r_b}
# \end{equation}

# ## Governing Equation

# The governing equations can be written in conservative form:
# 
# \begin{equation}
# \frac{\partial U}{\partial t} +
# \frac{\partial E}{\partial x} + 
# \frac{\partial F}{\partial r} +
# \frac{1}{r}R = 0
# \end{equation}

# \begin{equation}
# \frac{\partial U}{\partial t} +
# \frac{\partial E}{\partial \xi} + \eta_x\frac{\partial E}{\partial\eta} +
# \eta_r \frac{\partial F}{\partial \eta} + \frac{R}{r} = 0
# \end{equation}

# \begin{equation} 
# U = \begin{bmatrix}
# \rho \\
# \rho u \\
# \rho v
# \end{bmatrix}
# \end{equation}

# \begin{equation}
# \nabla\cdot V = \eta_x u_\eta +\eta_r v_\eta + \frac{v}{r}
# \end{equation}

# \begin{equation}
# \sigma_{xx} = \frac{2\mu}{Re}\eta_x u_\eta -\frac{2\mu}{3 Re}\nabla\cdot V
# \end{equation}
# 
# \begin{equation}
# \tau_{xr} = \frac{\mu}{Re}\eta_x u_\eta -\frac{2\mu}{3 Re}\nabla\cdot V
# \end{equation} 
# 
# \begin{equation}
# \sigma_{rr} = \frac{2\mu}{Re}\eta_x v_\eta -\frac{2\mu}{3Re}\nabla\cdot V
# \end{equation}
# 
# \begin{equation}
# \sigma_{\phi\phi} = -p +\frac{2\mu}{Re}\frac{v}{r} -\frac{2\mu}{3Re}\nabla\cdot V
# \end{equation}

# \begin{equation}
# E = \begin{bmatrix}
# \rho u \\
# \rho u^2 + p - \sigma_{xx} \\
# \rho u v - \tau_{xr}
# \end{bmatrix}
# \end{equation}
# 
# \begin{equation}
# F = \begin{bmatrix}
# \rho v \\
# \rho u v - \tau_{xr} \\
# \rho v^2 + p -\sigma_{rr}
# \end{bmatrix}
# \end{equation}
# 
# \begin{equation}
# R = \begin{bmatrix}
# \rho v \\
# \rho u v - \tau_{xr} \\
# \rho vv - \sigma_{rr} - \sigma_{\phi\phi}
# \end{bmatrix}
# \end{equation}

# ## Parabolized Navier-Stokes

# For the parabolized Navier-Stokes equations, we will drop the time-dependent term, and set up a marching solution in the $\xi$ direction. We will drop the shear stress terms from the $E_\xi$ term.

# \begin{equation}
# E^* = \begin{bmatrix}
# \rho u \\
# \rho u^2 + p \\
# \rho u v
# \end{bmatrix}
# \end{equation}

# ### Solving $E^*$ for primative variables
# 
# In the PNS test code, we assume constant total enthalpy. 
# 
# \begin{equation}
# p = \rho R^* T
# \end{equation}
# 
# \begin{equation}
# c_p = \frac{\gamma}{\gamma - 1}R^*
# \end{equation}
# 
# \begin{equation}
# h = c_p T
# \end{equation}
# 
# Therefore:
# 
# \begin{equation}
# p = \frac{\gamma - 1}{\gamma} \rho h
# \end{equation}
# 
# \begin{equation}
# H^* = h + \frac{u^2+ v^2}{2}
# \end{equation}
# 
# \begin{equation}
# c = \sqrt{(\gamma - 1)h}
# \end{equation}

# In[3]:


At this point, we have eight equations in eight unknowns. 


# \begin{equation}
# H = \frac{1}{2} + \frac{1}{(\gamma - 1)M^2}
# \end{equation}
# 
# \begin{equation}
# p = \frac{1}{\gamma M^2}
# \end{equation}

# Let:
# 
# \begin{equation}
# a = \rho u \\
# b = \rho u^2 + p \\
# c = \rho u v
# \end{equation}

# Let's let **sympy** solve this:

# In[11]:


c_p,R,gamma = sympy.symbols('c_p R gamma')


# In[13]:


eq1 = sympy.Eq(c_p,gamma/(gamma-1)*R)
eq1


# In[17]:


a, M_x, u, v, p, rho, T, H, h = \
    sympy.symbols('a M_x u v p rho T H h')


# In[20]:


eq2 = sympy.Eq(h + 0.5*(u**2 + v**2),H)
eq2


# In[23]:


eq3 = sympy.Eq(sympy.sqrt(gamma*R*T),a)
eq3


# In[25]:


eq4 = sympy.Eq(u,a * M_x)
eq4


# In[26]:


u.subs 


# In[ ]:




