#!/usr/bin/env python
# coding: utf-8

# # Fortran Solver

# 
# \begin{equation}
# H = H_\infty-\frac{(C/A)^2}{2}
# \end{equation}
# 
# \begin{equation}
# \phi=\frac{2(\gamma-1)KA^2}{\gamma B^2}
# \end{equation}
# 
# \begin{equation}
# phm = \frac{\gamma}{\gamma + 1}
# \end{equation}
# 
# \begin{equation}
# rad=\sqrt{1-\phi-\phi/\gamma}
# \end{equation}
# 
# \begin{equation}
# den=\gamma\phi-(\gamma - 1)
# \end{equation}
# 
# \begin{equation}
# M_x = \frac{1.-\phi+rad}{den}
# \end{equation}
# 
# \begin{equation}
# p = \frac{B}{1+\gamma*M_x}
# \end{equation}
# 
# \begin{equation}
# t = \frac{K}{1.+\frac{\gamma - 1}{2}*M_x}
# \end{equation}
# 
# \begin{equation}
# \rho = \frac{\gamma*p}{\gamma - 1)t}
# \end{equation}
# 
# \begin{equation}
# u = A/\rho
# \end{equation}
# 
# 

# ## Parabolized Navier Stokes Solver

# \begin{equation}
# U = \begin{bmatrix}
# \rho u \\
# \rho u^2 + p \\
# \rho uv \\
# \rho uw \\
# (\rho e_t + p)u
# \end{bmatrix} = 
# \begin{bmatrix}
# A \\
# B \\
# C \\
# D \\
# E 
# \end{bmatrix}
# \end{equation}
# 
# where:
# 
# \begin{equation}
# e_t = e + \frac{u^2 + v^2 + w^2}{2}
# \end{equation}
# 
# \begin{equation}
# p = (\gamma - 1)\rho e
# \end{equation}
# 
# \begin{equation}
# T = \gamma M_\infty^2 \frac{p}{\rho}
# \end{equation}

# This gives:
# 
# Unknowns: $ \rho\ u\ v\ w\ p\ e$

# In[1]:


from sympy import var, solve, simplify


# In[2]:


rho, u, v, w, p, e, gamma = var('rho u v w p e gamma')
A, B, C, D, E, K = var('A B C D E K')


# In[3]:


eq1 = rho*u - A
eq2 = rho*u**2 + p - B
eq3 = rho*u*v - C
eq4 = rho*u*w - D
e_t = K + (u**2)/2
eq5 = (rho*e_t + p)*u - E
eq6 = p*(gamma - 1)*rho*e - p


# In[4]:


sol = solve([eq1,eq2,eq3,eq4,eq5,eq6],[rho, u, v, w, p, e])
sol


# In[5]:


rho = sol[0][0]
rho


# In[6]:


u = sol[0][1]
u


# In[7]:


v = sol[0][2]
v


# In[8]:


w = sol[0][3]
w


# In[9]:


p = sol[0][4]
p


# In[10]:


e = sol[0][5]
e


# In[ ]:




