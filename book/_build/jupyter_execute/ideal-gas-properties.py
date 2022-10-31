#!/usr/bin/env python
# coding: utf-8

# # Ideal Gas Properties

# In this section, we will derive a few equations that describe the properties of an ideal gas. To assist in this work, we will use *Python* *SymPy* to validate the equations. Note that *SymPy* does not care if the equations look "neat" so I will need to jockey things to get the needed form. Remember that all equations in *SymPy* need to be of the form $exp = 0$. (The right hand side is assumed.)
# 
# To get started, we first set up the *SymPy* functions we will be using.

# In[1]:


from sympy import var, solve, simplify, expand, collect


# I wil be using [Principles of Ideal-Fluid Aerodynamics](karamcheta1966), as a reference for the equations for an ideal gas:

# ## State Equation
# 
# In dimensional form, the state equation is given by:
# 
# \begin{equation}                                                            \hat{p} = \hat{\rho} \hat{R} \hat{T}
# \end{equation}
# 
# The hat symbol indicates a dimensional property.
# 
# It is convenient to nondimensionalize the governing equations using reference property values, usually those found in the free-stream.

# ### Nondimensionalizing Variables
# 
# We nondimensionalize all variables in this study using free-stream reference values:
# 
# \begin{align}
# p &= \frac{\hat{p}}{\hat{p}_{ref}} \\
# \rho &= \frac{\hat{\rho}}{\hat{\rho}_\infty} \\
# u &= \frac{\hat{u}}{\hat{u}_\infty} \\
# T &= \frac{\hat{T}}{\hat{T}_{ref}} \\
# H &= \frac{\hat{H}}{\hat{H}_{ref}}
# \end{align}

# where
# 
# \begin{align}
# \hat{H}_{ref} &= \hat{U}^2_\infty \\
# \hat{T}_{ref} &= \frac{\hat{U}^2_\infty}{\hat{C_p}} \\
# \hat{p}_{ref} &= \hat{\rho}_\infty\hat{U}^2_\infty
# \end{align}

# Now we will let *SymPy* generate the nondimensional form of this equation:

# In[2]:


rho, rho_ref = var('rho rho_ref')
T, T_ref = var('T T_ref')
Cp_ref = var('Cp_ref')


# In[3]:


p, p_ref = var('p p_ref')
R = var('R_ref')
u, u_ref = var('u u_ref')


# Now, we create expressions for the dimensional properties:

# In[4]:


u_h = u*u_ref
T_h = T*u_ref**2/Cp_ref 
p_h = p*(rho_ref*u_ref**2)
rho_h = rho*rho_ref


# Now, let's use the state equation:

# In[5]:


state = p_h - rho_h*R*T_h
eq1 = state/(rho_ref*u_ref**2)
eq1.simplify()


# From ideal gas theory:
# 
# \begin{equation}
# c^*_p - c^*_v = R^*
# \end{equation}
# 

# 
# Dividing by $c^*_p$ gives
# 
# \begin{align}
# \frac{R^*}{c^*_p} &= 1 - \frac{c^*_v}{c^*_p}\\
# &= 1 - \frac{1}{\gamma}
# &= \frac{\gamma - 1}{\gamma}
# \end{align}

# In[6]:


gamma = var('gamma')
eq2 = eq1.subs(R_ref,Cp_ref*gamma/(gamma-1))
eq2.simplify()


# In[7]:


nstate = solve([eq2],[p])
nstate


# Getting *SymPy* to display the final result proved difficult, so I will just give it here:
# 
# ### Nondimensional State Equation
# 
# \begin{equation}
# p = \frac{\gamma}{\gamma - 1}\rho T
# \end{equation}

# Using these expressions we get these nondimensional variables:
# 
# \begin{align}
# p &= \frac{P^*}{\rho^*_{ref} {U^{*^2}}_\infty} \\
# \rho &= \frac{\rho^*}{\rho^*_\infty} \\
# T &= \frac{c_p^* T^*}{{U^{*^2}_\infty}} \\
# H &= \frac{H^*}{{U^{*^2}}_\infty}
# \end{align}}
# 
# ### Non-Dimensional State Equation
# 
# \begin{equation}
# \rho^*_\infty {{U^*}^2}_\infty(p) =
# \rho^*_\infty(\rho)R^*\frac{{{U^*}^2}_\infty}{c^*_p}(T)
# \end{equation}
# 
# \begin{equation}
# \rho^*_\infty {{U^*}^2}_\infty(p) =
# \frac{\rho^*_\infty R^*{{U^*}^2}_\infty}{c^*_p}(\rho T)
# \end{equation}
# 
# \begin{equation}
# p = \rho T \frac{R^*_\infty}{c^*_p} 
# \end{equation}

# Which gives this form for the state equation:
# 
# \begin{equation}
# p = \frac{\gamma - 1}{\gamma} \rho T
# \end{equation}

# ## Speed of Sound
# 
# \begin{equation}
# a = \sqrt{\frac{\gamma p}{\rho}}
# \end{equation}
# 
# but:
# 
# \begin{equation}
# \frac{p}{\rho} = \frac{\gamma - 1}{\gamma}T
# \end{equation}
# 
# \begin{equation}
# a = \sqrt{(\gamma-1)T}
# \end{equation}

# ## Mach Number
# 
# \begin{equation}
# M_x = \frac{u}{a}
# \end{equation}
# 
# \begin{equation}
# u^2 = (\gamma - 1)T M_x^2
# \end{equation}
# 
# From the definition of $B$ we get this:
# 
# \begin{align}
# B &= p + \rho u^2 \\
#  &= p + (\gamma - 1)\rho T M_x^2 \\
#  &= p + \gamma p M_x^2 \\
#  &= p(1 + \gamma M_x^2)
# \end{align}Therefore:
# 
# \begin{equation}
# p = \frac{B}{1 + \gamma M_x^2}
# \end{equation}

# We now have enough equations to look at how we will solve the governing equations in nondimensional form.
# 

# From the definition of total enthalpy we get this:
# 
# \begin{equation}
# T =  H -\frac{1}{2}\left( (\gamma - 1)T M_x^2 + v^2 + w^2 \right) 
# \end{equation}
# 
# Now define a new term $K$ as follows:
# 
# \begin{equation}
# K = H - \frac{1}{2}\left(v^2 + w^2\right)
# \end{equation}
# 
# Then:
# 
# \begin{equation}
# T = K -\frac{1}{2}(\gamma - 1)T M_x^2
# \end{equation}
# 
# or
# 
# \begin{equation}
# T\left(1 + \frac{(\gamma - 1)M_x^2}{2}\right) = K 
# \end{equation}
# 
# Rearranging:
# 
# \begin{equation}
# T = \frac{K}{\left(1 + \frac{(\gamma - 1)M_x^2}{2}\right)}
# \end{equation}

# From the state equation, we have:
# 
# \begin{equation}
# \rho = \frac{\gamma p}{(\gamma - 1)T}
# \end{equation}
# 
# Substituting for $p$ and $T$, we get this equation:
# 
# \begin{equation}
# \rho = \frac{\gamma B\left(1 + \frac{\gamma - 1}{2} M_x^2\right)}{K(\gamma - 1)(1 +\gamma M_x^2)}
# \end{equation}

# From the definition of $A$ we get:
#     
# \begin{equation}
# u = \frac{A}{\rho} = \frac{AK(\gamma - 1)(1 +\gamma M_x^2)}{\gamma B\left(1 + \frac{\gamma - 1}{2} M_x^2\right)}
# \end{equation}
# 
# Let's simplify this a bit by grouping terms:
# 
# \begin{equation}
# \alpha = \frac{AK(\gamma - 1)}{\gamma B}
# \end{equation}
# 
# With this, our reduced equation for $u$ becomes:
# 
# \begin{equation}
# u = \frac{\alpha(1 +\gamma M_x^2)}{\left(1 + \frac{\gamma - 1}{2} M_x^2\right)}
# \end{equation}

# Now, from our definition of the Mach Number:
# 
# \begin{equation}
# u^2 = (\gamma - 1)T M_x^2
# \end{equation}
# 
# Eliminating $T$, we get this equation:
# 
# \begin{equation}
# u^2 = \frac{(\gamma - 1)K M_x^2}{\left(1 + \frac{(\gamma - 1)}{2}M_x^2\right)}
# \end{equation}

# Combining these two equations for $u^2$, we get this:
# 
# \begin{equation}
# \frac{(\gamma - 1)K M_x^2}{\left(1 + \frac{(\gamma - 1)}{2}M_x^2\right)} =
# \Biggl\{
# \frac{\alpha(1 +\gamma M_x^2)}{\left(1 + \frac{\gamma - 1}{2} M_x^2\right)}
# \Biggr\}^2
# \end{equation}

# Rearranging:
# 
# \begin{equation}
# \alpha^2(1 + \gamma M_x^2)^2 = \biggl(1 + \frac{\gamma - 1}{2} M_x^2\bigg)(\gamma - 1)K M_x^2
# \end{equation}

# \begin{equation}
# \frac{\alpha^2}{(\gamma - 1)K}(1 + \gamma M_x^2)^2
# \end{equation}

# Let's let **sypy** reduce this to an equation for $M_x^2$:

# In[8]:


gamma, alpha, K, A, B, MX2 = var('gamma alpha K A B MX2')
eq10 = alpha**2*(1 + gamma*MX2)**2 - (1 + (gamma-1)/2*MX2)*(gamma-1)*K*MX2
eq11 = collect(eq10, MX2)
eq12 = eq11.expand()
eq13 = collect(eq12,MX2)
eq13


# In[9]:


A = eq13.coeff(MX2,2)
collect(A,K)


# Collecting terms of $M_x^2$ we get:
# 
# \begin{equation}
# (M_x^2)^2(K\gamma + \alpha^2\gamma^2 - \frac{\gamma^2 K}{2}) + 
# (M_x^2)(\gamma K + K + 2\alpha^2\gamma  +
# \alpha^2
# \end{equation}

# In[10]:


eqt1 = -((gamma - 1)**2)/2
eqt1.expand()


# In[11]:


B = eq13.coeff(MX2,1)
collect(B,K)


# In[12]:


res = solve([eq13], [MX2])
res[0][0].simplify()


# In[ ]:




