#!/usr/bin/env python
# coding: utf-8

# # Solving For Primative Variables

# The governing equations for parabolized Navier-Stokes are set to generate these variables as the solution proceeds:

# \begin{align}
# A &= \rho u\\
# B &= \rho u^2 + p \\
# C &= \rho u v \\
# D &= \rho u w \\
# E &= \rho u H
# \end{align}

# We need to set up a routine that will reduce these variables to the basic flow properties. To do this, we need some additional equations. We will use **sympy** to set up these equations.
# 
# A simple inspection will show that we can eliminate three of these variables immediately:
# 
# \begin{align}
# v &= \frac{C}{A} \\
# w &= \frac{D}{A} \\
# H &= \frac{E}{A}
# \end{align}
# 
# This leaves us with two equations in three remaining unknowns $\rho, u, p$.

# In[1]:


from sympy import var, solve, simplify, expand, collect


# When using **sympy**, we need to create equations where the right-hand side is always zero. That is simple enough is we just move any terms on the right over to the left side.

# In[2]:


A, B = var('A B')
rho, u, p = var('rho u p')


# Now, we create our equations:

# In[3]:


eq1 = rho*u - A
eq2 = rho*u**2 + p -B


# ## State Equation
# 
# From our previous work, here is the state equation in  nondimensional form:
#     
# \begin{equation}
# p = \frac{\gamma - 1}{\gamma} \rho T
# \end{equation}   

# Using this equation adds two more variables, one of which is a constant for our work: $\gamma = 1.4$. We need to add these terms to *SymPy*, then add the state equation to our working set:

# In[4]:


T, gamma = var('T gamma')
eq3 = (gamma/(gamma-1))*rho*T - p


# ## Definition of Total Enthalpy

# Next, we will use the total enthalpy to complete our equation set:
# 
# \begin{equation}
# H = h + \frac{1}{2}\left(u^2 + v^2 + w^2\right)
# \end{equation}
# 
# Define a new variable:
# 
# \begin{equation}
# K = H - \frac{1}{2}(v^2 + w^2)
# \end{equation}
# 
# Thus:
# 
# \begin{equation}
# K = h + \frac{u^2}{2}
# \end{equation}

# In[5]:


h, K = var('h K')
eq4 = h + (u**2)/2 - K


# ### Definition of Specific Heat Coefficient
# 
# \begin{equation}
# c_p = \frac{\gamma}{\gamma-1}R_{gas}
# \end{equation}

# In[6]:


c_p, R = var('c_p R')
eq5 = (gamma/(gamma-1))*R - c_p


# ### Definition of Static Enthalpy
# 
# \begin{equation}
# h = c_p T
# \end{equation}
# 
# or
# 
# \begin{equation}
# \rho c_p T = \frac{\rho R_{gas}}{\gamma - 1}
# \end{equation}

# In[7]:


eq6 = c_p * T - h


# In[8]:


eq7 = solve([eq1,eq2,eq3,eq4],[rho,u,p,T])


# In[9]:


rho_c = (eq7[0][0])
rho_c


# In[10]:


u_c = eq7[0][1]
u_c


# In[11]:


p_c = eq7[0][2]
p_c


# In[12]:


T_c = eq7[0][3]
T_c


# That gives us six equations in the six unknowns we need to evaluate:
#     
# \begin{equation}
# \rho, u, p, T, c_p, h
# \end{equation}

# *SymPy* could solve this set of equations, but the result would be pretty messy. WE can simplify the set somewhat with a little more work.

# Looking at the equation set, we can eliminate the simple variables $C_p, v, w, H$ since whey fall out directly. 
# 
# \begin{align}
# C_p &= \frac{\gamma}{(\gamma-1)}R_{gas} \\
# v &= C/A \\
# w &= D/A \\
# H &= E/A \\
# \end{align}
# 
# Furthermore, $h$ is related to $T$ by a simple equation. 
# 
# \begin{equation}
# h = c_p T
# \end{equation}
# 
# The definition of total enthalpy can be reduced a bit as well. We define a new variable, $K$ as:
# 
# \begin{equation}
# K = \frac{1}{2}\left(v^2 + w^2\right) - H
# \end{equation}
# 
# \begin{equation}
# K = \frac{1}{2}\left(
# \frac{C^2}{A^2} + \frac{D^2}{A^2}
# \right) - \frac{E}{A} 
# \end{equation}
# 
# Let's recreate our equation set with these reductions:

# In[13]:


R, K = var('R K')
Eq1 = rho*u - A
Eq2 = rho*u**2 + p - B
Eq3 = rho*R*T - p
Eq4 = c_p*T + (u**2)/2 + K


# In[14]:


res = solve(
        [Eq1,Eq2,Eq3,Eq4],
        [rho,u,p,T]
)
res


# Well, we still have some messy equations to work with. Let's try to simplify these. We will start by looking at the first results for $each of our primitives from the solution:

# In[15]:


rho_c = res[0][0] # rho
rho_c


# In[16]:


u_c = res[0][1] # u
u_c


# In[17]:


p_c = res[0][2] # p
p_c


# In[18]:


T_c = res[0][3] # T
T_c


# At this point, we could solve for all of our primitive variables using the above equations. However, we can simplify things by introducing the axial *Mach Number*:
# 
# \begin{equation}
# M_x = \frac{u}{\sqrt{(\gamma-1)T}}
# \end{equation}
# 
# Therefore:
# 
# \begin{equation}
# u^2 = (\gamma - 1)T M_x^2
# \end{equation}

# Using this equation, we get this equation for $B$:
#     
# \begin{equation}
# B = (\gamma - 1)\rho T M_x^2  + p
# \end{equation}

# From the state equation, we find:
# 
# \begin{equation}
# \rho T = \frac{\gamma p}{\gamma - 1}
# \end{equation}

# Substituting this, we now get:
#     
# \begin{equation}
# B = p ( 1 + \gamma M_x^2)
# \end{equation}
# 
# or:
# \begin{equation}
# p = \frac{B}{1 + \gamma M_x^2}
# \end{equation}

# From the definition of total enthalpy:
# 
# \begin{equation}
# H = c_p T + \frac{1}{2}\left((\gamma - 1)T M_x^2 + v^2 + w^2)\right) 
# \end{equation}
# 
# from the definition of $K$ above above:
# 
# \begin{equation}
# H = K - \frac{1}{2}\left(v^2 + w^2\right)
# \end{equation}
# 
# or:
# 
# \begin{equation}
# \left(v^2 + w^2\right) = 2(H - K)
# \end{equation}

# Plugging this back into our equation for $H$ above, we get this:
# 
# \begin{equation}
# c_p T + \frac{1}{2}\left((\gamma - 1)T M_x^2\right) - K = 0 
# \end{equation}
# 
# Rearranging, we get this equation for $T$:
# 
# \begin{equation}
# T = \frac{K}{\biggl\{c_p + \left(\frac{(\gamma - 1)}{2}M_x^2\right)\biggr\}}
# \end{equation}
# 

# From the state equation, we get this equation:
# 
# \begin{equation}
# \rho = \frac{\gamma p}{(\gamma - 1) T} =
# \end{equation}

# In[19]:


Mx2 = var('Mx2')
p2 = B/(1 + gamma * Mx2)
T2 = K/(c_p + (gamma - 1)/2*Mx2)
rho2 = gamma*p2/((gamma - 1)*T2)
rho2


# Substituting for $p$ and $T$, we get this:
# 
# \begin{equation}
# \rho = \frac{\gamma B\left(C_p + \frac{\gamma-1}{2}M_x^2\right)}
# {K(\gamma - 1)(1 + \gamma M_x^2)}
# \end{equation}

# In[20]:


usq2 = (A/rho2)**2
usq2


# In[21]:


usq3 = ((gamma - 1)*T2*Mx2**2)
usq3


# In[22]:


f1 = usq2/usq3
f1.simplify()


# Introduce a new term:
# 
# \begin{equation}
# \alpha = \frac{2A^2K(\gamma-1)}{B^2\gamma^2}
# \end{equation}

# In[23]:


alpha = var('alpha')
f2 = f1*(B**2*gamma**2*alpha)/(2*A**2*K*(gamma - 1))
f2.simplify()


# This expression is equal to one, so we can write:
# 
# \begin{equation}
# \alpha(M_x^2\gamma + 1)^2 = M_x^2(M_x^2(\gamma - 1) + 2c_p)
# \end{equation}

# In[24]:


f3 = f2 - 1
f3.simplify()


# since this expression is equal to zero, we can reduce it:

# In[25]:


eq11 = f3*(Mx2**2*(Mx2*(gamma-1) + 2*c_p))
eq12 = eq11.expand()
eq12.simplify()


# Here is that equation in "normal" form:
#     
# \begin{equation}
# {M_x^2}^3\left(1-\gamma\right) +
# {M_x^2}^2\left(\alpha\gamma^2-2c_p\right) +
# {M_x^2}\left(2\alpha\gamma\right) +
# \left(\alpha\right) = 0
# \end{equation}

# In[26]:


aa = (1-gamma)
aa


# In[27]:


bb = alpha*gamma**2 - 2*(c_p)
bb


# In[28]:


cc = 2*alpha*gamma
cc


# In[29]:


dd = alpha
dd


# This is a simple polynomial in $M^2_x$. Let's see what *SymPy* can do with such an equation.

# In[30]:


x = var('x')
poly = aa*x**3 + bb*x**2 + cc*x + dd
poly


# only one of these results has a real value:

# In[31]:


M = solve([poly],[x])
M


# \begin{equation}
# {M_x^2}^2(\gamma B) + M_x^2(\gamma A) + A = 0
# \end{equation}

# ## unverified code below ===============================================

# In[ ]:


phi = var('phi')
Ksub = gamma*B**2/(2*(gamma-1)*A**2*K)
Ksub


# In[ ]:


test = T_c.subs(K,Ksub)
test.simplify()


# ## Definition of Mach Number
# 
# We will focus only on the axial Mach Number $M_x$:
# 
# \begin{equation}
# M_x = \frac{u}{a}
# \end{equation}

# In[ ]:


a, M_x = var('a M_x')
Eq5 = M_x * a - u


# ## Speed of Sound
# 
# \begin{equation}
# a = \sqrt{\gamma R T}
# \end{equation}

# This result in $u$:
# 
# \begin{equation}
#  u =  M_x \sqrt{\gamma RT}
# \end{equation}
# 
# With this equation, we can evaluate $\rho u^2$:
# 
# \begin{equation}
# \rho u^2 = M_x^2(\gamma \rho R T)
# \end{equation}
# 
# Substituting the state equation, we can get this result:
# 
# \begin{equation}
# B = p(1 + M_x^2\gamma)
# \end{equation}
# 
# or
# 
# \begin{equation}
# p = \frac{B}{1+ M_x^2\gamma}
# \end{equation}

# In[ ]:


eq6 = B/(1+gamma*M_x**2) - p
eq6


# From the definition of the total enthalpy we can get this:
# 
# \begin{equation}
# c_p T = H - \frac{1}{2}(u^2 + v^2 + w^2) 
# \end{equation}
# 
# Let K be defined as follows:
# 
# \begin{equation}
# K = H - \frac{1}{2}(v^2 + w^2)
# \end{equation}
# 
# Now we have this equation:
# 
# \begin{equation}
# \frac{\gamma}{\gamma - 1}R T = K -\frac{1}{2}M_X^2\gamma R T
# \end{equation}
# 
# Rearranging:
# 
# \begin{equation}
# K = \left(\frac{1}{\gamma - 1} +\frac{1}{2}M_x^2\right)\gamma R T
# \end{equation}
# 
# \begin{equation}
# K = \left(\frac{\gamma}{\gamma - 1} + \frac{\gamma M_x^2}{2}\right) R T
# \end{equation}

# \begin{equation}
# 
# \end{equation}

# For isentropic flow of a perfect gas we have:
# 
# \begin{equation}
# \frac{p}{\rho^\gamma} = constant
# \end{equation}

# \begin{align}
# A &= \rho u\\
# B &= \rho u^2 + p \\
# \end{align}
# 
# Therefore:
# 
# \begin{equation}
# \frac{A^2K}{B^2} =  \frac{K\rho^2 u^2}{\rho^2 u^4 + 2 \rho u^2 p + p^2}
# \end{equation}
# 
# 

# In[ ]:


A = rho * u
B = rho * u**2 + p
K= (1/(gamma - 1)+M_x**2/2)*gamma*Rgas*T
p = B/(1 + gamma*M_x**2)

phi = 2*(gamma -1)*A**2*K/(gamma*B**2)
phi.simplify()


# In[ ]:




