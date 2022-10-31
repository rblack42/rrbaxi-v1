#!/usr/bin/env python
# coding: utf-8

# # Fluid Properties

# The various fluid properties we study have wildly different magnitudes and units of measure. Making sure all units correct in our calculations can be a challenge. Even worse, we can run into difficulties with numerical precision since computers are limited in the range of values they can handle.
# 
# The most common solution to these issues is to eliminate the units of measure, and create a set of *nondimensional* properties whose magnitudes should be more reasonable for our calculations. 
# 
# In this section, we will set up these nondimensional properties.

# In our study, we have these fluid properties:
#     
#     
# |Property|Symbol|Units|
# |:---------:|:--------:|:--------------:|
# | density  | $\rho$ | $\frac{slugs}{ft^3}$ |
# | pressure | $p$ | $\frac{lb_f}{ft^2}$ |
# | Temperature | $T$ | $^o R$ |
# | Velocity | $V$ | $\frac{ft}{sec}$ |
# | Thermal Conductivity | $k$ | $\frac{lb_f}{sec ^o R}$ |
# | Kinematic Viscosity | $\mu$ | $\frac{lb_f sec}{ft^2}$ |

# We can use the *Python* **pint** package to work with these units. Here is some code that will let us manage units:

# In[1]:


from pint import UnitRegistry


# The *UnitRegistry* is a class that will create an object that understands how units can be converted between systems of measure. Only one instance of this class can be created for any project. We will do that setup work next:

# In[2]:


u = UnitRegistry()


# We can demonstrate the use of this registry by defining the gas constant $R$ for air. From Wikipedia [Gas Constant](https://en.wikipedia.org/wiki/Gas_constant) we find that the *Specific Gas Constant* for air has a value of 1716.46 $\frac{ ft lb_f}{slugs ^oR}$. Here is how we express that value using **pint**:

# In[3]:


Rgas = 1716.46*u.ft*u.lbf/(u.slug*u.rankine)
Rgas


# Internally, **pint** uses "base" units, which are metric. We can generate the correct values for the properties in these base units using the **to_base_units()** method:

# In[4]:


Rgas.to_base_units()


# There are additional properties we will need. The *Specific Heat* properties are defined as follows:
# 
# \begin{equation}
# R_{gas} = c_p - c_v 
# \end{equation}
# 
# Where $c_p$ is the *specific heat at constant pressure* and $c_v$ is the *specific heat at constant volume. These two properties are commonly related as follows:
# 
# \begin{equation}
# \gamma = \frac{c_p}{c_v} = 1.4\ (for Air)
# \end{equation}

# From *Wikipedia [Specific Heat Capacity](https://en.wikipedia.org/wiki/Specific_heat_capacity#Imperial_engineering_units) we find this value: $c_p$ = 520.3 $\frac{joule}{kg ^oK}$.
# 
# This value is valid for temperatures at room temperature. Unfortunately, this value varies with temperature. 

# In[5]:


c_p = 520.3 * u.joule/(u.kg*u.kelvin)
c_p.to_base_units()


# Even though the units for $c_p$ did not look the same as those used for $R_{gas}, **pint** showed they are the same (only expressed in different measuring systems).

# ## Example Test Data

# In my original study, I used test data from wind-tunnel tests conducted at the *USAF Arnold Engineering Development Center* for flow over an ogive-cylinder. The test data included these values:
# 
# - Mach Number: 5.95
# - $T_0$: 830 $^oR$
# - $P_\infty$:  0.167 psi
# - $q_\infty$: 4.13 psi
# - $\frac{Re}{ft}$: 4.9e6

# In[6]:


p = 0.167 * u.lbf/u.inch**2
p


# ### Test Temperature

# We need to determine the working temperature for this test:
# 
# The stagnation temperature $T_0$ is defined as:
# 
# \begin{equation}
# \frac{T_0}{T} = 1 + \frac{\gamma - 1}{2} M^2
# \end{equation}
# 

# In[7]:


gamma = 1.4
Mach = 5.95
Tref = 830*u.rankine/(1 + (gamma - 1)/2 * Mach**2)
Tref


# To calculate the value for $c_p$ we need to do a bit of setup with **pint**:

# ### Equation of State

# All gases obey the *State Equation:
#     
# \begin{equation}
# p = \rho R_{gas} T
# \end{equation}
# 
# Since we have the working pressure and temperature, we can find the density:

# In[8]:


rho = p / (Rgas * Tref)
rho.to_base_units()


# ### Speed of Sound

# The speed of sound is defined as:
# 
# \begin{equation}
# c = \sqrt{\gamma\frac{p}{\rho}}
# \end{equation}

# In[9]:


c = (gamma*p/rho)**0.5
c.to_base_units()


# From the wind-tunnel *Mach Number*, we can now calculate the test velocity:

# In[10]:


Mach = 5.95
Vref = Mach * c
Vref.to_base_units()


# I wonder what that is in miles-per-hour!

# In[11]:


Vref.to('mph')


# The dynamic pressure $q_\infty$ is defined as:
# 
# \begin{equation}
# q = \frac{1}{2}\rho u^2
# \end{equation}

# From this eqution, we can calculate the velocity a different way as a check:

# In[12]:


Qref = 4.13* u.lbf/u.inch**2
Uref = (2*Qref/rho)**0.5
Uref.to_base_units()


# Those two calculations are pretty close.

# ### Specific Heat 

# In addition to these data values, the wind tunnel value for $c_p$ is given by this formula:
# 
# \begin{equation}
# c_p =  (3.15789e-5 T + 0.098947) \frac{btu}{lbm ^oR}
# \end{equation}

# In[13]:


btu = u.Quantity("1 BTU")
btu.to_base_units()


# In[14]:


lbm = 1*u.lbf/u.gravity
lbm.to_base_units()


# In[15]:


Cp =  (3.15789e-5 * Tref.magnitude + 0.098947)  * btu/(lbm*u.rankine)
Cp.to_base_units()


# ### Reynolds Number

# The equation defining the *Reynolds Number* is:
# 
# \begin{equation}
# Re = \frac{\rho U L_{ref}}{\mu}
# \end{equation}

# The test model is an ogive-cylinder, 50 inches long. According to the test data, that would produce a *Reynolds Number* of

# In[16]:


Re = 4.9e6 * 8.5/12
Re


# This value is actually higher than the value I found in the original study. That value was 
# $Re = 2179168.0$
# 
# TODO: Check this with AEDC

# ### Viscosity

# We can calculate the viscosity coefficient by using the definition of the *Reynolds Number*:
#     
# \begin{equation}
# Re = \frac{\rho u L}{\mu}
# \end{equation}

# In[17]:


L = 8.5/12*u.foot
mu = (rho * Vref * L)/Re
mu.to_base_units()


# From the original study, the reference viscosity was given as:

# In[18]:


Mu_ref = 7.65034e-7 * lbm/(u.ft*u.second)
Mu_ref.to_base_units()


# In[19]:


REc = 2179168.0
Muc = (rho * Vref * L)/REc
Muc.to_base_units()


# Again, the value calculated here does not match data used in the original study.
# 
# TOTO: Check with AEDC

# In[20]:


L_ref = REc*Mu_ref/(rho * Vref)
L_ref.to('inch')


# From that calculation, the reference length does not look right either. 
# 
# TOTO: Check with AEDC

# ## Dimensional Property Summary

# For reference, here is a summary of the values we calculated (subject to corrections based on new data later!)

# In[21]:


props = {}
props['density'] = rho
props['pressure'] = p
props['temperature'] = Tref
props['viscosity'] = mu
props['mach'] = Mach
props['RE'] = Re
props['Rgas'] = Rgas
props['c_p'] = c_p
props['gamma'] = gamma


# In[22]:


for key in props:
    try:
        print(key,'=',props[key].to_base_units())
    except:
        print(key,'=',props[key])


# In[ ]:




