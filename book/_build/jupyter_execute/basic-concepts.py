#!/usr/bin/env python
# coding: utf-8

# # Basic Concepts

# This study will involve analyzing the flow of air over a simple body of revolution. We will be describing that air by specifying some basic properties that we can measure. For this reason, it makes sense t start our study by looking at the*units of measure* we will be using.

# ## Units of Measure

# There are two measurement systems in common use today: the *English System* and the *Metric System*. (Also called the *British System* and the *Standard International System*.) The *units of measure* used in these two systems are shown below:
# 
# | Unit |English | Metric |
# |:----:|:------:|:------:|
# | Mass | slug | kg|
# | Length | ft | meter |
# | Time | second | second|
# | Force | lbf | Newton |
# 
# Force and mass are related by Newton's Second Law:
# 
# \begin{equation}
# force = mass * acceleration\ (f=ma)
# \end{equation}

# So, what is a "slug", and what is a "Newton"? Neither term is in common use, except by engineers and scientists. Let's consider how we measure the *weight* of an object.
# 
# According to Newton's Law, the *weight* of n object is a measure of the force that object exerts on a scale of some ind. The force is the result of the pull of gravity on all objects near the Earth. 
# 
# The acceleration due to gravity varies with the distance from the center of the Earth. For research purposes we use a [Standard Gravity](https://en.wikipedia.org/wiki/Standard_gravity):
# 
# | Unit | English | Metric |
# :-----:|:-------:|:------:|
# | acceleration of gravity | 32.17405 $\frac{ft}{sec^2}$ | 9.80665 $\frac{m}{sec^2}$ |

# Using Newton's Second Law and these values, we can figure out what a *slug* and a *Newton* are.
# 
# Suppose we have an object that weight 100 lbf in the English system. That gives us a mass of:
# 
# \begin{equation}
# mass = 100 \frac{100}{32.17405}\frac{lbf\ sec^2}{ft}
# \end{equation}

# A mass of one $\frac{lbf\ sec^2}{ft}$ is called a $slug$:

# In[1]:


english_mass = 100/32.17405
english_mass


# The mass is 3.11028 $\frac{lbf\ sec^2}{ft}$ or 3.11028 $slug$

# The conversion factor from $lbf$ to $Newton$ is 4.44822. Therefore, that same mass in the metric system is:

# In[2]:


metric_mass = mass * 4.44822
metric_mass


# So the mass of 3.11028 $slug$ is identical to a mass of 13.82549 $Newtons$.
# 
# All of these basic measures are standardized by an internation body of scientists. 

# ## Dimensional Analysis

# When performing engineering calculations, it is critically important that the equations used be consistent in the units of measure they use. Part of making sure things work correctly involves someting called *dimensional analysis*. This is a simple process where we tract the basic units attached to every property we involve in a calculation to make sure we are not doing math on inconsistent values. 
# 
# In this analysis we use simple symbols to represent the units involved: $MLT\Theta$. These represent:
# 
# - M - Mass
# - L - Length
# - T - Time
# - $\Theta$ - Temperature
# 
# Let's look at how these symbols can be used.

# ### Weight

# Acceleration is measured in distance per second squared: $\frac{L}{T^2}$. Using *Newton's Second Law* we know that a force is the product of mass times the acceleration due to gravity. Therefore, the weight of an object will have units of $\frac{ML}{T^2}$.

# ### Density

# The *density* of something is the amount of that substance per unit of volume. Thus, the units of *density* are $\frac{M}{L^3}$

# ### Energy

# Another important property we study is *energy*. There are several different forms of *energy*: 
# 
# - Potential Energy - the "potential" for an object to have motion, such as lifting a weight.
# - Kinetic Energy - the energy derived from an object in motion. The basic units of *energy* are weight times length: $\frac{ML^2}{T^2}$

# ### Temperature

# Another common property is the *temperature* of an object. At the atomic level, *temperature* is a measure of the *kinetic energy* of the molecules of the object. The higher the *temperature* the faster the molecules move around. 
