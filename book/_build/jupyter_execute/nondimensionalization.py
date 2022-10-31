#!/usr/bin/env python
# coding: utf-8

# # Nondimensionalization

# In order to simplify our calculations and generate property values in a more reasonable number range, we will reduce out fluid properties to a nondimensional form.
# 
# This is a fairly simple process. It is common to use *free-stream* conditions as reference values. WE them divide our fluid properties by the associated free-stream values to produce a nondimensional form. In this process, we will indicate dimensional properties using a "hat" symbol over the property. Free-stream, properties will use an $\infty$ subscript:
# 
# Let's start out with some obvious properties:

# \begin{align}
# \rho &= \frac{\hat{\rho}}{\hat{\rho}_\infty} \\
# T &= \frac{\hat{T}}{\hat{T}_\infty} \\
# U &= \frac{\hat{U}}{\hat{U}_\infty} \\
# \mu = \frac{\hat{\mu}}{\hat{\mu}_\infty}
# \end{align}

# The governing equations involve derivatives of property combinations with respect to some coordinate axis. We will nondimensionalize lengths usinf a reference length. For this study, that length is related to the test model, an ogive-cylinder tested at AEDC in the late 1970s.
# 
# \begin{equation}
# l = \frac{\hat{l}}{\hat{l}_{ref}}
# \end{equation}
# 
# A typical derivitive looks like this:
# 
# \begin{equation}
# \frac{\partial}{\partial x} = \hat{l}_{ref}\frac{\partial}{\partial \hat{x}}
# \end{equation}

# ## Governing Equation Property Terms

# In reviewing the governing equations, we see various combinations of basic fluid properties. We need to consider those combinations to come up with suitable nondimensional forms.

# ### Continuity 

# From the *Continuity Equation, we see terms that look like this:
# 
# \begin{equation}
# \frac{\partial(\hat{\rho}\hat{U})}{\partial \hat{x}}
# \end{equation}
# 
# Using free-stream properties for reference values, we get this nondimensional form:
# 
# \begin{equation}
# \frac{\partial \rho U}{\partial x} = 
# \frac{\partial(\hat{\rho}\hat{U})}{\partial\hat{x}} * \frac{\hat{l}_{ref}}{\hat{\rho}_\infty\hat{U}_\infty}
# \end{equation}

# ### Momentum

# The momentum equations are more complex. 
# 
# \begin{equation}
# \frac{\partial(\rho U^2)}{\partial x} =
# \frac{\partial(\hat{\rho}\hat{U}^2)}{\partial \hat{x}} * \frac{\hat{l}_{ref}}{\hat{\rho}_\infty\hat{U}^2_\infty}
# \end{equation}

# The pressure term found in the *momentum equations* has these units: 
# 
# \begin{equation}
# \frac{\partial\hat{p}}{\partial \hat{x}} \hat{=}
# \frac{\hat{p}_{ref}}{\hat{l}_{ref}}
# \end{equation}

# For unit consistency, we should use a different value for the reference pressure:
# 
# \begin{equation}
# \hat{p}_{ref} \hat{=} \hat{\rho}_\infty \hat{U}^2_\infty
# \end{equation}
# 
# Therefore:
# 
# \begin{equation}
# \frac{\partial p}{\partial x} = 
# \frac{\partial\hat{p}}{\partial\hat{x}} * 
# \frac{\hat{l}_{ref}}{\hat{\rho}_\infty \hat{U}_\infty}
# \end{equation}

# The shear stress terms have this form:
# 
# \begin{equation}
# \hat{\mu}\frac{\partial^2 \hat{U}}{\partial \hat{x}^2} \hat{=}
# \frac{\hat{\mu}_{ref}\hat{U}_\infty}{\hat{l}
# _{ref}^2}
# \end{equation}

# For consistency, we need to create a new parameter $R$:
# 
# \begin{equation}
# R \hat{=} \frac{\hat{\mu}_{ref}\hat{U_\infty}}{\hat{l}^2_{ref}} * 
# \frac{\hat{l}_{ref}}{\hat{\rho}_\infty \hat{U}^2_\infty}
# \end{equation}

# Simplifying:
# 
# \begin{equation}
# R = \frac{\hat{\mu}_{ref}}{\hat{\rho}_\infty \hat{U}^2_\infty \hat{l}_{ref}}
# \end{equation}

# From the definition of the *Reynolds Number*,  $Re = \frac{1}{R}$.

# Therefore:
# 
# \begin{equation}
# \mu\frac{\partial^2 U}{\partial x^2} =
# \hat{\mu}\frac{\partial^2 \hat{U}}{\partial \hat{x}^2}Re
# \end{equation}

# ### Energy

# The primary term used in the *energy equation* is the *total enthalpy*:
# 
# \begin{equation}
# \frac{\partial\hat{\rho} \hat{U} \hat{H}}{\partial x} \hat{=}
# \frac{\hat{\rho}_\infty \hat{U}_\infty\hat{H}_{ref}}{\hat{l}_{ref}}
# \end{equation}
# 

# The definition of the *total enthalpy* is:
# 
# \begin{equation}
# H = h + \frac{1}{2}(u^2 + v^2 + w^2)
# \end{equation}
# 
# where $h$ is given by:
# 
# \begin{equation}
# h = c_p T
# \end{equation}

# In[ ]:




