variables: p T rho c_p h T p u v c a b c

constants R H gamma

given R^*,gamma -> c_p = gamma/(gamma - 1)R^*		unknowns

given H^* = h + V^2/2								h,u,v

a = sqrt(gamma * p/rho) (speed of sound)			h,u,v,p,rho

h = gamma/(gamma-1) p/rho = c_p T					T,u,v,p,rho

state: p = rho R^* T -> p/rho = R^* T 

PNS: A = rho u
     B = rho u^2 + p
     C = rho u v

M_x = u/a = u/sqrt((gamma-1)*c_p*T)

u^2 = (gamma - 1)*c_p*T*M_x^2

rho*u^2 = (gamma - 1) c_p*rho*T*M_x^2
		= (gamma - 1) (gamma/(gamma - 1)*R^*)*rho*T*M_x^2
        = gamma*p*M_x^2

p = B - gamma*p*M_x^2

p = B/(1+gamma*M_x^2)

H = c_p*T + 1/2*((gamma-1)*c_p*T*M_x^2 + v^2)

K = H - 1/2*v^2 = c_p*t + 1/2*(gamma - 1)*c_p*T*M_x^2

K = c_p*T(1 + 1/2*(gamma - 1)*M_x^2)

T =  K/(c_p*(1 + 1/2*(gamma - 1)*M_x^2))

c_p*rho*T = K*rho/(1 + 1/2*(gamma-1)*M_x^2)

(gamma/(gamma-1)*R^* * rho * T = K*rho/(1 + 1/2*(gamma - 1)*M_x^2)

(gamma/(gamma - 1)p = k*rho/(1 + 1/2*(gamma - 1)*M_x^2)

