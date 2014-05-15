from cmath import sqrt


from mpmath import mpc
from numpy import mean,std
from numpy.random import rand,seed
import weierstrass_elliptic as we
import weierstrass_ellipticOLD as weOLD
import matplotlib.pyplot as pl
#seed(9)
flag = True
while flag:
	g2 = 1 - 2*rand()
	g3 = 1 - 2*rand()
	wp = we.weierstrass_elliptic(g2,g3)
	if wp.Delta < 0: 
		flag = False

sol=list()
nn = 1
g2,g3 = 0.00114640883045, -0.294460277819
wpOLD = weOLD.weierstrass_elliptic(g2,g3)
wp = we.weierstrass_elliptic(g2,g3)
z =  mpc(real='6.2799995306229794', imag='5.201543383152698')
for i in range(nn):
	z =  mpc(real='6.2799995306229794', imag='5.201543383152698')
	#a = wp.Pprime(z)
	#b = wpOLD.Pprime(z)
	a = wp.P(z)
	b = wpOLD.P(z)
	err = a-b
	err = sqrt((err*err.conjugate()).real).real
	print err
	sol.append(err)
	if err > 1e-2:
		print z,wp.Delta
		print " "

print g2,g3
print mean(sol), std(sol), float(sum([s<1e-2 for s in sol]))/nn

