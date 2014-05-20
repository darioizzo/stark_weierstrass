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
nn = 100
wpOLD = weOLD.weierstrass_elliptic(g2,g3)
print "g2: " + str(g2) + ", g3: " + str(g3)

for i in range(nn):
	z = mpc(10 - rand()*20,10 - rand()*20)
	a = wp.Pprime(z)
	b = wpOLD.Pprime(z)
	#a = wp.P(z)
	#b = wpOLD.P(z)
	err = a-b
	err = sqrt((err*err.conjugate()).real).real
	print err
	sol.append(err)
	if err > 1e-2:
		print z,wp.Delta
		print " "

print g2,g3
print mean(sol), std(sol), float(sum([s<1e-2 for s in sol]))/nn

