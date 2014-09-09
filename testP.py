from cmath import sqrt


from mpmath import mpc
from numpy import mean,std
from numpy.random import rand,seed
import weierstrass_elliptic as we
import weierstrass_ellipticOLD as weOLD
import matplotlib.pyplot as pl
import mpmath
import time


mpmath.mp.dps = 30

g2 = mpmath.mpf(1 - 2*rand())
g3 = mpmath.mpf(1 - 2*rand())

wp = we.weierstrass_elliptic(g2,g3)
wpOLD = weOLD.weierstrass_elliptic(g2,g3)

sol=list()
nn = 100

print "g2: " + str(g2) + ", g3: " + str(g3) + ", Delta: " + str(wp.Delta)
print "\nTesting Weierstrass P new implementation: a) study of the error"

for i in range(nn):
	z = mpc(10 - rand()*20,10 - rand()*20)
	mpmath.mp.dps = 15
	a = wp.P(z)
	mpmath.mp.dps = 30
	b = wpOLD.P(z)
	err = a-b
	err = sqrt((err*err.conjugate()).real).real
	sol.append(err)

print "Avg: " + str(mean(sol))
print "Std: " + str(std(sol))
print "Max: " + str(max(sol))

print "\nTesting Weierstrass P new implementation: a) study of the speed"
mpmath.mp.dps = 16
start_time = time.time()
for i in range(nn):
	z = mpc(10 - rand()*20,10 - rand()*20)
	a = wp.P(z)
print("New Implementation --- %s seconds ---" % (time.time() - start_time) )

start_time = time.time()
for i in range(nn):
	z = mpc(10 - rand()*20,10 - rand()*20)
	a = wpOLD.P(z)
print("Old Implementation --- %s seconds ---" % (time.time() - start_time) )
