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
sol=list()
sol2=list()
nn = 100


print("\nTesting Weierstrass P new implementation on %i random complex numbers" % nn)
print("\na) study of the error (ground truth is computed with 30 digits precision)")
print("   in parenthesis the error of the Jacobi based computations for P")

for i in range(nn):
	# We sample at random g2, g3 and z
	g2 = mpmath.mpf(1 - 2*rand())
	g3 = mpmath.mpf(1 - 2*rand())
	wp = we.weierstrass_elliptic(g2,g3)
	wpOLD = weOLD.weierstrass_elliptic(g2,g3)
	z = mpc(10 - rand()*20,10 - rand()*20)

	# We compute with 15 digits precision the P using the new method and the one based on Jacobi's functions
	mpmath.mp.dps = 15
	a = wp.P(z)
	c = wpOLD.P(z)

	# We compute with 30 digits precision the P to be used as ground truth
	mpmath.mp.dps = 30
	b = wpOLD.P(z)

	# We store the errors
	err = a-b
	err = sqrt((err*err.conjugate()).real).real
	sol.append(err)

	err = c-b
	err = sqrt((err*err.conjugate()).real).real
	sol2.append(err)

print("Avg: {} ({})".format(mean(sol),mean(sol2)) )
print("Std: {} ({})".format(std(sol),std(sol2)) )
print("Max: {} ({})".format(max(sol),max(sol2)) )

print("\nb) study of the execution speed (15 digits are used for both implementations)")
mpmath.mp.dps = 15
seed(123)
start_time = time.time()
for i in range(nn):
	g2 = mpmath.mpf(1 - 2*rand())
	g3 = mpmath.mpf(1 - 2*rand())
	z = mpc(10 - rand()*20,10 - rand()*20)
	wp = we.weierstrass_elliptic(g2,g3)
	a = wp.P(z)
print("New Implementation --- %s seconds ---" % (time.time() - start_time) )

seed(123)
start_time = time.time()
for i in range(nn):
	g2 = mpmath.mpf(1 - 2*rand())
	g3 = mpmath.mpf(1 - 2*rand())
	z = mpc(10 - rand()*20,10 - rand()*20)
	wpOLD = weOLD.weierstrass_elliptic(g2,g3)
	a = wpOLD.P(z)
print("Old Implementation --- %s seconds ---" % (time.time() - start_time) )
