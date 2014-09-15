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
n_trials = 3000

sol_dupl=list([0]*n_trials)
sol_theta=list([0]*n_trials)
sol_jacobi=list([0]*n_trials)
sol_refer=list([0]*n_trials)



print("\nTesting Weierstrass P new implementation on %i randomly sampled complex numbers" % n_trials)
print("\na) study of the error (ground truth is computed with 30 digits precision)")
print("   when the precision is set to 15 digits")

random_seed = 675637 
seed(random_seed)

for i in range(n_trials):
	# We sample at random g2, g3 and z
	g2 = mpmath.mpf(1 - 2*rand())
	g3 = mpmath.mpf(1 - 2*rand())
	wp = we.weierstrass_elliptic(g2,g3)
	wpOLD = weOLD.weierstrass_elliptic(g2,g3)
	z = mpc(10 - rand()*20,10 - rand()*20)

	# We compute with 15 digits precision the P using the new method and the one based on Jacobi's functions
	mpmath.mp.dps = 15
	a = wp.P(z)
	b = wp.Ptheta(z)
	c = wpOLD.P(z)

	# We compute with 30 digits precision the P to be used as ground truth
	mpmath.mp.dps = 30
	d = wpOLD.P(z)

	# We store the errors
	err = a-d
	err = sqrt((err*err.conjugate()).real).real
	sol_dupl[i] = err

	err = b-d
	err = sqrt((err*err.conjugate()).real).real
	sol_theta[i] = err

	err = c-d
	err = sqrt((err*err.conjugate()).real).real
	sol_jacobi[i] = err

print("Avg: Jacobi:{} Duplication:{} Theta:{}".format(mean(sol_jacobi),mean(sol_dupl),mean(sol_theta)) )
print("Std: Jacobi:{} Duplication:{} Theta:{}".format(std(sol_jacobi),std(sol_dupl),std(sol_theta)) )
print("Max: Jacobi:{} Duplication:{} Theta:{}".format(max(sol_jacobi),max(sol_dupl),max(sol_theta)) )

print("\nb) study of the execution speed (15 digits are used for both implementations)")
mpmath.mp.dps = 15

seed(random_seed)
start_time = time.time()
for i in range(n_trials):
	g2 = mpmath.mpf(1 - 2*rand())
	g3 = mpmath.mpf(1 - 2*rand())
	z = mpc(10 - rand()*20,10 - rand()*20)
	wpOLD = weOLD.weierstrass_elliptic(g2,g3)
	a = wpOLD.P(z)
print("Jacobi's formulas --- %s seconds ---" % (time.time() - start_time) )

seed(random_seed)
start_time = time.time()
for i in range(n_trials):
	g2 = mpmath.mpf(1 - 2*rand())
	g3 = mpmath.mpf(1 - 2*rand())
	z = mpc(10 - rand()*20,10 - rand()*20)
	wp = we.weierstrass_elliptic(g2,g3)
	a = wp.P(z)
print("Duplication Formulas --- %s seconds ---" % (time.time() - start_time) )



seed(random_seed)
start_time = time.time()
for i in range(n_trials):
	g2 = mpmath.mpf(1 - 2*rand())
	g3 = mpmath.mpf(1 - 2*rand())
	z = mpc(10 - rand()*20,10 - rand()*20)
	wp = we.weierstrass_elliptic(g2,g3)
	a = wp.Ptheta(z)
print("Thetas formulas --- %s seconds ---" % (time.time() - start_time) )
