# -*- coding: iso-8859-1 -*-
# Copyright (C) 2014 by Dario Izzo
# dario.izzo@gmail.com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

class weierstrass_elliptic(object):

	def __init__(self,g2,g3):
		from mpmath import mpf
		g2 = mpf(g2)
		#Handle negative g3.
		#if g3 < 0:
		#	g3 = -mpf(g3)
		#	self.__ng3 = True
		#else:
		#	g3 = mpf(g3)
		self.__ng3 = False
		#assert(g3 >= 0)
		# Upon construction we compute/store the discriminant, the roots, the periods and more
		self.__Delta = 4 * 4 * g2*g2*g2 - 27 * 4*4 * g3*g3
		self.__invariants = (g2,g3)
		#self.__user_invariants = self.__compute_user_invariants()
		self.__roots = self.__compute_roots()
		self.__omegas = self.__compute_omegas()
		self.__periods = self.__compute_periods()
		#self.__user_periods = self.__compute_user_periods()

	def __compute_roots(self):
		from mpmath import  sqrt
		from mpmath import mpf, polyroots
		p = - self.__invariants[0] / 4.0
		q = - self.__invariants[1] / 4.0
		u1=1.0
		u2=(-1.0+sqrt(3)*1j)/2.0
		u3=(-1.0-sqrt(3)*1j)/2.0
		Delta0 = -3.0*p
		Delta1 = 27.0*q
			
		C1 = ((Delta1+sqrt(Delta1**2.0-4.0*Delta0**3.0))/2.0 +0j)**(1.0/3.0)
		C2 = ((Delta1-sqrt(Delta1**2.0-4.0*Delta0**3.0))/2.0 +0j)**(1.0/3.0)
	
		#This avoids loosing numerical precision when C->0
		if abs(C1)>abs(C2):
			C=C1
		else:
			C=C2

		e1 = -1.0/3.0*(u1*C+Delta0/u1/C)
		e2 = -1.0/3.0*(u2*C+Delta0/u2/C)
		e3 = -1.0/3.0*(u3*C+Delta0/u3/C)

		if self.__Delta > 0:
			e1,e2,e3 = sorted([e1.real,e2.real,e3.real], reverse = True)
		else:
			e1,e2,e3 = sorted([e1,e2,e3], key = lambda x: x.imag, reverse = True)
		return e1,e2,e3

	def __compute_omegas(self):
		# A+S 18.9.
		from mpmath import sqrt, ellipk, mpc, pi, mpf
		Delta = self.Delta
		e1, e2, e3 = self.__roots
		if Delta > 0:
			m = (e2 - e3) / (e1 - e3)
			Km = ellipk(m)
			Kpm = ellipk(1 - m)
			om = Km / sqrt(e1 - e3)
			omp = mpc(0,1) * om * Kpm / Km
			om2 = om + omp
			om2p = omp - om
		elif Delta < 0:
			# NOTE: the expression in the sqrt has to be real and positive, as e1 and e3 are
			# complex conjugate and e2 is real.
			H2 = (sqrt((e2 - e3) * (e2 - e1))).real
			assert(H2 > 0)
			m = mpf(1) / mpf(2) - 3 * e2 / (4 * H2)
			Km = ellipk(m)
			Kpm = ellipk(1 - m)
			om2 = Km / sqrt(H2)
			om2p = mpc(0,1) * Kpm * om2 / Km
			om = (om2 - om2p) / 2
			omp = (om2 + om2p) / 2
		else:
			g2, g3 = self.__invariants
			if g2 == 0 and g3 == 0:
				om = mpf('+inf')
				omp = mpc(0,'+inf')
			else:
				# NOTE: here there is no need for the dichotomy on the sign of g3 because
				# we are already working in a regime in which g3 >= 0 by definition.
				c = e1 / 2
				om = 1 / sqrt(12 * c) * pi()
				omp = mpc(0,'+inf')
			self.__om = om
			self.__omp = omp
			om2 = om + omp
			om2p = omp - om
		return om,omp,om2,om2p

	def __compute_periods(self):
		om,omp,om2,om2p = self.__omegas
		if self.Delta >= 0:
			return 2 * om, 2 * omp
		else:
			return 2 *(om + omp), 2 * omp

	def __compute_user_periods(self):
		from mpmath import mpc
		Delta = self.Delta
		# NOTE: here there is no need to handle Delta == 0 separately,
		# as it falls under the case of om purely real and omp purely imaginary.
		if Delta >= 0:
			T1, T3 = self.__periods
			if self.__ng3:
				return T3.imag, T1.real * mpc(0,1)
			else:
				return T1, T3
		else:
			T1, T3 = (sum(self.__periods)).real, self.__periods[1]
			if self.__ng3:
				return 2 * T3.imag, mpc(-T3.imag,T3.real)
			else:
				return T1, T3
	def __compute_user_invariants(self):
		if self.__ng3:
			return self.__invariants[0],-self.__invariants[1]
		else:
			return self.__invariants
	
	@property
	def invariants(self):
		from copy import deepcopy
		return deepcopy(self.__user_invariants)
	@property
	def Delta(self):
		from copy import deepcopy
		return deepcopy(self.__Delta)
	@property
	def periods(self):
		from copy import deepcopy
		return deepcopy(self.__user_periods)
	#@property
	#def omegas(self):
	#	from copy import deepcopy
	#	return deepcopy((self.__om,self.__omp,self.__om2,self.__om2p))
	@property
	def roots(self):
		from copy import deepcopy
		if self.__ng3:
			from functools import cmp_to_key
			# g3 < 0 means that all roots change sign.
			retval = [-x for x in self.__roots]
			# Sort by imaginary part, then real.
			return sorted(retval,key = cmp_to_key(lambda z1, z2: z1.imag - z2.imag if z1.imag != z2.imag else z1.real - z2.real),reverse=True)
		else:
			return deepcopy(self.__roots)
	def __repr__(self):
		retval = 'Invariants:\t' + str(self.invariants) + '\n'
		# NOTE: the Delta does not change for differences in sign of g_3.
		retval += 'Delta:\t\t' + str(self.Delta) + '\n'
		retval += 'Periods:\t' + str(self.periods) + '\n'
		retval += 'Roots:\t\t' + str(self.roots) + '\n'
		return retval

	def P(self,z):
		# 0 - Number of iterations is fixed
		# this could be adapted together with the number of terms retained
		# in the Laurent series to increase performance .. but how?		
		N = 4

		# 1 - we reduce to the fundamental period parallelogram (the one in Abramowitz)
		z_r = self.reduce_to_fpp2(z)

		# 2 - we reduce z_r first to 1/4 FPP to further reduce the error of the iterations
		# and then to the fundamental rectangle

		x = z_r.real
		y = z_r.imag
		om,omp,om2,om2p = self.__omegas
		
		if self.Delta >= 0:
			if (x>=om and y>=omp.imag): 	#R3
				z_r = 2.0*self.__om2-z_r
			elif (x>som): 				#R4
				z_r = (2.0*om-z_r).conjugate()
			elif (y>omp.imag): 			#R2
				z_r = (z_r-2.0*omp).conjugate()
			else: 						#R1
				z_r = z_r
		else:
			if (y>=0 and x<=om2.real):		#R1
				z_r = z_r
			if (y>=0 and x>om2.real):		#R2
				z_r = z_r.conjugate()
			if (y<=0 and x<=om2.real):		#R3
				z_r = 2.0*om-z_r
			if (y<=0 and x<=om2.real):		#R4
				z_r = (2.0*om-z_r).conjugate()
		
			if z_r.imag > om2p.imag / 2.0:		
				z_r = 2*omp - z_r
		
		# 3 - We half z_r N times
		for i in range(N):
			z_r = z_r/2.0
		# 4 - we now compute the Laurent series to the 9th term in the N-halved z_r
		g2, g3 = self.__invariants
		c2 = g2/20.0
		c3 = g3/28.0
		c23 = c2*c2*c2
		c32 = c3*c3
		c4 = c2*c2/3.0
		c5 = 3.0*c2*c3/11.0
		c6 = (2.0*c23+3.0*c32)/39.0
		c7 = 2.0*c4*c3/11.0
		c8 = 5.0*c2*(11.0*c23+36.0*c32)/7239.0
		c9 = c3*(29.0*c23+11.0*c32)/2717.0
		z02=z_r*z_r
		z04=z02*z02
		z06=z04*z02
		z08=z04*z04
		P = 1.0/z02 + c2*z02 + c3*z04 + c4*z04*z02 + c5*z04*z04 + c6*z06*z04 + c7*z06*z06 + c8*z08*z06 + c9*z08*z08

		# 5 - We apply the duplication formula N times
		for j in range(N):
			P2 = P*P
			P = - 2.0 * P + (6*P*P-1.0/2.0*g2)*(6.0*P2-1.0/2.0*g2)/(4.0*(4.0*P2*P-g2*P-g3))

		# 6 - We return the appropriate value for P (according to the 1/4 FPP reduction formulas)
		if self.Delta >=0:
			if (x>om.real and y>omp.imag): 	#R3
				return P
			elif (x>om.real): 			#R4
				return P.conjugate()
			elif (y>omp.imag): 			#R2
				return P.conjugate()
			else: 						#R1
				return P
		else:
			if (y>=0 and x<=om2.real):		#R1
				return P
			if (y<0 and x <= om2.real):		#R2
				return P.conjugate()
			if (y<=0 and x>om2.real):		#R3
				return P
			if (y>0 and x>om2.real):			#R4
				return P.conjugate()
		
	def Pprime(self,z):
		# 1 - We compute P
		P = self.P(z)
		# 2 - From P we compute Pprime
		return self.__Pprime_from_P(z,P)

	def __Pprime_from_P(self,z,P):
		from mpmath import sqrt,sign,phase,pi
		txt = ''
		# 1 - We start by computing Pprime without knowing its sign
		g2, g3 = self.__invariants
		Pprime = sqrt(4*P*P*P - g2*P - g3)

		# 2 - We now resolve the sign ambiguity using the method in (Abramowitz 18.8)
		# 2.1 - First we reduce to studying the behaviour in the fundamental parallelogram
		z_r = self.reduce_to_fpp2(z)

		om,omp,om2,om2p = self.__omegas
		if g3<0: #DOES NOT WORK .... REDUCTION TO FPP PROBABLY FAILED
			Pprime_r = Pprime*1j
			z_r = z_r*1j
			om,omp,om2,om2p = -1j*omp, 1j*om,-1j*om2p,1j*om2
		else:
			Pprime_r = Pprime


		# 2.2 - We then proceed separately according to the sign of Delta
		if self.Delta > 0:
			# 2.2.1 - We reduce the study of the sign to the fundamental rectangle (in this case 1/4 FPP)
			x = z_r.real
			y = z_r.imag
			if (x>om.real and y>omp.imag): 			#R3
				z_r = 2*om2-z_r
				Pprime_r =  - Pprime
			elif (x>om.real): 				#R4
				z_r = (2*om-z_r).conjugate()
				Pprime_r = -Pprime.conjugate()
			elif (y>omp.imag): 				#R2
				z_r = (z_r-2*omp).conjugate()
				Pprime_r = Pprime.conjugate()
			else: 						#R1
				Pprime_r = Pprime

			# 2.2.2 - We apply Abramowitz criteria
			z_r = z_r/om
			a = omp.imag/om.real
			x = z_r.real
			y = z_r.imag
			if (y<=0.4 and x<=0.4): 				#region B
				zmac = -2.0/z_r/z_r/z_r
				flag = (sign(zmac.real)==sign(Pprime_r.real))
			elif (y>=0.4 and x<=0.5): 				#region A (first case)	
				flag = (Pprime_r.real >= 0)		
			else: 							#elswhere
				flag = (Pprime_r.imag >= 0)
			if not flag:
				Pprime = -Pprime
		else:
			# 2.2.1 - We reduce the study of the sign to the fundamental rectangle 
			# using formulas 18.2.22 - 18.2.33 from Abramowitz
			x = z_r.real
			y = z_r.imag
			if (x<=om2.real and y>=0.0): 		#R1
				Pprime_r =  Pprime
				txt = 'R1'
			elif (x<=om2.real and y<0.0): 		#R2
				z_r = z_r.conjugate()
				Pprime_r = Pprime.conjugate()
				txt = 'R2'
			elif (x>om2.real and y<=0.0): 		#R3
				z_r = 2*om2 - z_r
				Pprime_r = - Pprime
				txt = 'R3'
			else: 					#R4
				z_r = (2*om2 - z_r).conjugate()
				Pprime_r = - Pprime.conjugate()
				txt = 'R4'
			if z_r.imag > om2p.imag / 2.0:		#We move from 1/4 FPP to fundamental rectangle	
				z_r = 2*omp - z_r
				Pprime_r = -Pprime_r
				txt = txt + ' - FR'
			

			# 2.2.2 - We then apply Abramowitz criteria by first reducing to the case
			# where the real half-period is unity
			z_r = z_r/om2
			
			# and by then applying the criterion
			a = om2p.imag/om2.real
			x = z_r.real
			y = z_r.imag

			if (y<=0.4 and x<=0.4): 				#region B
				zmac = -2.0/z_r/z_r/z_r
				flag = (sign(zmac.real)==sign(Pprime_r.real))
				txt = txt + ' - Region B'
			elif a>=1.05:
				if (y>=0.4 and x<=0.5): 			#as region A of Delta>0
					flag = (Pprime_r.real >= 0)
					txt = txt + ' - Region A real'				
				else: 		
					flag = (Pprime_r.imag >= 0)
					txt = txt + ' - Region A imag'
				
			elif a>1:							
				if(y>=0.4 and x<=0.4):
					flag = (Pprime_r.real >= 0)
					txt = txt + ' - case 1'
				elif (x >0.4 and x <=0.5 and y >0.4 and y <=0.5):
					arg = phase(Pprime_r)
					flag = (arg > -pi/4 and arg < 3.0/4.0*pi)
					txt = txt + ' - case 2 (arg)'
				else:
					flag = (Pprime_r.imag >= 0)
					txt = txt + ' - case 3 (elsewhere)'
			else:
				print "BOH"			

			if not flag:
				Pprime = -Pprime
		print txt
		return Pprime
	
	def Pinv(self,P):
		from mpmath import ellipf, sqrt, asin, acos, mpc, mpf
		Delta = self.Delta
		e1, e2, e3 = self.__roots
		if self.__ng3:
			P = -P
		if Delta > 0:
			m = (e2 - e3) / (e1 - e3)
			retval = (1 / sqrt(e1 - e3)) * ellipf(asin(sqrt((e1 - e3)/(P - e3))),m=m)
		elif Delta < 0:
			H2 = (sqrt((e2 - e3) * (e2 - e1))).real
			assert(H2 > 0)
			m = mpf(1) / mpf(2) - 3 * e2 / (4 * H2)
			retval = 1 / (2 * sqrt(H2)) * ellipf(acos((e2-P+H2)/(e2-P-H2)),m=m)
		else:
			g2, g3 = self.__invariants
			if g2 == 0 and g3 == 0:
				retval = 1 / sqrt(P)
			else:
				c = e1 / 2
				retval = (1 / sqrt(3 * c)) * asin(sqrt((3 * c)/(P + c)))
		if self.__ng3:
			retval /= mpc(0,1)
		alpha, beta, _, _ = self.reduce_to_fpp(retval)
		T1, T2 = self.periods
		return T1 * alpha + T2 * beta

	def reduce_to_fpp2(self,z):
		from math import floor
		T1, T2 = self.__periods
		a, b = T1.real, T2.real
		d, c = T1.imag, T2.imag
		assert(d == 0)
		b1,b2 = z.real,z.imag
		x2 = b2 / c
		x1 = (b1 - b*x2) / a
		N = int(floor(x1))
		M = int(floor(x2))
		alpha = x1 - N
		beta = x2 - M
		assert(alpha >= 0 and beta >= 0)
		retval =  alpha*T1+beta*T2
		if self.Delta > 0:
			return retval
		else:
			om,omp,om2,om2p = self.__omegas
			if retval.real < om2.real:
				return retval
			if retval.real > 2*om2.real:
				return retval-2*omp
			zz = retval-2*omp
			if zz.imag/zz.real >= om.imag/om.real:
				return zz
			else:
				return retval
	def reduce_to_fpp3(self,z):
		from math import floor
		T1, T2 = self.__periods
		a, b = T1.real, T2.real
		d, c = T1.imag, T2.imag
		assert(d == 0)
		b1,b2 = z.real,z.imag
		x2 = b2 / c
		x1 = (b1 - b*x2) / a
		N = int(floor(x1))
		M = int(floor(x2))
		alpha = x1 - N
		beta = x2 - M
		assert(alpha >= 0 and beta >= 0)
		return  alpha*T1+beta*T2

				
				

	def reduce_to_fpp(self,z):
		# TODO: need to check what happens here when periods are infinite.
		from mpmath import floor, matrix, lu_solve
		T1, T2 = self.__periods
		R1, R2 = T1.real, T2.real
		I1, I2 = T1.imag, T2.imag
		A = matrix([[R1,R2],[I1,I2]])
		b = matrix([z.real,z.imag])
		x = lu_solve(A,b)
		N = int(floor(x[0]))
		M = int(floor(x[1]))
		alpha = x[0] - N
		beta = x[1] - M
		assert(alpha >= 0 and beta >= 0)
		return alpha,beta,N,M

	def plot_infpp(self,z):
		import matplotlib.pyplot as pl
		om,omp,om2,om2p = self.__omegas
		pl.plot(2* om.real,2* om.imag,'ok')
		pl.plot(2* omp.real,2*omp.imag,'ok')
		pl.plot(2* om2.real,2*om2.imag,'ok')
		pl.plot(z.real,z.imag,'or')
		pl.plot(0,0,'ok')

		pl.show()
