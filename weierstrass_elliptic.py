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
		g3 = mpf(g3)

		self.__ng3 = False
		self.__Delta = 4 * 4 * g2*g2*g2 - 27 * 4*4 * g3*g3
		self.__invariants = (g2,g3)
		self.__roots = self.__compute_roots(mpf('4.0'),-g2,-g3)
		self.__omegas = self.__compute_omegas()
		self.__periods = self.__compute_periods()
		self.__maclaurin = self.__compute_maclaurin()
		#self.__etas = self.__compute_etas()

	def __compute_etas(self):
		N=3
		from mpmath import mpc,pi
		g2, g3 = self.__invariants
		c = self.__maclaurin
		om,omp,om2,om2p = self.__omegas
		z = omp

		for i in range(N):
			z = z/2

		Delta = self.__Delta 
		z02=z*z
		z03=z*z02
		z04=z02*z02
		z05=z03*z02
		z06=z04*z02
		z07=z05*z02
		z08=z06*z02
		z09=z07*z02
		z10=z08*z02
		z11=z09*z02
		etap = 1/z - c[2]*z03/3 - c[3]*z05/5-c[4]*z07/7-c[5]*z09/9-c[6]*z11/11-c[7]*z11*z02/13  -c[8]*z11*z02*z02/15-c[9]*z11*z02*z02*z02/17
		P = 1/z02 + c[2]*z02   + c[3]*z04 + c[4]*z06 + c[5]*z08 + c[6]*z06*z04 + c[7]*z06*z06 + c[8]*z08*z06      + c[9]*z08*z08
		Pprime = -2/z03 + 2*c[2]*z   + 4*c[3]*z03 + 6*c[4]*z05 + 8*c[5]*z07 + 10*c[6]*z09 + 12*c[7]*z11 + 14*c[8]*z07*z06  #    + 16*c[9]*z08*z07

		#  We apply the duplication formulas N times
		for j in range(N):
			P2 = P*P
			Ppprime = 6*P2-g2/2
			etap = 2*etap+Ppprime/(2*Pprime)
			Pprime = (-4*Pprime*Pprime*Pprime*Pprime+12*P*Pprime*Pprime*Ppprime-Ppprime*Ppprime*Ppprime)/(4*Pprime*Pprime*Pprime)
			P = - 2 * P + (6*P2-g2/2)*(6*P2-g2/2)/(4*P2*P-g2*P-g3)/4

		if Delta <0:
			eta = etap.conjugate()
		else:
			eta = (pi*mpc(0,1)/2+etap*om)/omp #Legendre Relation

		return eta,etap

	def zeta(self,z):
		N=3
		from mpmath import mpc,pi
		g2, g3 = self.__invariants
		c = self.__maclaurin
		om,omp,om2,om2p = self.__omegas

		for i in range(N):
			z = z/2

		Delta = self.__Delta 
		z02=z*z
		z03=z*z02
		z04=z02*z02
		z05=z03*z02
		z06=z04*z02
		z07=z05*z02
		z08=z06*z02
		z09=z07*z02
		z10=z08*z02
		z11=z09*z02
		zeta = 1/z - c[2]*z03/3 - c[3]*z05/5-c[4]*z07/7-c[5]*z09/9-c[6]*z11/11-c[7]*z11*z02/13  -c[8]*z11*z02*z02/15-c[9]*z11*z02*z02*z02/17
		P = 1/z02 + c[2]*z02   + c[3]*z04 + c[4]*z06 + c[5]*z08 + c[6]*z06*z04 + c[7]*z06*z06 + c[8]*z08*z06      + c[9]*z08*z08
		Pprime = -2/z03 + 2*c[2]*z   + 4*c[3]*z03 + 6*c[4]*z05 + 8*c[5]*z07 + 10*c[6]*z09 + 12*c[7]*z11 + 14*c[8]*z07*z06  #    + 16*c[9]*z08*z07

		#  We apply the duplication formulas N times
		for j in range(N):
			P2 = P*P
			Ppprime = 6*P2-g2/2
			zeta = 2*zeta+Ppprime/(2*Pprime)
			Pprime = (-4*Pprime*Pprime*Pprime*Pprime+12*P*Pprime*Pprime*Ppprime-Ppprime*Ppprime*Ppprime)/(4*Pprime*Pprime*Pprime)
			P = - 2 * P + (6*P2-g2/2)*(6*P2-g2/2)/(4*P2*P-g2*P-g3)/4
		return zeta

	def __compute_maclaurin(self):
		c = list([0]*10)
		g2, g3 = self.__invariants
		c[2] = g2/20
		c[3] = g3/28
		c[4] = c[2]*c[2]/3
		c23 = c[2]*c[2]*c[2]
		c32 = c[3]*c[3]
		c[5] = 3*c[2]*c[3]/11
		c[6] = (2*c23+3*c32)/39
		c[7] = 2*c[4]*c[3]/11
		c[8] = 5*c[2]*(11*c23+36*c32)/7239
		c[9] = c[3]*(29*c23+11*c32)/2717
		return c

	def __compute_roots(self,a,c,d):
		from mpmath import mpf, mpc, sqrt, cbrt
		assert(all([isinstance(_,mpf) for _ in [a,c,d]]))
		Delta = self.__Delta
		# NOTE: this was the original function used for root finding.
		# proots, err = polyroots([a,0,c,d],error=True,maxsteps=5000000)
		# Computation of the cubic roots.
		# TODO special casing.
		u_list = [mpf(1),mpc(-1,sqrt(3))/2,mpc(-1,-sqrt(3))/2]
		Delta0 = -3 * a * c
		Delta1 = 27 * a * a * d
		C1 = cbrt((Delta1 + sqrt(Delta1 * Delta1 - 4 * Delta0 * Delta0 * Delta0)) / 2)
		C2 = cbrt((Delta1 - sqrt(Delta1 * Delta1 - 4 * Delta0 * Delta0 * Delta0)) / 2)
		if abs(C1) > abs(C2):
			C = C1
		else:
			C = C2
		proots = [(-1 / (3 * a)) * (u * C + Delta0 / (u * C)) for u in u_list]
		# NOTE: we ignore any residual imaginary part that we know must come from numerical artefacts.
		if Delta < 0:
			# Sort the roots following the convention: complex with negative imaginary, real, complex with positive imaginary.
			# Then assign the roots following the P convention (e2 real root, e1 complex with positive imaginary).
			e3,e2,e1 = sorted(proots,key = lambda c: c.imag)
		else:
			# The convention in this case is to sort in descending order.
			e1,e2,e3 = sorted([_.real for _ in proots],reverse = True)
		return e1,e2,e3

	def __compute_omegas(self):
		# A+S 18.9.
		from mpmath import sqrt, ellipk, mpc, pi, mpf
		Delta = self.__Delta
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
		if self.__Delta >= 0:
			return 2 * om, 2 * omp
		else:
			return 2 *(om + omp), 2 * omp

	def __compute_user_periods(self):
		from mpmath import mpc
		Delta = self.__Delta
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
		return deepcopy(self.__compute_user_invariants())
	@property
	def etas(self):
		from copy import deepcopy
		return deepcopy(self.__compute_etas())
	@property
	def Delta(self):
		from copy import deepcopy
		return deepcopy(self.__Delta)
	@property
	def periods(self):
		from copy import deepcopy
		return deepcopy(self.__compute_user_periods())
	@property
	def omegas(self):
		from copy import deepcopy
		return deepcopy(self.__omegas)
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
		retval += 'Omegas:\t\t' + str(self.omegas) + '\n'
		return retval

	def P(self,z):
		# 0 - Number of iterations (halvings) is fixed
		# this could be adapted together with the number of terms retained
		# in the Laurent series to increase performance .. but how?		
		N = 4

		# 1 - we first reduce z to the fundamental period parallelogram (the one defined in Abramowitz)
		# and we call the new value z_r. P(z)=P(z_r) by definition of the FPP
		T1,T2 = self.__periods
		z_r = self.reduce_to_fpp2(z,T1,T2)

		# 2 - we then reduce z_r to 1/4 FPP 
		x = z_r.real
		y = z_r.imag
		om,omp,om2,om2p = self.__omegas
		
		if self.__Delta >= 0:
			if (x>=om and y>=omp.imag): 	#R3
				z_r = 2*om2-z_r
			elif (x>om): 					#R4
				z_r = (2*om-z_r).conjugate()
			elif (y>omp.imag): 				#R2
				z_r = (z_r-2*omp).conjugate()
			else: 							#R1
				z_r = z_r
		else:
			if (y>=0 and x<=om2.real):		#R1
				z_r = z_r
			if (y>=0 and x>om2.real):		#R2
				z_r = z_r.conjugate()
			if (y<=0 and x<=om2.real):		#R3
				z_r = 2*om-z_r
			if (y<=0 and x<=om2.real):		#R4
				z_r = (2*om-z_r).conjugate()
		
			# 2a - We reduce z_r to the fundamental rectangle
			if z_r.imag > om2p.imag / 2:		
				z_r = 2*omp - z_r
		
		# 3 - We half z_r N times (can this be done efficiently with bit operators?)
		for i in range(N):
			z_r = z_r/2

		# 4 - we compute the Laurent series to the 9th term in the N-halved z_r
		g2, g3 = self.__invariants
		c = self.__maclaurin

		z02=z_r*z_r
		z04=z02*z02
		z06=z04*z02
		z08=z04*z04
		P = 1/z02 + c[2]*z02 + c[3]*z04 + c[4]*z06 + c[5]*z08 + c[6]*z06*z04 + c[7]*z06*z06 + c[8]*z08*z06# + c[9]*z08*z08

		# 5 - We apply the duplication formula N times
		for j in range(N):
			P2 = P*P
			P = - 2 * P + (6*P2-g2/2)*(6*P2-g2/2)/(4*P2*P-g2*P-g3)/4

		# 6 - We return the appropriate value for P (according to the 1/4 FPP reduction formulas)
		if self.__Delta >=0:
			if (x>om.real and y>omp.imag): 	#R3
				return P
			elif (x>om.real): 				#R4
				return P.conjugate()
			elif (y>omp.imag): 				#R2
				return P.conjugate()
			else: 							#R1
				return P
		else:
			if (y>=0 and x<=om2.real):		#R1
				return P
			if (y<0 and x <= om2.real):		#R2
				return P.conjugate()
			if (y<=0 and x>om2.real):		#R3
				return P
			if (y>0 and x>om2.real):		#R4
				return P.conjugate()

	def Ptheta(self,z):
		from mpmath import pi, jtheta, exp, mpc
		e1,e2,e3 = self.__roots
		om,omp,om2,om2p = self.__omegas
		Delta = self.__Delta
		if Delta > 0:
			tau = omp/om
			q = exp(mpc(0,1)*pi*tau)
			v = (pi * z) / (2*om)
			omega = om
			e = e1
		else:
			tau2 = om2p/(2*om2)
			q = mpc(0,1)*exp(mpc(0,1)*pi*tau2)
			v = (pi * z) / (2*om2)
			omega = om2
			e = e2

		retval = e+pi**2/(4*omega*omega)*( jtheta(n=1,z=0,q=q,derivative=1) * jtheta(n=2,z=v,q=q,derivative=0) 
			/ ( jtheta(n=2,z=0,q=q,derivative=0) * jtheta(n=1,z=v,q=q,derivative=0) ) )**2
		return retval

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
	
		om,omp,om2,om2p = self.__omegas
		if g3<0: #DOES NOT WORK .... REDUCTION TO FPP PROBABLY FAILED
			Pprime_r = Pprime*1j
			z_r = z*1j
			om,omp,om2,om2p = -1j*omp, 1j*om,-1j*om2p,1j*om2
			if self.Delta >= 0:
				T1,T2 = 2 * om, 2 * omp
			else:
				T1,T2 = 2 *(om + omp), 2 * omp
			z_r = self.reduce_to_fpp2(z_r,T1,T2)
		else:
			z_r = z
			T1,T2 = self.__periods
			z_r = self.reduce_to_fpp2(z_r,T1,T2)
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

	def reduce_to_fpp2(self,z,T1,T2):
		from math import floor
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
