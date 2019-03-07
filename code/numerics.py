### Numerical Methods
from __future__ import print_function
from numpy.polynomial import Polynomial as P
import numpy as np
import matplotlib.pyplot as plt

class Matrix:
	def __init__(self, r, c):
		self.rows = r
		self.cols = c
		self.data = np.zeros((r,c))
	
	def pivot(self, r):
		"""
		Pivots in row r, for each row after r, first el turns to 0.
		changes self.data
		returns matrix P s.t. self.data = P*self.data
		"""
		m = Matrix(self.rows, self.cols)
		I = np.eye(self.rows, self.cols)
		for i in range(r+1, self.rows):
			I[i,r] = -1.0*self.data[i,r]/self.data[r,r]
			self.data[i,:] = self.data[i,:] + I[i,r] * self.data[r,:]
		m.set_data(I)
		return m
		
		
	def print_data(self):
		for i in range(self.rows):
			for j in range(self.cols):
				print(self.data[i,j], end=',')
			print("")
	
	def set_data(self, data):
		self.data =  np.array(data)

	def get_cols(self):
		return self.cols

	def get_rows(self):
		return self.rows

	def get_data(self):
		return np.array(self.data)

	def mult(self, A):
		"""
		returns B = self * A , does not change self or A
		"""
		r,c  = self.rows, A.get_cols()
		B = Matrix(r,c)
		D1 = np.array(self.get_data())
		D2 = np.array(A.get_data())
		B.set_data(np.matmul(D1, D2))
		return B

	def invert(self):
		"""
		Inverts self, changes values
		"""
		D = np.array(self.get_data())
		self.set_data(np.linalg.inv(D))


# fixed point method
def findFixedPt(f, start, iterations=10, tolerance=0.01):
	#apply f to p0 iterations times, or until error is within tolerance
	#assume there is a fixed point
	prevpt = start
	newpt = f(prevpt)
	count = 1
	error = -1
	while (count < iterations):
		error = abs(newpt - prevpt)
		if (error < tolerance):
			return newpt
		newpt, prevpt = f(newpt), newpt
		count = count + 1
		if (abs(newpt) > 100000 ):
			# give up 
			print("Give up!")
			return newpt
	return newpt

# binary search for zero
def findZero(f, interval, subintervals=100, tolerance=0.01):
	# if there is no zero, create subintervals intervals and check for a zero
	a = interval[0]
	b = interval[1]
	intervals = [[a + 1.0*i*(b-a)/subintervals, a + 1.0*(i+1)*(b-a)/subintervals] for i in range(subintervals)] 
	for sub in intervals:
		a, b = sub
		fa, fb = f(a), f(b)
		if ((fa < 0) ^ (fb < 0) or abs(fa) < tolerance or abs(fb) < tolerance):
			#theres a zero in this subinterval, so find via binary search
			while( (b - a) > tolerance):
				if ( abs(fa) < tolerance ):
					return a
				elif ( abs(fb) < tolerance ):
					return b
				else:
					c = (a+b) /2.0
					fc = f(c)
					if ((fa < 0) ^ (fc < 0)):
						b = c
						fb = f(b)
					else:
						a = c
						fa = f(a)
			return a
	#found no zeros if here
	print("Give up!")
	return 0

# find poly interpolation of function
def findPolynomialInterpolationf(f, pts):
	xdata = pts
	ydata = [ f(p) for p in xdata ]
	# use nevilles method
	fcns = [ P(y) for y in ydata] 
	# initialise j, increment each time to indicate gij = p_ i,i+1,...i+j
	j = 0
	# keep combining until we get one function
	while (len(fcns) > 1):
		j = j + 1
		newfcns = []
		print("step " + str(j))
		print(fcns)
		for i in range(len(fcns)-1):
			#method to combine 2 polys into one, to interpolate 1 extra point
			a = xdata[i]
			b = xdata[i+j]
			p1 = P([-a,1])
			p2 = P([-b,1])
			combined = (p1*fcns[i+1] - p2*fcns[i] ) / (b-a)
			newfcns.append(combined)
		fcns = newfcns
	return fcns[0]

# same as above but now pts is 2d of x,y data
def findPolynomialInterpolation(pts):
	xdata = [ pts[i][0] for i in range(len(pts)) ]
	ydata = [ pts[i][1] for i in range(len(pts)) ]
	# use nevilles method
	fcns = [ P(y) for y in ydata] 
	# initialise j, increment each time to indicate gij = p_ i,i+1,...i+j
	j = 0
	# keep combining until we get one function
	while (len(fcns) > 1):
		j = j + 1
		newfcns = []
		for i in range(len(fcns)-1):
			#method to combine 2 polys into one, to interpolate 1 extra point
			a = xdata[i]
			b = xdata[i+j]
			p1 = P([-a,1])
			p2 = P([-b,1])
			combined = (p1*fcns[i+1] - p2*fcns[i] ) / (b-a)
			newfcns.append(combined)
		fcns = newfcns
	return fcns[0]

#finds the linear spline of f
def findLinearSplinef(f, pts):
	# assume pts sorted
	def g(x):
		if x < pts[0]:
			print ("Lower than range!")
			i = 0
		elif  x > pts[-1]:
			print ("Greater than range!")
			i = -1
		else: 
			i = 0
			while (i < len(pts) and x >= pts[i]):
				i = i + 1
		# x is in [ pts[i-1], pts[i] ]
		if (i == len(pts)):
			i = i - 1
		a, b = pts[i-1], pts[i]
		fa, fb = f(a), f(b)
		return 1.0*(x - a)/(b-a)*(fb-fa) + fa
	return g

# same as above but now pts is 2d of x,y data 
def findLinearSpline(pts):
	# assume pts sorted
	xdata = [ pts[i][0] for i in range(len(pts)) ]
	ydata = [ pts[i][1] for i in range(len(pts)) ]
	def g(x):
		if x < xdata[0]:
			print ("Lower than range!")
			i = 0
		elif  x > xdata[-1]:
			print ("Greater than range!")
			i = -1
		else: 
			i = 0
			while (i < len(xdata) and x >= xdata[i]):
				i = i + 1
		# x is in [ pts[i-1], pts[i] ]
		if (i == len(xdata)):
			i = i - 1
		a, b = xdata[i-1], xdata[i]
		fa, fb = ydata[i-1], ydata[i]
		return 1.0*(x - a)/(b-a)*(fb-fa) + fa
	return g

# finds the polynomial derivative of f
def findPolyDerivativef(f, pts):
	pn = findPolynomialInterpolation(f, pts)
	coeffs = pn.coef
	pnprime = [ i*coeffs[i] for i in range(len(coeffs)) ]
	return P(pnprime[1:])
	
# same as above but for pts having x,y data
def findPolyDerivative(pts):
	pn = findPolynomialInterpolation(pts)
	coeffs = pn.coef
	pnprime = [ i*coeffs[i] for i in range(len(coeffs)) ]
	return P(pnprime[1:])
	
# Finds the function poly function Li(xj) = (i==j)
def findPolyKronecker(pts, pt):
	if pt not in pts:
		print ("Error: pt not in pts, returning 0")
		return P([0])
	denoms = [ (pt-ptj) for ptj in pts if ptj != pt ]
	denominator = np.prod(denoms)
	const = 1.0/ denominator
	polys = [ P([-ptj, 1]) for ptj in pts if ptj != pt  ]
	combined = np.prod(polys) * const
	return combined

# pts is x,y1, x+h,y2, x+2h,y3; finds derivative at x, x+h, x+2h
def findThreePointDerivative(pts):
	h1 = pts[1][0] - pts[0][0] 
	h2 = pts[2][0] - pts[1][0] 
	if (h1 != h2):
		print ("Low accuracy since data is not formatted as: (x,y1), (x+h,y2), (x+2h,y3)")
	h = (h1+h2) /2.0
	dx0 = (-3*pts[0][1] + 4*pts[1][1] - pts[2][1]) / (2*h)
	dx1 = (pts[2][1] - pts[0][1]) / (2*h)
	dx2 = (3*pts[2][1] - 4*pts[1][1] + pts[0][1]) / (2*h)
	return [dx0, dx1, dx2]
	
def dividedDifferencesf(f, pts):
	xdata = pts
	ydata = [ f(p) for p in pts ]
	# use nevilles method
	fs = [ y for y in ydata] 
	# initialise j, increment each time to indicate gij = p_ i,i+1,...i+j
	j = 0
	# cache contains everything we need
	cache = [fs[0]]
	# keep combining until we get last f
	while (len(fs) > 1):
		j = j + 1
		newfs = []
		print("Step " + str(j))
		print(fs)
		for i in range(len(fs)-1):
			#method to combine 2 fs into one
			a = xdata[i]
			b = xdata[i+j]
			combined = (fs[i+1] - fs[i] ) / (b-a)
			newfs.append(combined)
		fs = newfs
		cache.append(fs[0])
	return cache

def integratePoly(coeffs):
	newCoeffs = [0]
	for i in range(len(coeffs)):
		newCoeffs.append(1.0*coeffs[i]/(i+1))
	return newCoeffs

def evaluatePoly(coeffs, a, b=None):
	val = 0
	if (b is None):
		for i in range(len(coeffs)):
			val += coeffs[i]*(a**i)
	else:
		for i in range(len(coeffs)):
			val += coeffs[i]*(b**i) - coeffs[i]*(a**i)
	return val

def integrateTrap(pts):
	xdata = [ pt[0] for pt in pts]
	ydata = [ pt[1] for pt in pts]
	val = 0
	for i in range(len(xdata)-1):
		intervalLength = xdata[i+1] - xdata[i]
		avgHeight = (ydata[i+1] + ydata[i])/2.0
		val += intervalLength*avgHeight
	return val

def integrateSimpsonsEquidistant(pts):
	xdata = [ pt[0] for pt in pts]
	h = xdata[1] - xdata[0] # assume h fixed for x_(i+1) - x_i
	ydata = [ pt[1] for pt in pts]
	val = -ydata[0] + ydata[-1] # - here  bc we add 2y_0 in the for loop
	evens = range( int((len(xdata)-1)/2) )
	evens = np.array(evens) * 2
	for i in evens:
		val += 4*ydata[i+1] + 2*ydata[i]
	return val/3.0*h

def createSubintervals(a, b, N):
	step = 1.0*(b-a)/N
	curr = a
	nextpt = curr+step
	subints = []
	while (nextpt <= b ):
		subints.append([curr, nextpt])
		curr = nextpt
		nextpt = curr + step
	if(len(subints) != N):
		subints.append([curr, b])
	return subints
	
#test f
def f(x):
	#return x**2+x-1
	return np.sin(x+10)

def fprime(x):
	#return 2*x+1
	return np.cos(x+10)

def invx(x):
	return 1.0/x


"""
x = findFixedPt(f, 0.5, 10)
print("fixed point is " + str(x))
print("test: f(x) is " + str(f(x)))

z = findZero(f, [-1, 1], tolerance=0.0001)
print("zero is " + str(z))


print("\n")
sp = findLinearSpline(f, range(5))
for i in range(-1, 6):
	print ("spline is " + str(sp(i) ))
	print ("f is " + str(f(i) ))

pn = findPolynomialInterpolation(f, range(10))
pncoeffs = pn.coef
print("approximate poly is: ")
print (str(pncoeffs[0]), end="")
for i in range(1, len(pncoeffs)):
	print ("+" + str(pncoeffs[i]) + "x**" + str(i), end="")

pnprime = findDerivative(f, range(10))
pncoeffs = pnprime.coef
print("approximate poly derivative is: ")
print (str(pncoeffs[0]), end="")
for i in range(1, len(pncoeffs)):
	print ("+" + str(pncoeffs[i]) + "x**" + str(i), end="")
"""
"""p = findPolynomialInterpolationf(invx, [1, 2, 3])
print (p)
for i in range(3):
	print (p(i+1))
	"""

Adata = [[1, 2, 3, 4], [2, 5, 7, 9], [3, 7, 11, 14], [4, 9, 14, 19] ]
A = Matrix(4, 4)
A.set_data(Adata)
B = Matrix(4, 4)
B.set_data(Adata)
L1 = A.pivot(0)
L2 = A.pivot(1)
L3 = A.pivot(2)
L3.print_data()
print("L^ Av")
A.print_data()
Linv = L3.mult(L2).mult(L1)
Linv.invert()
Linv.print_data()
L1.invert()
L1.print_data()
L2.invert()
L3.invert()
L1.mult(L2).mult(L3).print_data()
