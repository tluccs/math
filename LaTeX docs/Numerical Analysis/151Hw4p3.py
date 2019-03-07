import numpy as np

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

def f(x):
	return np.cos(x)

a = 0
b = np.pi
N = 200 # N must be even for simpsons to work. 
# Since we have N intervals, we get N+1 (x,f(x)) tuples, and composite simpsons works 
# For odd number of inputs

intervals =  createSubintervals(a, b, N)
pts = [ [x[0], f(x[0]) ]  for x in intervals]

pts.append( [intervals[-1][1], f(intervals[-1][1])] )
# Comments are outputs
print(integrateTrap(pts)) #1.15220333274e-14

print(integrateSimpsonsEquidistant(pts)) #5.85031809399e-15
