""" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 'ErrorFormulae.py'
 Copyright (C) 2017 Jessica Kirkby-Kent
 Ver. 1.000 (For python 2.7)
 --Part of pyebop.py--
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 To show how most of the errors in pyebop are calculated....

"""
import math

def calFractionalErr(value, err):
	""" 
	Will work how much the error is as a fraction. e.g error/value for something
	 written as value +/- error(err)
	 
	"""
	fracErr = err / value
	return abs(fracErr)
	
def amountBasedOnFract(value, fracErr):
	""" 
	the opposite of calFractionalErr. Will calculate an amount of error, given
	the fractional error, and the value to which this error applies.
		
		PARAMETER:
			value - 
			fracErr - 
		RETURN:
			amount - to resulting amount of value that is error
	"""
	amount = fracErr * value
	return amount

def combineInQuad(x, y):
	""" 
	Will combine two values in quadridure
	
		RETURN - ans - the resulting answer.
	"""
	
	ans = math.sqrt(math.pow(x,2)+math.pow(y,2))
	return ans

def ErrWhenDiv(x, xErr, y, yErr):
	""" 
	Based on the assumption that calculating error for x/y
	
	return a list [value, err] where the error is not fractional or percentage
	"""
	xFrac = calFractionalErr(x, xErr)
	yFrac = calFractionalErr(y, yErr)
	quadFracErr = combineInQuad(xFrac, yFrac)
	value = x/y
	amountErr = amountBasedOnFract(value, quadFracErr)
	
	
	return [value, amountErr]

def powErr(fracErr, pow):

	""" 
	
	"""
	ans = fracErr * pow
	return ans
	
def atanErr(x,xErr,y,yErr):
	""" 
	Takes to parameter values and their associated error, and calculates the overall value
	 and the combined error. The error for the maximum of the answer +/- their uncertainies
	 and is given in degrees. Ans= atan(x/y) in degrees
	 
	 RETURNS:
	 	
	 	[Ans, OverallErr] - list of two values, the is the resulting answer and OverallErr
	 						it's associated error	
	"""
	Ans = math.degrees(math.atan2(x,y))
	Max = math.degrees(math.atan2((x+xErr),(y-yErr)))
	Min = math.degrees(math.atan2((x-xErr),(y+yErr)))	

	OverallErr = max(abs(Max-Ans),abs(Min-Ans))

	if OverallErr > 180:
		OverallErr = min(abs(Max-Ans),abs(Min-Ans))
		
	return [Ans, OverallErr]	
 
def sqrtxsqysq(x, xErr, y, yErr):
	""" 
	to calulate the error on something with the function:
		sqrt(x^2+y^2)
	"""
	value = math.sqrt(x**2+y**2)
	
	xErrFrac = calFractionalErr(x,xErr)
	yErrFrac = calFractionalErr(y,yErr)
	xPowErr = powErr(xErrFrac,2)
	yPowErr = powErr(yErrFrac,2)
	
	inQuad = combineInQuad(xPowErr, yPowErr)
	combinePow = powErr(inQuad, 0.5)
	
	AnsErr = amountBasedOnFract(value, combinePow)
	
	return [value, AnsErr]

def r1Err(rsum, rsumErr, kratio, kratioErr):
	""" 
	A specific calculation for the analyseLC code 
	
	"""
	value = rsum/(1+kratio)
	
	xErrFrac = calFractionalErr(rsum, rsumErr)
	yErrFrac = calFractionalErr(kratio, kratioErr)
	
	inQuad = combineInQuad(xErrFrac, yErrFrac)
	AnsErr = amountBasedOnFract(value, inQuad)
	
	return [value, AnsErr]
	
def r2Err(rsum, rsumErr, kratio, kratioErr):
	""" 
	A specific err calc for r2
	
	"""
	value = rsum/((1/kratio)+1)
	
	xErrFrac = calFractionalErr(rsum, rsumErr)
	yErrFrac = calFractionalErr(kratio, kratioErr)
	
	oneOverKErr = combineInQuad(0, yErrFrac)
	aboveNotFrac = amountBasedOnFract((1.0/kratio), oneOverKErr)
	abovePlus1 = combineInQuad(aboveNotFrac, 0)
	aboveFrac = calFractionalErr(((1/kratio)+1), abovePlus1)
	
	totCombineErr = combineInQuad(aboveFrac, xErrFrac)
	AnsErr = amountBasedOnFract(value, totCombineErr)
	
	return [value, AnsErr]
