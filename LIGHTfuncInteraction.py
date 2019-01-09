""" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 'LIGHTfuncInteraction.py'
 Copyright (C) 2017 Jessica Kirkby-Kent
 Contact: j.kirkbykent@gmail.com
 Ver. 1.000 (For python 2.7)
  --Part of pyebop.py--
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 This script will be used to request the initial setup from a user. It also has
  the main function that interaxts with the fortarn 'dlight' script. This 
  fortran script is what generates a synthetic lightcurve model given a set of 
  parameters.
 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 FUNCTIONS:
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 checkUser_YN_Input(userInput): 
   Takes a string and checks whether it matches
    suitable 'y', 'n', etc. Will continue to loop until an input that matches
    is entered. 
 
 float_check(string):
   Attempts to convert 'string' into a float, if this is unsuccessful, asks 
    user for a new input.
 
 YN_to_bool(user_ans)
   Takes a string from and if the string matches "y", "Y", "yes", "Yes", then
    it with return "True", otherwise it will return "False"
 
 chk_val_in_range(testVal, lowerLimit=0, upperLimit=None)
   Thake three value, the 1st is the number that you want to check to see if
    it is in a certain range. The 2nd and 3rd number defines the lower and 
    upper limits of the range respectively. If the 1st number is not in the
    range, will ask user for a new value. upperLimit=None means there is no
    upper limit.
    
 bigger_than_chk(testVal, lowerLim=0)
   Will check whether or not testVal is larger than lowerLim (which equals 0
    by default). If testVal does not pass then user must enter a new value.
 
 calcMag(v, phase)
   Calculate a model value for the phase, based on the parameters sorted in
    the list v.
    
 phase_mag_plot(phase, mag)
   Will plot the first value 'phase' against the second 'mag' to produce a 
    lightcurve, with appropriate labels applied to the axes
    
 generate_values(incl, e_cosw, e_sinw, u_p, u_s, y_grav_p, y_grav_s, srefl_p, 
 					srefl_s, qratio, tangle, thirdlight, phasecor, magzero
 					dgamma, phasemin, phasemax, nph)
   Will go through the process of asking the user for initial values to 
    generate a model lightcurve. Most parameters have a default if no specific
    value is given.
   					

 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 PARAMETERS:
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 
 USER SET:
 
 sbratio    -   CENTRAL SURFACE BRIGHTNESS
 		    -   VALUES: must be greater than 0
 
 rsum       -   SUM OF FRACTIONAL RADII
            -   .....
 		    -   VALUES: must be greater than 0
           
 kratio     -   RATIO OF RADII
            -   .....
 		    -   VALUES: must be greater than 0
 
 
 OPTIONAL USER SET:
 
 incl       -   INCLINATION
            -   VALID RANGE = 0-90
 
 e_cosw     -   e*cos(omega) and e*sin(omega)
 e_sinw     -   Used for eccentric orbit, where e is the eccentricity and omega
                 is the longitude of periastron.
            -   DEFAULT = 0 for both, representing a circular orbit
            -   VALID RANGE = 0-1 for both

 
 OTHER OPTIONALS:
            
 u_p        -   LINEAR LIMB DARKENING COEFFICENTS
 u_s        -   u_p is the coefficient for the primary, and u_s is coefficient
                 for the secondary.
            -   DEFAULT = 0.50 for both
            -   VALID RANGE = 0-1 for both
                              
 ygrav_p    -   GRAVITY DARKENING EXPONENTS
 ygrav_s    -   ygrav_p is the exponent for the primary, and ygrav_s is the
                 exponent for the secondary
            -   DEFAULT = 0.08
 
 srefl_p    -   REFLECTED LIGHT FACTORS
 srefl_s    -   srefl_p is factor for the primary and srefl_s is the factor for 
                 the secondary
            -   DEFAULT = 0
           
 qratio     -   MASS RATIO
            -   DEFAULT = 0.5
 		    -   VALUES: must be greater than 0
           
 tangle     -   TIDE LEAD/LAG ANGLE
            -   ......
            -   DEFAULT = 0
 		    -   VALUES: must be greater than or equal to 0
           
 thirdlight -   THIRD LIGHT
            -   DEFAULT = 0
 		    -   VALUES: must be greater than or equal to 0
 
 phasecor   -   PHASE CORRECTION
            -   ............
            -   DEFAULT = 0
			
 magzero    -   LIGHT SCALE FACTOR
            -   ..........
            -   DEFAULT = 0
            
 dgamma     -   INTERGRATION RING SIZE (degrees)
            -   The size of the rings used when calculating the limb dark.
                 and brightness etc.
            -   DEFAULT = 1
 		    -   VALUES: must be greater than 0
 
 
 PHASE OPTIONS:
 
                PHASE will be displayed from 'phasemin' to 'phasemax'          
 
 phasemin   -   Starting point for the phase plot
            -   DEFAULT = -0.2
 
 phasemax   -   End point for the phase plot
            -   DEFAULT = 0.8
 
 nph        -   No. of bins used to split up the phase.
            -   DEFAULT = 1001

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""

# Set out all of the parameters that will be used with their default values if
#  necessary (#default). (some will be test values at the miunte, marked with '#test')

#User set
#sbratio = 0.4307 #test
#rsum = 0.0947 #test
#kratio = 0.3927 #test

#optional user set
incl = 89.5# default
e_cosw = 0 # 0 #default
e_sinw = 0 # 0 #default

# Addtional options
u_p = 0.5 # default
u_s = 0.5 # deafult

y_grav_p = 0.08 #default = 0.08
y_grav_s = 0.08 #default = 0.08

srefl_p = 0 # default = 0
srefl_s = 0 # default = 0

qratio = 0.5 # default = 0.5
tangle = 0 #default = 0
thirdlight = 0 # default = 0

phasecor = 0 #default = 0
magzero = 0 #default = 0

dgamma = 5 #default, 5 for ground based, 1 for space based.

#Phase options
phasemin = -0.2 #default
phasemax = 0.8 #default
nph = 1001 #default

""" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
import numpy as np
import dlight as l2
import matplotlib.pyplot as plt
""" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

def checkUser_YN_Input(userInput):
	""" 
	
	This is a simple function to check whether or not the user has entered a 
	 sensible value for the answer to yes/no style question. It will accept
	 values of "y", "Y", "yes", "Yes" and "n", "N", "no", "No". Anything else
	 and the user will be asked to enter another value, and the code will not
	 continue until a valid answer is entered.
	 
	PARAMETERS:
	 
	 userInput  -  The string that needs to be tested
	  
	 
	RETURN:
	 
	 userInput  -  If the input was valid it will return the original 
	                 userInput string, however if the user had to enter a new
	                 input, then this will be returned instead.
	
	"""
	
	# A boolean type variable which tracks whether or not the user has entered 
	#  a sensible value
	suitableInput = False
	# Loop to contiune checking the users input, until a sensible value is
	#  entered
	while suitableInput ==False:
		# Acceptable "yes" values
		if (userInput == "y") or (userInput == "Y") or \
		        (userInput == "yes") or (userInput == "Yes"):
			
			suitableInput = True
		
		# Acceptable "no" values
		elif (userInput == "n") or (userInput == "N") or \
		        (userInput == "no") or (userInput == "No"):
			
			suitableInput = True
		
		# If anything thing else is entered
		else: 
			userInput = raw_input("INVALID INPUT - Please try again: ")
	
	return userInput
	
def chk_inputLength(arr):
	""" 
	To check that the two values have been entered, more than 2 ask for new input,
	 if 1 assume it is an observation and take error to be 1. has to be greater 
	 than zero.
	 
	 	PARAMETERS:
	 		
	 		arr: the array that needs to be checked
	 		
	 		
	 	RETURN:
	 	
	 		arr: Will the be the original array if the length is ok, or other wise
	 		 will have been corrected.
	 
	"""
	
	suitableFormat = False
	while suitableFormat == False:
	
		length = len(arr)
				
		if length == 2 or length == 1:
			suitableFormat = True
			
			if length == 1:
				arr = np.append(arr, 1.0)
				print 'Single Value warning: Taken to be observation, error = 1.0'

		if length > 2 or length <= 0:
			line = raw_input('ERROR - Incorrect length format - Please try again: ')
			arr = np.array(map(LIGHTfuncInteraction.float_chk, line.split()))


	return arr
	


def float_chk(string):
	""" 
	
	This purpose of this function is the test whether or not a string inputted
	 by the user can be successfully converted into a float. It attemps to 
	 convert 'string' into a float, if successful the converted value is return.
	 If unsuccessful, it will produce a ValueError, and request a new input from
	 the user. This request will be repeat until a value that can be converted 
	 to a float is entered.
	
	PARAMETERS:
	
	 string  -  The string that you want to turn into a float
	 
	RETURN:
	 
	 convertedVal  -  The successfully converted float version of the string
	
	"""
	
	# Variable that controls whether or not the conversion has been succesful,
	suitableInput = False
	
	# A loop to ensure a reasonable value is eventually entered
	while suitableInput==False:
		# Try to convert the value to a float, if successful, return converted
		#  value and allow the loop to be exited.
		try:
			convertedVal = float(string)
			suitableInput = True
			return convertedVal
		# If unsuccessful, flag an error and ask user for a new input.
		except ValueError:
			print string, "has incorrect format." 
			string = raw_input("Please enter a float: ")
	
def YN_to_bool(user_ans):
	""" 
	
	Takes a yes/no answer from the user an converts it into a boolean true/false
	 type. This is more flexiable than working with individual strings of chars
	 that could be different depend on what the user entered.
	 
	 PARAMETERS:
	 	
	 	user_ans - the reponse the user will have given to a yes/no question
	 	
	 RETURN:
	 	
	 	bool - the true/false answer the yes/no reponse has been converted into
	 	
	"""
	# If user answer is yes, etc. bool it true
	if (user_ans == "y") or (user_ans == "Y") or \
	        (user_ans == "yes") or (user_ans == "Yes"):
			
		bool = True
		
	# If used answer is no, etc. bool is false
	else:
		bool = False
	
	return bool

def chk_val_in_range(testVal, lowerLimit=0, upperLimit=None):
	""" 
	
	This function will check to see if the number 'testVal' is outside of the
	 limits 'lowerLimit' and 'upperLimit'. testVal can be equal to these values
	 If it is outside the range, will ask user for a new value. It will not
	 continue until a suitable value is entered.
	 
	 PARAMETERS:
	 	
	 	testVal - The value that needs to be checked
	 	lowerLimit - The lowest number in the range that is still valid
	 				  (Default = 0)
	 	upperLimit - The highest number in the range that is still valid
	 				  (Default = infinity)
	 	
	 RETURN:
	 
	 	testVal - Will return the original value if value was with the range,
	 			   if not it will return the new value.
	 			   
	"""
	
	suitableVal = False
	while suitableVal == False:
	
		if (testVal < lowerLimit 
			or (upperLimit is not None and testVal>upperLimit)):
			
			if upperLimit == None:
				upperLimit = "Infinity"
			
			print "This value is outside the range " + str(lowerLimit) + "-" + \
				str(upperLimit) + "."
			testVal = float_chk(raw_input("Please enter a value within this"
						" range: "))
		else:
			suitableVal = True
	
	return testVal
	

def bigger_than_chk(testVal, lowerLim=0):
	""" 
	
	This function will check whether the number 'testVal' is bigger than the 
	 value 'lowerLim', which is set to zero by defalut. If the value doesn't 
	 pass the test, then user is requested to enter a new value. NOTE: the
	 test value cannot equal 'lowerLim'
	 
	 PARAMETERS:
	 
	 	testVal   -  The number to be checked
	 	
	 	lowerLim  -  'testVal' must be bigger than this number.
	 	
	 RETURN:
	 
	 	testVal   -  Returns the orginal value if testVal is larger than LowerLim
	 				   otherwise it will return the new testVal entered by the
	 				   user
	
	"""
	suitableVal = False
	while suitableVal == False:
		
		if testVal > lowerLim:
			suitableVal = True
		else:
			testVal = float_chk(raw_input("INVALID INPUT - Please enter a value"
						" that is larger than "+ str(lowerLim)+": "))
	
	return testVal


def extreme_chk(testVal, lowerLim, upperLim):
	""" 
	 This function looks at whether or not a value is outside of the limits,
	  'lowerLim' and 'upperLim', and gives a warning if this is found to be
	  true. The user then has the option of changing the value if they wish
	  or the can continue with the original value they entered.
	  
	  PARAMETERS:
	  
	  	testVal   -  The number that is going to be checked
	  	lowerLim  -  The smallest number allowed in the range
	  	upperLim  -  The largest number allowed in the range
	  	
	  RETURN:
	  
	  	testVal   -  If the value was found to be in the allowed range, or
	  					the user wanted to continue with the value they entered
	  					then the original value will be returned. Otherwise
	  					the users new value will be entered.
	"""
	suitableVal = False
	while suitableVal == False:
		#If outside the limits, warn user that this is the case, and ask user 
		# if they want to use the values anyway
		if testVal > upperLim or testVal < lowerLim:
			user_ans = raw_input("This value seems a bit extreme, are you sure"
						" you want to continue? (y/n)")
			user_ans = checkUser_YN_Input(user_ans)
			bool_user_ans = YN_to_bool(user_ans)
			# If user doesn't want to use the value, ask for new value
			if bool_user_ans == False:
				testVal = float_chk(raw_input("Please enter a new value: "))
			else:
				suitableVal=True
			#if user is happy to use the testVal
		else:
			suitableVal=True
	
	return testVal

def calcMag2(v, phase):
	""" 
	Calc. magnitude for each value in phase using the fortran code 'dlight.f'
	
	 PARAMETERS:
	 	
	 	v - The list of current value for the parameters defining the model
	 	phase - the phase for which the magnitude needs to be calculated
	 	
	 RETURN:
	 
	 	mag - the resulting magnitude values
	 	
	"""
	mag = l2.dlight(v,phase)
	
	return mag

	
def generate_values(incl, e_cosw, e_sinw, u_p, u_s, y_grav_p, y_grav_s, srefl_p,
						 srefl_s, qratio, tangle, thirdlight, phasecor, magzero,
						 dgamma, phasemin, phasemax, nph):

	""" 

	 Code to get the user to enter sensible values for the parameters. Parameter 
	  entered as follows:
		-  sbratio
		-  rsum
	 	-  kratio
	 These 3 parameter will need to be set every time the script is run.

	"""

	# raw_input is a string, will want most of the actual value to behave as floats
	sbratio = float_chk(raw_input("Please enter a value for the central surface" 
	           " brightnes ratio: "))
	# Make sure entered number isn't negative or =0	
	sbratio = bigger_than_chk(sbratio)
	
	#Same for rsum and kratio         
	rsum = float_chk(raw_input("Please enter a value for the fractional sum of"
			   " the two stars radii: "))
	rsum = bigger_than_chk(rsum)
			   
	kratio = float_chk(raw_input("Please enter a value for the ratio of the two"
			   " stars' radii: "))
	kratio = bigger_than_chk(kratio)


	""" 
	 Next the script will ask the user whether or not they would like to edit the
	  change the values used for the inclination, e_cosw and esinw
	"""

	# Ask user abou changing the inclination, whislt checking for sensible answer
	user_ans = checkUser_YN_Input(raw_input("Change inclination? (y/n)"
				" (Default = "+str(incl)+" degrees): "))
	# Convert answer to boolean
	bool_user_ans = YN_to_bool(user_ans)
	# If user want to change the value, ask for new input, otherwise use default
	if bool_user_ans == True:
		suitableVal = False
		incl = float_chk(raw_input("Enter a value for the inclination: "))
		# Check the entered value is between 0 and 90
		incl = chk_val_in_range(incl, 0, 90)
			



	# Ask the user about changing ecosw and esinw
	user_ans = checkUser_YN_Input(raw_input("Change ecos(w) or esin(w)?"
				" (y/n) (Defaults = "+str(e_cosw)+" and "+str(e_sinw)+" respctively: "))
	# Convert answer to boolean
	bool_user_ans = YN_to_bool(user_ans)
	if bool_user_ans == True:
		e_cosw = float_chk(raw_input("Enter a value for ecos(w): "))
		e_cosw = chk_val_in_range(e_cosw, -1, 1)
		e_sinw = float_chk(raw_input("Enter a value for esin(w): "))
		e_sinw = chk_val_in_range(e_sinw, -1, 1)

	""" 

	 Now to give the user the option of setting some additional options, such as the
	  limb darkening coefficients, gravity darkening coefficients, etc. It should not
	  be a requirement to set these, but a choice to set only the only that are 
	  necessary.

	"""

	print "\n~~~~~~~~~~~~~~~~~~~"
	print "ADDITIONAL OPTIONS:"
	print "~~~~~~~~~~~~~~~~~~~\n"
	print("Please select the appropriate letter to set a new value for a parameter."
			 " \nThe default settings will be used if no new value is entered."
			 " \nEnter 'x' to continue once done.\n")

	# Design 'table' to hold the selection letter, parameter name and parameter value
	a = ["a","b","c","d","e","f","g","h","i","j","k","l"]
	b = ["Pri. Limb Darkening","Sec. Limb Darkening", "Pri. Gravity Exponent", 
			"Sec. Gravity Exponent", "Pri. Reflected Light", "Sec. Reflected Light",
			"Mass Ratio", "Lag Angle", "Third Light", "Phase Correction",
			"Light Scale Factor", "Int. Ring Size" ]
	c = [u_p,u_s, y_grav_p, y_grav_s, srefl_p, srefl_s, qratio, tangle, thirdlight,
			phasecor, magzero, dgamma]
	# Create an array of other smaller arrays. The smaller arrays contain the value
	#  from the same elemnent no. from the lists 'a', 'b', and 'c'
	params = zip(a, b, c)

	# This is some formatting to get the 'headers', the "-" means justify on left,
	#  the no. is the length and "s" is string. "-"*no. will repeat "-" no times 
	print "%-6s %-25s %s" % ("Letter", "Parameter Name", "Default")
	print "%6s %25s %9s" % ("-" * 6, "-" * 25, "-" * 9)

	for param in params:
		#"g" is the shortest way to represent the float
		print "%-6s %-25s %9.4g" % (param[0], param[1], param[2])


	""" 
	 The code that will decide what happens depending on the value entered by the
	  user.
	"""

	enterX = False
	valSelect = False
	# While the user hasn't decided to contiune with the program..
	while enterX == False:

		# Store user input, will want to make sure at some point a single letter is 
		#  entered
		if valSelect == False:
			selectedLetter = raw_input("\nPlease select a letter: ")
	
		if selectedLetter == "x":
			enterX = True

		# Set Limb Darkening coefficients
		elif selectedLetter == "a":
			u_p = float_chk(raw_input("Enter a value for the primary's linear limb"
					" darkening coefficient: "))
			u_p = chk_val_in_range(u_p, 0, 1)
			valSelect = False
			continue
		elif selectedLetter == "b":
			u_s = float_chk(raw_input("Enter a value for the secondary's linear limb"
					" darkening coefficient: "))
			u_s = chk_val_in_range(u_s, 0, 1)
			valSelect = False
			continue


		# Set gravity darkening exponents
		elif selectedLetter == "c":
			y_grav_p = float_chk(raw_input("Enter a value for the primary's gravity"
					" darkening exponent: "))
			y_grav_p = chk_val_in_range(y_grav_p, 0, 1)
			valSelect = False
			continue
		elif selectedLetter == "d":
			y_grav_s = float_chk(raw_input("Enter a value for the secondary's gravity"
					" darkening exponent: "))
			y_grav_s = chk_val_in_range(y_grav_s, 0, 1)
			valSelect = False
			continue

	
		# Set reflected light factors
		elif selectedLetter == "e":
			srefl_p = float_chk(raw_input("Enter a value for the light reflected"
						" from the primary star: "))
			srefl_p = chk_val_in_range(srefl_p, 0, 1)
			valSelect = False
			continue
		elif selectedLetter == "f":
			srefl_s = float_chk(raw_input("Enter a value for the light reflected"
						" from the secondary star: "))
			srefl_s = chk_val_in_range(srefl_s, 0, 1)
			valSelect = False
			continue

		#Set mass ratio, lag angle and third light
		elif selectedLetter == "g":
			qratio = float_chk(raw_input("Enter a value for the mass ratio: "))
			qratio = bigger_than_chk(qratio)
			valSelect = False
			continue
		elif selectedLetter == "h":
			tangle = float_chk(raw_input("Enter a value for the lag angle: "))
			tangle = chk_val_in_range(tangle)
			valSelect = False
			continue
		elif selectedLetter == "i":
			thirdlight = float_chk(raw_input("Enter a value for the amount of"
							" thirdlight present: "))
			thirdlight = chk_val_in_range(thirdlight)
			valSelect = False
			continue
	
	
		#Set phase correction, magzero and dgamma
		elif selectedLetter == "j":
			phasecor = float_chk(raw_input("Enter a value for the amount of"
							" phase correction: "))
			phasecor = chk_val_in_range(phasecor, -1, 1)
			valSelect = False
			continue
		elif selectedLetter == "k":
			magzero = float_chk(raw_input("Enter a value for the magnitude's"
						"  zeropoint: "))
			magzero = extreme_chk(magzero, -10, 10)
			valSelect = False
			continue
		elif selectedLetter == "l":
			dgamma = float_chk(raw_input("Enter a value for the size of the "
						"integration ring for the limb darkening (in degrees): "))
			dgamma = bigger_than_chk(dgamma)
			valSelect = False
			continue

	
		# If something else is enter that isn't a valid selection, ask for new 
		#  input
		else:
			valSelect = False
			while valSelect == False:
				selectedLetter = raw_input("Invalid selection - Please try again: ")
				if selectedLetter == "x":
					valSelect = True
				for val in a:
					if selectedLetter == val:
						valSelect = True
				continue



	print "\n~~~~~~~~~~~~~~"
	print "PHASE OPTIONS:"
	print "~~~~~~~~~~~~~~\n"

	# Give the option of change the limit from which a phase plot is generated
	print"Currently a plot will be generated with the phase plotted between", \
			phasemin, "and", phasemax
	# Ask the user about changing phasemin and phasemax
	user_ans = checkUser_YN_Input(raw_input("Change these values? (y/n) "))
	# Convert answer to boolean
	bool_user_ans = YN_to_bool(user_ans)


	maxLessThanMin = True
	# If the user does want to change from default values
	if bool_user_ans == True:
		# maxLessThanMin mean the user will be contiunally aske for new values if 
		#  the value entered for the phasemax is less than the phasemin. Once 
		#  addressed the script will continue.
		while maxLessThanMin == True:
			phasemin = float_chk(raw_input("Please enter a value for minimum phase (-1 - 1): "))
			phasemin = chk_val_in_range(phasemin, -1, 1)
			print "The maximum phase will be the minimum +1"
			phasemax = phasemin + 1.0
			maxLessThanMin = False

	# Ask the user about changing the no of bins use for splitting up the phase
	user_ans = checkUser_YN_Input(raw_input("Change the no. of bins (Default = "+
					str(nph)+") (y/n): "))
	# Convert answer to boolean
	bool_user_ans = YN_to_bool(user_ans)
	if bool_user_ans == True:
		nph = float_chk(raw_input("Please enter a value for number of bins (min 50): "))
		# no max value, really so have just used a large number
		nph = chk_val_in_range(nph, 50, 100000000)


	print "\n~~~~~~~~~~~~~~~~~~"
	print "PARAMETER SUMMARY:"
	print "~~~~~~~~~~~~~~~~~~\n"

	# new b and c lists containing all parameters
	a1 = ["Brightness Ratio", "Sum of Frac. Radii", "Ratio of Radii", 
			"Pri. Limb Darkening", "Sec. Limb Darkening", "Inclination",
			"e cos(omega)", "e sin(omega)", "Pri. Gravity Exponent",
			"Sec. Gravity Exponent", "Pri. Reflected Light", 
			"Sec. Reflected Light","Mass Ratio", "Lag Angle",
			"Third Light", "Phase Correction","Light Scale Factor", "Int. Ring Size",
			"Min. Phase", "Max Phase", "No. of Bins"]
	b1 = [sbratio, rsum, kratio, u_p, u_s, incl, e_cosw, e_sinw, y_grav_p, y_grav_s,
			srefl_p, srefl_s, qratio, tangle, thirdlight, phasecor, magzero, dgamma,
			phasemin, phasemax, nph]
		
	# Create an array of other smaller arrays. The smaller arrays contain the value
	#  from the same elemnent no. from the lists 'a', 'b', and 'c'
	params1 = zip(a1, b1)

	# This is some formatting to get the 'headers', the "-" means justify on left,
	#  the no. is the length and "s" is string. "-"*no. will repeat "-" no times 
	print "%-25s %s" % ("Parameter Name", "Value")
	print "%25s %9s" % ("-" * 25, "-" * 9)

	for param in params1:
		#"g" is the shortest way to represent the float
		print "%-25s %-9.4g" % (param[0], param[1])
	print "%25s %9s" % ("-" * 25, "-" * 9)

	return b1

""" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 Use these parameters in the fortan code file 'dlight.f', to
  generate a synthetic lightcurve based on these parameters.

"""

def startRunning(incl=89.5, e_cosw=0.0, e_sinw=0.0, u_p=0.5, u_s=0.5, y_grav_p=0.08, 
			y_rav_s=0.08, srefl_p=0.0, srefl_s=0.0, qratio=0.5, tangle=0,
			thirdlight=0.0, phasecor=0.0,magzero=0.0, dgamma=1, phasemin =-0.2,
			phasemax = 0.8, nph=1001):

	print "\npyebop - Copyright (C) 2017  Jessica Kirkby-Kent\n"
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	print "This program comes with ABSOLUTELY NO WARRANTY."
	print "This is free software, and you are welcome to redistribute it"
	print "under certain conditions; see GNUlicense.txt file for details."
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
	
	print "\nWelcome to pyebop...\n"
	
	b1 = generate_values(incl, e_cosw, e_sinw, u_p, u_s, y_grav_p, y_grav_s, srefl_p,
			srefl_s, qratio, tangle, thirdlight, phasecor, magzero, dgamma,
			phasemin, phasemax, nph)

	#Set up list of parameter in format needed to be read into the fortran code
	v = b1[0:-3]

#	sbratio = b1[0]
#	rsum = b1[1]
#	kratio = b1[2]

	setMinPhase(b1[-3])
	setMaxPhase(b1[-2])
	setnph(b1[-1])

	#Split phase into bins
	phase = genPhase(phasemin,phasemax,nph)

	mag = calcMag2(v, phase)
	
	return v

def genPhase(phasemin, phasemax, nph):
	""" 
	 Use max and min phase to calcualte the phase bins, which depends on how 
	  many bins 'nph' are used. Returns value as an array
	"""
	phase = np.arange(phasemin,phasemax,(phasemax-phasemin)/nph)
	
	return phase

def setMaxPhase(val):
	""" 
	Use to set a new value for phasemax. Val is the new value
	"""
	global phasemax
	phasemax = val
	
def setMinPhase(val):
	""" 
	Use to set a new value for phasemin. Val is the new value
	"""
	global phasemin
	phasemin = val
def setnph(val):
	""" 
	Use to set a new value for nph. Val is the new value
	"""
	global nph
	nph = val
	
