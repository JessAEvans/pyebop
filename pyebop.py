""" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 'pyebop.py'
 Copyright (C) 2017 Jessica Kirkby-Kent
 Ver. 1.000 (For python 2.7)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 This code will use the fortran based subroutine LIGHT from the well known
  lightcurve code EBOP (Nelson & Davis, 1972 [1972ApJ...174..617N]; Popper 
  & Etzel, 1981 [1981AJ.....86..102P]) to fit lightcurves for detached eclipsing 
  binary systems. It incorporates an MCMC parameter-space search (through emcee
  by Foreman-Mackey et al. 2013 [2013PASP..125..306F]) to look for  multiple 
  solutions to the fit. It can also calculate uncertainties for the lightcurve 
  parameters using a prayer-bead algorithm as described in Southworth (2008) 
  [2008MNRAS.386.1644S] and based on the algorithm developed by Jenkins 
  et al. (2002) [2002ApJ...564..495J]. 
  
  This code is very similar to JKTEBOP (Southworth, 2004 [2004MNRAS.351.1277S]) 
  in what it does, however this code is written in Python :) and produces the 
  same results (at least for the systems I have tested it on. PYEBOP does not 
  fit period, P, or the time of primary eclipse, t_min. Also, only a linear 
  limb darkening law is available here, I believe JKTEBOP will also fit a 
  quadratic law.
  
  ##############################################################################
     
    IF YOU FIND THIS CODE USEFUL IN YOUR RESEARCH PLEASE CITE:
         
       Kirkby-Kent et al., 2016, A&A, 591, A124
       [http://adsabs.harvard.edu/abs/2016A%26A...591A.124K]
                 
 # #############################################################################

  OTHER REQUIRED PACKAGES:
 
   - numpy
   - astropy
   - matplotlib
   - emcee [http://dfm.io/emcee/current/#]
   - corner [https://corner.readthedocs.io/en/latest/]
   - mpfit (a copy of the version I used is included in the tar file)

  **Note emcee and corner are only needed for the MCMC runs.
 
  This code uses astropy.tables to read in the lightcurve data, and it assumes
   the file contains the headers "Magnitude", "Mag_Err" and "Phase", but not 
   necessarily in that order. It assumes the file is a .fits file, but will 
   accept other file types, e.g. ascii. For this, the line below the datafile
   line i.e.
   
        dataLC  = Table.read(datafile, format="fits")

   needs to be changed. Just alter the format to the appropriate type.
 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 TO SELECT WHICH PARAMETERS WILL VARY IN MINIMISATION:
 	
 	-  Select appropriate numbers, from the following list and set them in
 		'varyTheseParameters'. (13 - lag angle  and 17 - intergration ring size
 		 not included as it's best if these aren't varied)
 		
 		0  - central surface brightness
 		1  - sum of radii
 		2  - ratio of radii
 		3  - prim. limb darkening coefficient
 		4  - sec. limb darkening coefficient
 		5  - inclination
 		6  - e_cos w
 		7  - e_sin w
 		8  - prim. gravity darkening exponent
 		9  - sec. gravity darkening exponent
 		10 - prim. reflected light factor
 		11 - sec. reflectedlight factor
 		12 - mass ratio
 		14 - Third light
 		15 - phase correction
 		16 - lightscale factor
  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
from astropy.table import Table, Column, join
import LIGHTfuncInteraction
from math import pow, atan2, fabs, sqrt, atan, degrees, pi, exp
import numpy as np
import matplotlib.pyplot as plt
import mpfit as mp
import random
import os
import ErrorFormulae as err
""" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
# Choose the parameters that are include in the fit. Parameter details can be 
#  found in the comments at the start of the code
varyTheseParameters = [0,1,2,5,6,7]#[0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16]#


# For prayer-bead error analysis, how times do you want to shift the errors. See
#  documentation for more information. **Note: number must be a float.
noOfShifts = 10.0 #must be float


#      FOR emcee MCMC runs:

# MCMCthreads is the level of parallelisation used in the MCMC. Use "4" for 4-cores,
#  "2" for 2-cores or "1" for no parallelisation.
MCMCthreads = 2

#MCMCwalkers is the number of emcee 'walkers' or chains that ar run.
MCMCwalkers = 20

#MCMCsteps is the number steps in the MCMC chains
MCMCsteps = 100 # **Note minimum number of burn-in steps is 50

# Do you want to store a flatten version of the MCMC position
storeMCMCdata = False
#This is where the positions will be stored if storeMCMCdata is set to 'True'.
storeMCMCdataLoc = "/home/"

# When set to 'True' each parameter produces two plots for checking the MCMC runs
#  1) Value vs step for each walker - to look at how the walkers are exploring
#  2) A plot showing the running mean averaging over 'RMCsteps' steps.
# If a large number of steps is used consider increasing RMCsteps.
GenMCMCcheckPlots = False
RMCsteps = 50

errAdj = 1 #Parameter used to scale lightcurve uncertainties, if required.

""" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
# Read in the table containing the cleaned WASP data. Change which is marked as
#  'datafile' to change which table is used.

datafile = "/Users/Jessica/pythonScripts/pyth3/data/Wasp_EB32Fast.dat"

dataLC  = Table.read(datafile, format="fits")
dataLC.sort(["Phase"])

""" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PARAMETER CLASSES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

def fileExist(directory, name):
	""" 
	Python doesn't like trying to write a 'fits' file to a location  if the file
	 already exists so this checks to see if the file already exists, if it does
	 the file will be deleted, so a new file can be saved without issue.
    
	ARGUMENTS
		directory - the directory path to where the file should be saved
		name      - string containing the name of the actual file
        
    """
	filepath = directory+str(name)
	fileExists = os.path.isfile(filepath)


	while fileExists == True:

		overwriteFile = LIGHTfuncInteraction.checkUser_YN_Input(raw_input("Overwrite existing file? (y/n) "))
		# Convert answer to boolean
		bool_overwrite = LIGHTfuncInteraction.YN_to_bool(overwriteFile)

		if bool_overwrite == True:
			os.remove(filepath)
			fileExists = False

		else:
			new_name = raw_input("Please enter a new file name: ")
			filepath = directory + str(new_name)
			fileExists = os.path.isfile(filepath)
	
	return filepath

def checkElem(varyingElementList):
	""" 
	
	Want to avoid varying 'tangle' and 'dgamma' so if the elements for these 
	 parameters are entered to vary, this will remove them for the list of 
	 varying parameters. If they're not entered then nothing occurs.
	
	"""
	i=0
	while i<len(varyingElementList):
		if varyingElementList[i] == 13:
			varyingElementList = np.delete(varyingElementList, i)
			print "It's best not to vary the lag angle - removing parameter"
		elif varyingElementList[i] == 17:
			varyingElementList = np.delete(varyingElementList, i)
			print "Varying integration ring size will not help find a best fit model" \
						" - Removing parameter"
		elif varyingElementList[i] < 0 or varyingElementList[i] > 17:
			raise InvalidInputError('Please use number between 0 and 17')
		else:
			i+=1
		 
	
	return varyingElementList

def askEC2():

	""" 
	Ask the user to enter a values for ecos(w), in the format 'obs err', this is to
	 check to make sure its a float, and check to make sure it is a senible format
	 to use in the residual calculations. Note it will accept '0' for the error, but
	 it will throw an error later
	
		RETURN:
			ecosw - A 2D array, where the first value in a row in the observational
					 value and the second value is its error.
	"""

	print "\nEnter values for ecosw in format: 'observation' 'error' "
	print "Enter 'x' after final value."
	endLine = 'x'
	ecosw = np.array([])
	for line in iter(raw_input, endLine):
		numArr = np.array(map(LIGHTfuncInteraction.float_chk, line.split()))
		
		ecosw=np.append(ecosw, LIGHTfuncInteraction.chk_inputLength(numArr),axis=0)
	ecosw = np.reshape(ecosw, (-1,2))
	                                
	return ecosw

def askES2():

	""" 
	Ask the user to enter a values for esin(w), in the format 'obs err', this is to
	 check to make sure its a float, and check to make sure it is a senible format
	 to use in the residual calculations. Note it will accept '0' for the error, but
	 it will throw an error later
	
		RETURN:
			esinw - A 2D array, where the first value in a row in the observational
					 value and the second value is its error.
	"""

	print "\nEnter values for esinw in format: 'observation' 'error' "
	print "Enter 'x' after final value."
	endLine = 'x'
	esinw = np.array([])
	for line in iter(raw_input, endLine):
		numArr = np.array(map(LIGHTfuncInteraction.float_chk, line.split()))
		esinw=np.append(esinw, LIGHTfuncInteraction.chk_inputLength(numArr),axis=0)
		
	esinw = np.reshape(esinw, (-1,2))
	                                
	return esinw


def NOTvaryWhichParams(completeList, onesToVary):
	""" 
	If onesToVary is a list of elements of the parameters that are to be 
	 varied, want to create a list of parameters that won't be varied. It 
	 will simply delete the varying parameter from completeList.
	 
	 PARAMETERS:
	 
	 	completeList - A list of all the parameters used in order for the
	 					fortan light.f
	 	onesToVary   - A list of elements of the parameter that will be 
	 					allowed to vary.
	 RETURN:
	 
	 	
	 	elemList - a list of the parameters that are not to be varied for the 
	 			optimization/least square fit part.
	 			
	"""
	
	endRange=len(completeList)
	elemList = range(0,endRange)

	elemList = np.delete(elemList, onesToVary)
		
	return elemList

def lookupBounds(elementList, returnLists=False):
	""" 
	This function will take a list of elements between 0 and 17, and will lookup
	 this element in the Max and Min bounds arrays to find the maximum and 
	 minimum allowed values for the corresponding parameter. The element 
	 representing each parameter is displayed below. Will return an array with 
	 the bounds for the selected params.
	 
	 PARAMETERS:
	 	
	 	elementList	-  an array of elements that will be looked up
	 	returnLists -  If True will alter what is returned 
	 						(the max,min list instead)
	 	
	 RETURN:
	 
	 	paramsBounds - an array of bound for the parameters selected by 
	 					elementList.
	 					Has form [(min[0], max[0]), (min[1],max[1]),....]
	 					
	 	minVal,maxVal- arrays containing the lower and upper boound values for 
	 					each of the parameters selected by the elements in 
	 					elementList
	
	
	FOR REFERENCE:
		sbratio = AllBoundsList[0] 		0-None
		rsum = AllBoundsList[1]			0-None
		kratio = AllBoundsList[2]		0-None
		u_p = AllBoundsList[3]			0-1
		u_s  = AllBoundsList[4]			0-1
		incl = AllBoundsList[5]			0-90
		e_cosw = AllBoundsList[6]		-1-1
		e_sinw = AllBoundsList[7]		-1-1
		y_grav_p = AllBoundsList[8]		0-1	
		y_grav_s = AllBoundsList[9]		0-1	
		srefl_p = AllBoundsList[10]		0-1	
		srefl_s = AllBoundsList[11]		0-1
		qratio = AllBoundsList[12]		0-None
		tangle = AllBoundsList[13]		NOT VARYING(0-360)
		thirdlight = AllBoundsList[14]	0-10.0
		phasecor = AllBoundsList[15]	-1.0-1.0
		magzero = AllBoundsList[16]		-10-10
		dgamma = AllBoundsList[17]		NOT VARYING(1,10)
	"""
	
	# Upper bounds
	MaxBounds = ['None', 'None', 'None', 1.0, 1.0, 90.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
					'None', 360.0, 10.0, 1.0, 10.0, 10.0]	
	# Lower bounds
	MinBounds = [0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0 ,0.0, 0.0, 0.0,
					0.0, 0.0, -1.0, -10.0, 1.0]
	# Pick out the pairs of bounds that are needed for the elements of the parameters
	#  in elementList
	maxVal = np.take(MaxBounds, elementList)
	minVal = np.take(MinBounds, elementList)
	
	# Return the two min and max arrays. Used for the MyBounds class
	if returnLists==True:
		return [minVal, maxVal]
	# A continuous list of upper and lower bounds, used for the bound for 
	#  minimize part
	else:
		paramsBounds = zip(minVal, maxVal)
		return paramsBounds

def lookupParamNames(elementList, completeNames=False):
	""" 
	Function that will take a list of elements and look them up amonugst the 
	 array of parameter names and return the approppriate one as a list.
	 
	When completeNames==True the full parameter name will be returned, otherwise
	 it's shorter name will be used

	FOR REFERENCE:
		sbratio = completeList[0]
		rsum = completeList[1]
		kratio = completeList[2]
		u_p = completeList[3]
		u_s  = completeList[4]
		incl = completeList[5]
		e_cosw = completeList[6]
		e_sinw = completeList[7]
		y_grav_p = completeList[8]
		y_grav_s = completeList[9]
		srefl_p = completeList[10]
		srefl_s = completeList[11]
		qratio = completeList[12]
		tangle = completeList[13]
		thirdlight = completeList[14]
		phasecor = completeList[15]
		magzero = completeList[16]
		dgamma = completeList[17]
	
	"""
	
	if completeNames == True:
		completeList = ['surface_brightness', 'sum_of_radii', 'ratio_of_radii',
						 'limb_p', 'limb_s', 'inclination','e_cosw', 'e_sinw',
						 'y_grav_p', 'y_grav_s', 'refl_p', 'refl_s',
						  'mass_ratio', 'lag_angle', 'thirdlight', 
						'phase_correction', 'light_scale_factor', 'dgamma']
	else:
		
		completeList = ['sb', 'rsum', 'kratio', 'u_p', 'u_s', 'incl', 'e_cosw',
						'e_sinw', 'y_grav_p', 'y_grav_s', 'srefl_p', 'srefl_s',
						 'qratio',  'tangle', 'thirdlight', 'phasecor', 
						 'magzero', 'dgamma']
	
	selectedNames = np.take(completeList, elementList)
	
	return selectedNames

def findStepSize(num):
	""" 
	Set the stepsize to fixed stepsize initially...
	The pre-defined stepsizes for each of the parameters
		
							If fixed	percentage adjusted
		sbratio 	= ss[0]  - 0.00005	1%
		rsum 		= ss[1]  - 0.00001	1
		kratio		= ss[2]  - 0.00005	1
		u_p			= ss[3]  - 0.001	1
		u_s			= ss[4]  - 0.005	1
		incl		= ss[5]  - 0.01		.5
		e_cosw		= ss[6]  - 0.001	1
		e_sinw		= ss[7]  - 0.001	1
		y_grav_p	= ss[8]  - 0.001	1
		y_grav_s	= ss[9]  - 0.001	1
		srefl_p		= ss[10] - 0.0001	2
		srefl_s		= ss[11] - 0.0001	2
		qratio		= ss[12] - 0.0001	2
		tangle		= ss[13] - 0.01		0 #not varying
		thirdlight	= ss[14] - 0.001	1
		phasecor	= ss[15] - 0.0001	.1
		magzero		= ss[16] - 0.0001	.1
		dgamma		= ss[17] - 0.5		0 #not varying
	"""
	stepArray = [0.00005, 0.00001, 0.00005, 0.001, 0.005, 0.01, 0.001, 0.001, 
					0.001, 0.001, 0.0001, 0.0001, 0.001, 0.01, 0.001, 0.0001, 
					0.0001, 0.5]
	
	return stepArray[num]
	
	
def add2Params(num, vtest, varying=True):
	""" 
	Create a 'PARAMETER' class using value stored in the lists above for
	 name, bounds, and take the value from the array 'v'
	 
	 PARAMETERS:
	 	
	 	num - a number representing which element of the 'v' array, it will
	 			be obtaining values/arguments for.
	 	vtest - The list of parameter values entered by the user
	 	varying - boolean, if true the optimization will vary this parameter
	 				is not it will not.
	 				
	 RETURN:
	 
	 	para - a parameter class contain the appropriate information for the
	 			element 'num' in 'v'.
	 
	"""
	# if fix/'fixed' ==1 the that parameter WILL NOT vary. Can be used if want to
	#  add one parameter, and delare it to be fixed from the very beginning.
	if varying==False:
		fix = 1
	else:
		fix = 0
		
	# Limited if first/second value in this array is set to 1 then the 
	# lower/upper side of the parameter WILL be limited. Also, need to
	# give the 'limits'/bounds array at the same time. Lookup the parameters
	# bounds. To setup the 'limited' array, assume that the parameter has bounds
	# the if either of the upper or lower bounds are set to 'None' then change
	# the value in the limited array to 0, to note that it doesn't have a
	# bound it that particular case
	bounds = lookupBounds(num, True)
	
	#for now will leave the 'mpside' as the default 0 for now, and i won't set a
	# maximum stepsize.
	
	param = {'value': vtest[num], 'fixed':fix, 'limited':[1,1], 'limits':[0.,0.],
				'parname':lookupParamNames(num, True), 'step':findStepSize(num)}
	
	#It only ike assigning one bound at a time, and a bound can't be 'None', so 
	# in that case change it so the part is not limited, but se the other bound.
	if bounds[0] == 'None':
		param['limited'][0] = 0
		param['limits'][1] = bounds[1]	
	
	if bounds[1] == 'None':
		param['limited'][1] = 0
		param['limits'][0] = bounds[0]	
	
	else:
		param['limits'][0] = bounds[0]	
		param['limits'][1] = bounds[1]	

	#There not much point in printing the value for the tangle or dgamma. Have it
	# set that there will not be printed
	if num==13 or num==17:
		param['mpprint'] = 0
		 
	        
	return param


def residual2(p, fjac=None, x=None, data=None, data_err=None, constraints=None):
	""" 
	 'p' contains all the details about the parameters to be used in the
	  fitting. First will need to get the parameter values into a format that
	  can be use in the fortran code. 'x' are the values used to generate
	  the model, 'data' is the magnitudes from the lightcurve, and 'data_err' are
	  their associated errors. 
	  
	  fjac if == None, partial derivative shouldn't be calculated. else the
	   derivatives must be computed. is an array of len(p), where each entry
	   is 1 if that parameter is free and 0 if it is fixed. 
	  
	  need to also include a 'status flag' in this redisuals function and
	   and an optional partial derivatives array
	  
	  Calculates the chisq at each value of 'x'. For the inclination, will add
	   add in an a 'dummy' data point with will effective have a very high 
	   chi-sq if the parameter is venturing too close to the boundary of 90,
	   Otherwise this point's chi-sq is zero.	
	"""
	Ecos = constraints[0]
	Esin = constraints[1]

	#Calc model magnitude values for phases from the actual LC data
	model = LIGHTfuncInteraction.calcMag2(p,x)
		
	# Calc the difference between the actual magnitude and the models magnitude
	diff = data - model
	#error weighted residuals
	chiTerm = diff/(data_err*errAdj)


	if -100.0 not in Ecos:
		ECchiTerm = (p[6] - Ecos[:,0])/Ecos[:,1]
		chiTerm = np.append(chiTerm, ECchiTerm)

	if -100.0 not in Esin:
		ESchiTerm = (p[7] - Esin[:,0])/Esin[:,1]
		chiTerm = np.append(chiTerm, ESchiTerm)

	elif (np.array(val==-100.0 for val in Ecos).all() and np.array(val==-100.0 for val in Esin).all()):
		pass
	
#	# Something to try to keep the inclination, p[5] from venturing over 90 degrees
	#  Nothing happens until it reaches 89.95, after which there is a power law 
	#  which is related to difference between the value and 90. The smaller difference
	#  the larger the term that is added to the chisq value. This is done by effectively
	#  adding an extra data point with this large value. If the inclination value 
	#  isn't a problem then a data point with a residual of 0 is added
	if p[5] >= 90:
		chiTerm = np.append(chiTerm, 1000000000)
	elif p[5] > 89.995:
		InclCloseTo90 = 90.0-p[5]
		extraChi = 10**((0.01 - InclCloseTo90)*50)
		chiTerm = np.append(chiTerm, extraChi)
	else:
		chiTerm = np.append(chiTerm,0)

	#Non-negative status value means MPFIT should continue, negative means stop the 
	# calculation
	status =0
	
	return [status, chiTerm]
		
	
""" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

# Run pyebop script to ask the user for starting parameter values
v1 = LIGHTfuncInteraction.startRunning()


#For mpfit, all the information about the parameters in contained within the
# dictionary which I'm going to call 'params'
# This will be used to set up the list that contains a dictionary for each 
# parameter tangle and dgamma are included so they can be easily referenced 
# when getting values to be sent to the fortran code
params = []
for num in range(0,18):
	params.append(add2Params(num, v1))
	
#set up the array the will be used for models
v=[]
for num in range(0,18):
	parNam = lookupParamNames(num)
	v=np.append(v,params[num]['value'])

#From the user defines elements at the top of the script, check the elements
# are sensible
varyEle=checkElem(varyTheseParameters)
# indices of v, that we don't want to vary
NOTvaryEle = NOTvaryWhichParams(v, varyEle) 

# Generate magnitudes for a model based on these initially inputted parameter
#  values, will be used later, when generating plots
OrigModelMag = LIGHTfuncInteraction.calcMag2(v,dataLC["Phase"])


# For any elements not mentioned in the users list of varying parameters, will
#  want to set their 'vary' argument to False. Also set it so the non-varying
#  parameters aren't print on the screen
for num in NOTvaryEle:
	params[num]['fixed'] = 1
	params[num]['mpprint'] = 0

""" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 CONSTRAIN ECOSW AND ESINW?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
# Ask the user if they would like to use observer values to constrain ecosw and esinw.
#  These booleans are set to True with the question that follow, if the user 
#  answers Yes or Y
bool_constrainBoth = False
bool_constrainEC = False
bool_constrainES = False
EC=np.array([-100.0, -100.0])
ES=np.array([-100.0, -100.0])		

print "\n"

if (6 in varyEle and 7 in varyEle):
	# Ask to constrain ecosw and esinw
	constrainBoth = LIGHTfuncInteraction.checkUser_YN_Input(raw_input("Place constrains on both ecosw and esinw with observations? (y/n) "))
	 #Convert answer to boolean
	bool_constrainBoth = LIGHTfuncInteraction.YN_to_bool(constrainBoth)
 
	if bool_constrainBoth == True:
		
		EC = askEC2()
		ES = askES2()
				
if (6 in varyEle and 7 not in varyEle) or (bool_constrainBoth==False and 6 in varyEle and 7 in varyEle):
	constrainEC = LIGHTfuncInteraction.checkUser_YN_Input(raw_input("Constrain ecosw with observations? (y/n) "))
	 #Convert answer to boolean
	bool_constrainEC = LIGHTfuncInteraction.YN_to_bool(constrainEC)
	if bool_constrainEC == True:
		
		EC = askEC2()

if (7 in varyEle and 6 not in varyEle) or (bool_constrainBoth == False and 6 in varyEle and 7 in varyEle):
	constrainES = LIGHTfuncInteraction.checkUser_YN_Input(raw_input("Constrain esinw with observations? (y/n) "))
	bool_constrainES = LIGHTfuncInteraction.YN_to_bool(constrainES)
	if bool_constrainES == True:
	
		ES = askES2()

#List that contains the constraints for ecosw and esinw
ECEScon = [EC,ES]

""" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SETUP PARAMETERS DICTIONARY & RUN INITIAL FIT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

# Dict. of arguments that will be passed to the minimisation procedure.
funArg = {'x':np.array(dataLC['Phase']), 'data':np.array(dataLC['Magnitude']),
			 'data_err':np.array(dataLC['Mag_Err']), 'constraints':ECEScon}
# second argument in mpfit 'xall' should be an array, containing values for the 
#  starting parameters
print "\nCalculating best-fit parameters...."
result = mp.mpfit(residual2, v, functkw=funArg, parinfo=params, xtol=0.00001, 
			quiet=True)

""" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 FIT RESULTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
overParam = zip(result.params, result.perror)

print "%25s %12s %12s" % ("-" * 25, "-" * 12, "-" * 12)
print "%-25s %-12s %s" % ("Parameter", "Value", "1sig error")
print "%25s %12s %12s" % ("-" * 25, "-" * 12, "-" * 12)

for num in range(0, len(overParam)):
	print "%-25s %12.6g %12.6g" % (lookupParamNames(num, True), overParam[num][0], overParam[num][1])
print "%25s %12s %12s" % ("-" * 25, "-" * 12, "-" * 12)

print "\n**Note: 1-sigma error presented here generated from covariance matrix." 
print "\tIf parameter is fixed or it's value touches a boundary it's error is given as zero.\n"

if (6 in varyEle and 7 in varyEle):
	ecc = err.sqrtxsqysq(result.params[7],result.perror[7],result.params[6],result.perror[6])
	print "Orbital Eccentricity, e:\t\t", ecc[0], "+/-", ecc[1]
	omega = err.atanErr(result.params[7],result.perror[7],result.params[6],result.perror[6])
	if omega[0] <0:
		omega[0] = 360+omega[0]
	print "Periastron longitude, omega:\t\t", omega[0], "+/-", omega[1]

if (1 in varyEle and 2 in varyEle):	
	r1= err.r1Err(result.params[1], result.perror[1],result.params[2],result.perror[2])
	r2= err.r2Err(result.params[1], result.perror[1],result.params[2],result.perror[2])
	print "Fractional radius of primary, r_1:\t", r1[0], "+/-", r1[1]
	print "Fractional radius of secondary, r_2:\t", r2[0], "+/-", r2[1] 


print "\nNo. of iterations:", result.niter
print "Status = ", result.status
if (result.status <= 0): print 'error message: ', result.errmsg


TotalData = len(dataLC)+len(EC)+len(ES)
DOF = TotalData-len(varyEle)
print "Total no of data:\t\t\t", TotalData
print "No. of fitted parameters:\t\t", len(varyEle)
print "No. Degrees of freedom: \t\t", DOF
print "Total chi_squared of fit:\t\t", result.fnorm
print 'Reduced chisq of fit:\t\t\t', result.fnorm/DOF

trialRes = dataLC['Magnitude']-(LIGHTfuncInteraction.calcMag2(result.params, dataLC['Phase']))
LCrms = sqrt(np.sum(trialRes**2)/len(dataLC))

LCChisq = np.sum((trialRes/(dataLC['Mag_Err']*errAdj))**2)
LCredChisq = LCChisq/len(dataLC)

print "\nTotal no of LC data points:\t\t", len(dataLC)
print "rms of the LC residuals (mmag):\t\t", LCrms *1000
print "Reduced chisq of LC data:\t\t", LCredChisq 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
if bool_constrainBoth==True or bool_constrainEC==True or bool_constrainES==True:

	if bool_constrainBoth==True or bool_constrainEC == True:

		ECres = result.params[6] - EC[:,0]
		ECrms = sqrt(np.sum(ECres**2)/len(EC))
		ECchisq = np.sum((ECres/EC[:1])**2)
		ECredChisq = ECchisq/len(EC)
	
		print "rms of ecosw residuals:\t\t\t", ECrms
		print "Reduced chisq of the ecosw values:\t", ECredChisq
	if bool_constrainBoth==True or bool_constrainES == True:

		ESres = result.params[7] - ES[:,0]
		ESrms = sqrt(np.sum(ESres**2)/len(ES))
		ESchisq = np.sum((ESres/ES[:1])**2)
		ESredChisq = ESchisq/len(ES)

		print "rms of esinw residuals:\t\t\t", ESrms
		print "Reduced chisq of the esinw values:\t", ESredChisq
	
	print "\n"
	print "%52s" % ("-" * 52)
	print " Results for observed quantities"
	print "%52s" % ("-" * 52)
	print "%5s %14s %9s %12s %6s" % ("Type", "Measurement", "Error", "Model", "Sigma")

	if bool_constrainBoth==True or bool_constrainEC==True:
		for obsNo in range(0,len(EC)):
			print "%5s %14.6g %9.6g %12.6g %6.3g" % ("EC", EC[obsNo,0], EC[obsNo,1], result.params[6], ((result.params[6]-EC[obsNo,0])/EC[obsNo,1]))
	if bool_constrainBoth==True or bool_constrainES==True:
		for obsNo in range(0,len(ES)):
			print "%5s %14.6g %9.6g %12.6g %6.3g" % ("ES", ES[obsNo,0], ES[obsNo,1], result.params[7], ((result.params[7]-ES[obsNo,0])/ES[obsNo,1]))

	print "%52s" % ("-" * 52)
	print "\n"
	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# Copy of the resulting fitted parameters, for plotting and calculating the 
#  residuals
vFitted = np.copy(result.params)


""" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TO PLOT THE MODELLED CURVE OVER THE ORIGINAL DATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
 
def plotShiftArrayPoints(tableLength, phaseVals, maxPhase, minPhase):
	""" 
	 Need to get the phase data points plotted between the right points, along
	  with the magnitude plots.
	  For point that are above the maximumn phase value the user wants, -1 from 
		 the phase of each and add to the shiftedPaseArray.
	
	"""
	# Find the index elements where the phase is outside the max and min 
	# phase
	i=0
	while (i < tableLength):
		if phaseVals[i] > maxPhase:	
			break
		i+=1
	j=0
	while (j < tableLength):
		if phaseVals[j] > minPhase:	
			break
		j+=1

	return [i,j]
	
def shiftArray(array, maxMinVal):
	""" 
	For a given array, will modulate the position of the values with the array
	 so they are in appropriate positions for a specific phase range. 
	
		PARAMETERS
			array - an array that need to be shifted so that the values
						will be at the appropriate point for the phase
						defined by the value from plotSHiftArrayPoints
			maxMinArray - a two value array, the first contain the array
						element where the array value would be positioned with
						a phase greater that maxPhase, and second value
						is the same but for the minPhase
		RETURN:
			shiftedArray - the resulting shifted array
		
	"""
	i = maxMinVal[0]
	j = maxMinVal[1]

	shiftedArray = array[i:]
	shiftedArray = np.append(shiftedArray, array[:i])
	shiftedArray = np.append(shiftedArray, array[:j])

	return shiftedArray
	
	
ijVals = plotShiftArrayPoints(len(dataLC), dataLC['Phase'], LIGHTfuncInteraction.phasemax,
				 LIGHTfuncInteraction.phasemin)
i= ijVals[0]
j= ijVals[1]

shiftedPhaseArray = []
for val in dataLC['Phase'][i:]:
	shiftedPhaseArray = np.append(shiftedPhaseArray, val-1)
# Then want all the point with phase greater than the minimum phase,
#  these don't need to be shifted, just added to the shiftedPhaseArray
#  with the same phase as they had originally.	

shiftedPhaseArray = np.append(shiftedPhaseArray, dataLC['Phase'][:i])

#Finally add any point with phase less than the minimum plotted phase 
#  to the far end of shiftedPhaseArray, by adding 1 to their value
for val in dataLC['Phase'][:j]:
	shiftedPhaseArray = np.append(shiftedPhaseArray, val+1)


#Need to shift the Magnitude data points so they match the new 
#  shifted phase
shiftedMagArray = shiftArray(dataLC["Magnitude"], ijVals)

# generate a model with the fitted parameters
fittedModel = LIGHTfuncInteraction.calcMag2(vFitted, shiftedPhaseArray)


genPlot = LIGHTfuncInteraction.checkUser_YN_Input(raw_input("\nPlot best-fit model against data? (y/n) "))
# Convert answer to boolean
bool_genPlot = LIGHTfuncInteraction.YN_to_bool(genPlot)

if bool_genPlot == True:

	print "\nGenerating plot........\n"


	#Generate residual plot striaght line
	straightLine = [0]*len(shiftedPhaseArray)
	resid = shiftedMagArray-fittedModel
	
	#Generate a figure showing the model plotted against the data, and underneath
	# a plot of the residual between the two.
	f = plt.figure()
	# Setup the two plot, rowspan=2 means that 'ax1' will be twice the depth of 
	#  ax2
	ax1 = plt.subplot2grid((3,1), (0,0), rowspan=2)
	ax2 = plt.subplot2grid((3,1), (2,0), sharex = ax1)
	
	#Plot details for model/data plot
	ax1.scatter(shiftedPhaseArray, shiftedMagArray, marker="o", s=3, 
					edgecolor='none')
	ax1.plot(shiftedPhaseArray, fittedModel, color='r')
	ax1.axis([LIGHTfuncInteraction.phasemin, LIGHTfuncInteraction.phasemax, max(shiftedMagArray)+0.05,
			 min(shiftedMagArray)-0.05])
	ax1.set_ylabel("Magnitude")

	#Plot details for residuals pot
	ax2.scatter(shiftedPhaseArray, resid, marker="o", s=2, edgecolor='none')
	ax2.plot(shiftedPhaseArray, straightLine, color='r')
	ax2.axis([LIGHTfuncInteraction.phasemin, LIGHTfuncInteraction.phasemax, max(resid)+0.005,
		 min(resid)-0.005])
	ax2.set_xlabel("Phase")
	ax2.set_ylabel("O-C")
	
	# set the space between plots to zero, make it so the ticklabel don't
	# show in ax1
	f.subplots_adjust(hspace=0)
	plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
	f.show()	

""" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RUN MCMC STUFF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

def lnlike(residPara, theta, varyEle):
	""" 
	 For calculating the likelihood for the MCMC.
	 
	  residPara contains all the information need to generated a model.
	  theta is an array with the new test parameters.
	  varyEle - is a list of index, indicating which of the orignal model
	    parameters are varying.
	"""

	#The parameters needed for the residual calc are sent in the list 'residPara'
	para = residPara[0]
	fjacVal = residPara[1]
	x = residPara[2]
	mag = residPara[3]
	err = residPara[4]
	const=residPara[5]
  
  	# Places the test params in the list of model LC parmaeters.
	NewParam = np.copy(residPara[0])
	for no in range(0, len(theta)):
		NewParam[varyEle[no]] = theta[no]


	# Calc. array of residuals for the parameters#
	errcode, residArr = residual2(NewParam, fjacVal, x, mag, err, const)
		
	#If there are constraints on ecosw and esinw
	if -100.0 not in const[0]:
		err = np.append(err, const[0][:,1])

	if -100.0 not in const[1]:
		err = np.append(err, const[1][:,1])
	elif (np.array(val==-100.0 for val in const[0]).all() and np.array(val==-100.0 for val in const[1]).all()):
		pass
   	
   	# calc 1/Err for the errors associated with the data magnitudes
	errSqd = err**2
	invSig = 1.0/errSqd
	#Need to remove the extra term that could be added to to the residual calc
	#  if the inclination gets very close to 90.
	
	if len(residArr) > len(invSig):
		residArr1 = residArr[:-1]**2
	
	#This has come from the fitting example on the emcee web page, but residArr will already
	# be weighted by the errors.
	lnlikeNo = -0.5*(np.sum(residArr1 - np.log(2*pi*invSig)))


	return lnlikeNo
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#~~~~~~~~~~~~~~   
def lnprior(theta, varyEle, PriorInfo):
	""" 
	
	Uses the parameter bounds to constrain the parameters. The bound are collected before 
	   the MCMC is run in the infoForPrior 2d array, and then passed here with the other 
	   arguments. Theta has the length of the number of parameter selected to vary, and 
	   each value is a selected parameter step for the MCMC to test. For every parameter 
	   that is allowed to vary, need to be checked to make sure that its value is sensible.
       
	 """
	testCrit = True
	# check the test value is valid
	for val in range(0,len(theta)):
		if PriorInfo[val][1] > theta[val]:
			testCrit = False
		
		if PriorInfo[val][2] < theta[val]:
			testCrit = False
    
	if testCrit==True:
		return 0.0

	else:
		return -np.inf
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
def lnprob(theta, residPar, varyEle, infoForPrior):
	""" 
	 Combine the likelihood probability from lnlike, with that calculated by the priors in lnprior
	"""
	lp = lnprior(theta, varyEle, infoForPrior)
	if not np.isfinite(lp):
		return -np.inf
	overallProb= lp+ lnlike(residPar, theta, varyEle)
	
	return overallProb
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def getPriorInfo(varyEle, params):
	""" 
	Will select appropiate bound values for the MCMC prior, depending what parameters are selected to
	  vary. If a parameter in params has no upper or lower bounds, then it will be set to +/- 10 as
	  appropriate.
	"""
	infoForPrior = []
	for value in varyEle:
		paramInfoList=[]
		# add name
		paramInfoList.append(params[value]['parname'])
			
		# add prior limits
		if params[value]['limited'][0] == 0:
			paramInfoList.append(-10.0)
		else:
			paramInfoList.append(params[value]['limits'][0])
			
		if params[value]['limited'][1] == 0: 
			paramInfoList.append(10.0)
		else:
			if value ==5:
				paramInfoList.append(90.0)
			else: paramInfoList.append(params[value]['limits'][1])
	
		infoForPrior.append(paramInfoList)
	return infoForPrior	

def stepPlots(MCMCsampler, paramNames):
	""" 
	takes a smapler from the runMCMC and displays the steps for each walker, with each free parameter
		plotted separately
	"""
	noOfParams = len(MCMCsampler.chain[1,1,])
	print 'Number of parameters/plot =',noOfParams
	for paramNo in range(0, noOfParams):
		stepPlot = plt.figure()
		for walkerNo in range(0,len(MCMCsampler.chain[:,:,paramNo])):
			plt.plot(MCMCsampler.chain[:,:,paramNo][walkerNo])
			plt.xlabel('Step no.')
			plt.ylabel(paramNames[paramNo])
		stepPlot.show()

def runningMeanCheck(sample,paramNames, step=20):
	""" 
	Run this on a sample to see if the values converge, using a running mean technique
	"""
	for i in range(len(sample[0])):
		runningMean = []
		runningMeanFig = plt.figure()
		for num in range(0,len(sample), step):
			runningMean = np.append(runningMean, np.mean(sample[:,i][:num]))
		plt.plot(np.arange(0,len(sample),step), runningMean)
		plt.title(paramNames[i])
		runningMeanFig.show()


MCMCrun = LIGHTfuncInteraction.checkUser_YN_Input(raw_input("\nWould you like to run MCMC on the parameters? (y/n) "))
 #Convert answer to boolean
bool_MCMCrun = LIGHTfuncInteraction.YN_to_bool(MCMCrun)
 
if bool_MCMCrun == True:
	print "\n Running MCMC on parameters....\n"
	import emcee

	# No of dimension = no of parameters being fitted
	ndim =len(varyEle)
	
	#Set up a number of starting position equal to the number of MCMCwalkers. Currently set to original
	# starting position with room for a spread around this value.

#	pos = [result.params[varyEle]+result.perror[varyEle]+1e-2*np.random.randn(ndim) for i in range (MCMCwalkers)]
	pos = [result.params[varyEle]+1e-2*np.random.randn(ndim) for i in range (MCMCwalkers)]
 
	# Extra arguments need for the model generation in lnprior part of the MCMC
	resParams = [np.array(result.params), None, np.array(dataLC['Phase']), np.array(dataLC['Magnitude']), np.array(dataLC['Mag_Err']), ECEScon]

	
	#pickput bound values, for parameters used in MCMC. If there is no limit on a parameter
	#  might need to set it to 10 or something
	priorBounds = getPriorInfo(varyEle, params)


 	# Check to make sure the choosen starting points are actually valid.
	print "Checking starting points ...."
	for i in range(len(pos)):
		ValidStart = np.isfinite(lnprior(pos[i], varyEle, priorBounds))
		while ValidStart == False:
			pos[i] = result.params[varyEle]+1e-2*np.random.randn(ndim)
			#pos[i] = result.params[varyEle]+result.perror[varyEle]+1e-2*np.random.randn(ndim)
			ValidStart = np.isfinite(lnprior(pos[i],varyEle, priorBounds))
	
	print "Running sampler..."
	#no of walkers and steps declared at the top of script
	sampler =emcee.EnsembleSampler(MCMCwalkers, ndim, lnprob, threads=MCMCthreads, 
						args = (resParams, varyEle, priorBounds))
	sampler.run_mcmc(pos,MCMCsteps)
	
	if MCMCsteps > 250:
		burnIn = 200
	elif MCMCsteps > 150:
		burnIn = 100
	else:
		burnIn = 50 # note this isn't really enough, but is just here so the MCMC can be run for a few
					#  steps for testing
	print "Using", burnIn, "burn-in steps...\n"
	samples = sampler.chain[:,burnIn:,:].reshape((-1,ndim))
	
	
	#Make pretty param vs param plot using corner.py
	import corner
	labelList = lookupParamNames(varyEle) # FOR AXIS LABELS
	fig = corner.corner(samples, labels=labelList, truths=result.params[varyEle])#, plot_contours=False)
	fig.show()
	if GenMCMCcheckPlots == True:
		stepPlots(sampler, labelList)
		runningMeanCheck(samples, labelList, RMCsteps)
	
	# Store the sampler info in a file, specified by the location 'storeMCMCdataLoc' at the top of 
	#  the script, with the filename 'MCMCoriginal.txt'. If the file already exist there is the 
	#  option to save it as something else.
	if storeMCMCdata == True:
		print 'Saving MCMC walker positons and ln probabilities to file...\n'
		MCMCtable = Table(sampler.flatchain, names = labelList)
		col_prob = Column(sampler.flatlnprobability, name='lnprob')
		MCMCtable.add_column(col_prob)
		
		storedataPath = fileExist(storeMCMCdataLoc, 'MCMCoriginal.txt')
		MCMCtable.write(storedataPath, delimiter = "\t", format='ascii')


	print 'Mean acceptance fraction =', np.mean(sampler.acceptance_fraction)
	print 'Acceptance fraction range (min, max) =', np.min(sampler.acceptance_fraction), np.max(sampler.acceptance_fraction)
	try:
		print "Auto-correlation times for parameters:\n", sampler.acor
	except:
		print "Can't calculate auto correlation time"	


	print "%15s %15s %15s %15s %15s" % ("-" * 15, "-"*15, "-"*15, "-"*15, "-"*15)
	print "%-15s %-15s %-15s %-15s %-15s" % ("Parameter", "50th", "Stan. Dev.","50-15.9", "84.1-50")
	print "%15s %15s %15s %15s %15s" % ("-" * 15, "-"*15, "-"*15, "-"*15, "-"*15)
		
	# Fitted parameter values
	for num in range(0, len(varyEle)):
		per50 = np.percentile(sampler.flatchain[burnIn:, num], 50)
		lowPer = per50 - np.percentile(sampler.flatchain[burnIn:, num], 15.9)
		highPer = np.percentile(sampler.flatchain[burnIn:, num], 84.1) - per50

		print "%-15s %15.6g %15.6g %15.6g %15.6g" % (labelList[num], per50, np.std(sampler.flatchain[burnIn:, num]
			), lowPer, highPer)
	print "%15s %15s %15s %15s %15s" % ("-" * 15, "-"*15, "-"*15, "-"*15, "-"*15)

		
	
	print "\n"
	
	
""" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CALCULATE UNCERTAINTIES THROUGH PRAYER-BEAD ALGORITHM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

calErr = LIGHTfuncInteraction.checkUser_YN_Input(raw_input("\nDo you wish to calculate the uncertainties through the prayer bead algorithm? (y/n) "))
# Convert answer to boolean
bool_calErr = LIGHTfuncInteraction.YN_to_bool(calErr)

if bool_calErr == True:
	
	print '\nCalculating uncertainties....'

	if bool_MCMCrun == True:
		print "...using starting points from the MCMC..."
	
	# Model magnitudes with adjusted parameters
	modelMag = LIGHTfuncInteraction.calcMag2(vFitted,dataLC['Phase'])

	#Calculate the residuals
	res = np.array(dataLC['Magnitude']) - modelMag
	
	#Cyclic shift of residual
	shiftBy = int(round(len(dataLC['Magnitude'])/noOfShifts))

	storeErr=Table(names=lookupParamNames(varyEle))
	# Loop over the number of shifts need to cover entire data set
	totalShifts = 0
	lenNo = len(dataLC['Magnitude'])
	while totalShifts < lenNo:
		

		# Currently the paramter values are reset to the original starting values
		#  before the next shift in residuals. The very first minimisation is the
		# same as the one above.
		noShiftsComplete = totalShifts/shiftBy

		# Shift the residual
		shiftRes = np.roll(res, shiftBy*noShiftsComplete)
		
		# Add shifted residual to model
		synthData = dataLC['Magnitude'] + shiftRes
		
		#Perform 'least sq fit with synthetic data in same way as actual data
		funArg2 = {'x':np.array(dataLC['Phase']), 'data':np.array(synthData), 'data_err':np.array(dataLC['Mag_Err']), 'constraints':ECEScon}
		
		#If MCMC has been run, use random values out of samples. If MCMC hasn't been run just use the
		#  original starting values 'v'
		
		if bool_MCMCrun == True:
			
			MCMCv=np.copy(v)
			for val in range(0,len(varyEle)):
				MCMCv[varyEle[val]] = samples[random.randint(burnIn,len(samples)-1)][val]
			result2 = mp.mpfit(residual2, MCMCv, functkw=funArg2, parinfo=params, xtol=0.00001, quiet=True)
		
		else:
			result2 = mp.mpfit(residual2, v, functkw=funArg2, parinfo=params, xtol=0.00001, quiet=True)
		
		
		#Create an array to store parameter values from a particular minimisation
		#  ready to be added to a table.
		dataForErrTab = []
		dataForErrTab = np.append(dataForErrTab, result2.params[varyEle])

		totalShifts += shiftBy #part of loop

		#print 'totalshift', totalShifts,"/", lenNo, "fun=", result2.fnorm, "iter:", result2.niter
		
		#store each parameter in a table (storeErr)
		storeErr.add_row(dataForErrTab)

	errColData =[]
	for num in range(0, len(varyEle)):
		col = storeErr.colnames[num]
		colStd = np.std(storeErr[col])
		errColData = np.append(errColData, colStd)
	
	storeErr.add_row(errColData)	
	
	# parameters with associated (prayer bead) errors only if the esinw and ecosw are both selected to 
	# to be fitted parameters
	
	if (6 in varyEle and 7 in varyEle):
		ecc = err.sqrtxsqysq(result.params[7],storeErr[params[7]['parname']][-1],result.params[6],storeErr[params[6]['parname']][-1])
		omega = err.atanErr(result.params[7],storeErr[params[7]['parname']][-1],result.params[6],storeErr[params[6]['parname']][-1])
		if omega[0] <0:
			omega[0] = 360 + omega[0]

	#if both rsum and kratio have been selected to be fitted parameters can use them to workout
	#  r_1 and r_2, with associated errors.
	if (1 in varyEle and 2 in varyEle):	
		r1= err.r1Err(result.params[1], storeErr['rsum'][-1],result.params[2],storeErr['kratio'][-1])
		r2= err.r2Err(result.params[1], storeErr['rsum'][-1],result.params[2],storeErr['kratio'][-1])
	

	""" 
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	OVERALL SUMMARY OF PARAMS AND ERRORS
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	"""
	# Make a nice table if Prayer-bead has been rn
	# If the prayer bead analysis has been run, summarise all params and errors here, 
	#  including ecc, omega, r_1 and r_2
		
	# Condense the values and errors from prayer-bead run
	OverallBothParam = zip(result.params[varyEle],storeErr[-1])
		
	#Header
	print "\n\n"
	print "PRAYER-BEAD formatted results..." 
		
	print "%25s %31s" % ("-" * 25, "-"*31)
	print "%-25s %-15s %-15s" % ("Parameter", "Value", "Stan. Dev.")
	print "%25s %31s" % ("-" * 25, "-"*31)
		
	# Fitted parameter values
	for num in range(0, len(OverallBothParam)):
		print "%-25s %15.6g %15.6g" % (lookupParamNames(varyEle[num], True), OverallBothParam[num][0], OverallBothParam[num][1])
		
	#Add in calc. ecc and omega, if appropriate params have been selected
	if (6 in varyEle and 7 in varyEle):
		print "%25s %15s %15s" % ("-" * 25, "-"*15, "-"*15)
		print "%-25s %15.6g %15.6g" %("Orbital Eccentricity", ecc[0], ecc[1])
		print "%-25s %15.9g %15.7g" %("Periastron longitude", omega[0], omega[1])
	#Add in r_1 and r_2 if appropriate params selected
	if (1 in varyEle and 2 in varyEle):
		if (6 not in varyEle or 7 not in varyEle):
			print "%25s %31s" % ("-" * 25, "-"*31)
				
		print "%-25s %15.6g %15.6g" %("Frac. Radius r_1", r1[0], r1[1])
		print "%-25s %15.6g %15.6g" %("Frac. Radius r_2", r2[0], r2[1])
	print "%25s %31s" % ("-" * 25, "-"*31)
		

""" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
print '\n\nDone!\n\n'
print "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "\tIF YOU FIND THIS CODE USEFUL IN YOUR RESEARCH PLEASE CITE:"
print "\n"
print "\t\tKirkby-Kent et al., 2016, A&A, 591, A124"
print "\t  [http://adsabs.harvard.edu/abs/2016A%26A...591A.124K]"
print "\n"
print "\t\t\t\tThanks :)"
print "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
