# pyebop.py

 This code will use the fortran based subroutine LIGHT from the well known lightcurve code EBOP (Nelson & Davis, 1972 [1972ApJ...174..617N](http://adsabs.harvard.edu/abs/1972ApJ...174..617N); Popper & Etzel, 1981 [1981AJ.....86..102P](http://adsabs.harvard.edu/abs/1981AJ.....86..102P)) to fit lightcurves for detached eclipsing binary systems. It incorporates an MCMC parameter-space search (through emcee by Foreman-Mackey et al. 2013 [2013PASP..125..306F](http://adsabs.harvard.edu/abs/2013PASP..125..306F)) to look for  multiple solutions to the fit. It can also calculate uncertainties for the lightcurve parameters using a prayer-bead algorithm as described in Southworth (2008) [2008MNRAS.386.1644S](http://adsabs.harvard.edu/abs/2008MNRAS.386.1644S) and based on the algorithm developed by Jenkins et al. (2002) [2002ApJ...564..495J](http://adsabs.harvard.edu/abs/2002ApJ...564..495J).   
  
  This code is very similar to JKTEBOP (Southworth, 2004 [2004MNRAS.351.1277S](http://adsabs.harvard.edu/abs/2004MNRAS.351.1277S)) 
  in what it does, however this code is written in Python and produces the 
  same results (at least for the systems I have tested it on). PYEBOP does not 
  fit period, P, or the time of primary eclipse, t_min. Also, only a linear 
  limb darkening law is available here, I believe JKTEBOP will also fit a 
  quadratic law.  
  
  
  Note this code was written for Python 2.7. Most of the code will tranclate to Python 3
   but I have not found an equivalent of f2py for Python 3.
 
     
   ** IF YOU FIND THIS CODE USEFUL IN YOUR RESEARCH PLEASE CITE:  **  
         
       [Kirkby-Kent et al., 2016, A&A, 591, A124](http://adsabs.harvard.edu/abs/2016A%26A...591A.124K)  
                 
 
###  OTHER REQUIRED PACKAGES:
 
   * numpy  
   * astropy  
   * matplotlib  
   * [emcee*] (http://dfm.io/emcee/current/#)  
   * [corner*] (https://corner.readthedocs.io/en/latest/)  
   * mpfit there is a version included or [here](https://code.google.com/archive/p/astrolibpy/downloads)  

  * Note emcee and corner are only needed for the MCMC runs.
  
###  FILES THAT SHOULD BE INCLUDED:

	* pyebop.py
	* LIGHTfunctionInteration.py
	* ErrorFormulae.py
	* dlight.f
	* mpfit.py
	* GNUlicense.txt
	* README.md
	
## Getting started

	* **Compile the fortran** using:  
		```
		f2py -c -m dlight dlight.f
		```
		
		The python script expects the Fortran code to be called `dlight`. In principle it could be named something else, but you would then need to edit the module name that `LIGHTfuncInteraction.py` looks for.  
		
	* **Check some settings:** These are the very basics needed to get it to run. More detailed options will follow later.  
	
		If you open up the `pyebop.py` in a text editor and look for line 137 with the code:  
		```
		datafile = "/Users/.../data/Wasp_EB32Fast.dat"
		```
		You need to change this so it points to your own lightcurve file. Note the file needs to have the headers 'Phase', 'Magnitude' and 'Mag_Err', but not necessarily in that order. The code assumes it will be a `.fits` file. If this is not the case you can change line 139  
		```
		dataLC  = Table.read(datafile, format="fits")
		```
		so that the format is `ascii` or whatever. It should accept anything that astropy.table accepts.    
			
	* **To run:** Either start python in a terminal where `pyebop.py` is located and then use  
		```
		>>> execfile(“pyebop.py”)
		```
		and follow the instructions  
		
		**OR**  
		
		in a terminal where `pyebop.py` is located, run   
		```
		python pyebop.py
		```
		and follow the instructions.  
		
		Assuming the code has managed to load the lightcurve data and is running, the first few question are there to provide a starting point for fitting the parameters. If you end up running the code lots with the same parameters it is worth keeping a note of the inputted parameters and then copy/pasting them in to the running code.  
		
		Options to run MCMC and prayer-bead error analysis are given after the initial fit.  
		
	* **Priors on ecosw/esinw:**
		Can choose to put priors on both, just ecosw, or just esinw. They should be entered as  
		```obs err```
		with no brackets, commas etc. Enter 'x' when done.  
		

## Troubleshooting

1) **Problem:** If you get the following error when trying to run the code, please follow the instructions below. The mpfit is module is trying to import a scipy code that doesn't exist.  

	```
	Traceback (most recent call last):  
		File "pyebop.py", line 88, in <module>
			import mpfit as mp
		File "/Users/Jessica/pythonScripts/pyebop/betterVersion/mpfit.py", line 413, in <module>
			import scipy.lib.blas
	ImportError: No module named lib.blas
	```
		
	**Solution:** Open up the mpfit.py code. On line 413 change `import scipy.lib.blas` to `import scipy.linalg.blas` and on lines 599 and 600, where it says `= scipy.lib.blas.get_blas_funcs(` change it to `= scipy.linalg.blas.get_blas_funcs(`  
		
		
## Other Parameters

In the script, below the section which imports various python modules, there are a set of variable that you might want to change. Here are some notes on each of them.

* *varyTheseParameters – (list)* This sets which parameters are fitted. Each number is associated with a parameter as follows: 		
	```
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
 	11 - sec. reflected light factor
 	12 - mass ratio
 	14 - Third light
 	15 - phase correction
 	16 - lightscale factor
	```
	Although parameters in the dlight subroutine, 13 (the tide lag angle) and 17 (integration ring size) are not included for fitting. A description on how each parameter affects the lightcurve will be included later.
		
* *noOfShifts – (float)* The number of times the residuals are shifted where carrying out the prayer-bead algorithm, and calculating the uncertainties on the lightcurve parameters. Really depends on how much data you have and how long you're prepared to wait. I normally use about 50 for a data set with about 100 observations, and normally limit it to 500 for really larger data sets of a few tens of thousands.	

* *MCMCthreads – (integer)* The level of parallelisation used in the MCMC. Usually set to 2 or 4.

* *MCMCwalkers – (integer)* The number of emcee `walkers` that are used in the MCMC. These are essentially (not-quite) independent chains used to explore the parameter-space. Using more walkers can reduce the number of steps needed. I normally have 100-200.

* *MCMCsteps – (integer)* The number of steps taken by each emcee walker. Again depends on how complicated the parameter-space is, but about 1000 is normally enough. The number of burn-in step varies with the number of steps. If MCMCsteps > 250, burn-in = 200, if MCMCsteps >150, burn-in = 100, otherwise burn-in = 50.

* *storeMCMCdata – (boolean)* Set to 'True' if you want to store the steps and loglikelihoods generated during an emcee run.

* *storeMCMCdataLoc – (string)* This where MCMC data will be stored if storeMCMCdata is set to True.

* *GenMCMCcheckPlots – (boolean)* When set to 'True', each parameter produces two plots for checking the MCMC runs..
   1) Value vs step for each walker - to look at how the walkers are exploring
   2) A plot showing the running mean averaging over 'RMCsteps' steps.
   
* *RMCsteps – (integer)* How often the running mean in the GenMCMCcheckPlots is sampled. For 1000 MCMCsteps, RMCstep=50 is a reasonable number.  If a large number of steps is used consider increasing RMCsteps.

* *errAdj – (float)* factor to scale the lightcurve uncertainties.


## The Lightcurve Parameters

Up to sixteen parameters can be included in the lightcurve modelling, and passed to the dlight subroutine. Some of these parameters have a larger affect than others. Below is a brief description of each parameter.

* **Central surface brightness**, *J* - The brightness ratio given by J<sub>2<\sub>/J<sub>1<\sub> of the two stars, where J<sub>i<\sub> is the surface brightness of star i at the centre of the disk used to model the star. Affects how deep the eclipses are relative to each other.
			
* **Sum of the radii**, *r*<sub>sum<\sub> - Defined as *r1 + r2*, where r1 and r2 are the fractional radii for the primary and secondary, respectively. The fractional radius of a star is given by its radius R divided by the semi-major axis of the system a, so that ri = Ri /a. Influences the width of the eclipse.
			
Ratio of the radii, k - Defined as r2/r1. Multiplied with the surface brightness ratio, it defines the depth of a total eclipse.
			
Inclination, i - The angle between the orbital plane of the binary system and the observes line of sight. For eclipsing binaries the inclination has to be close to 90° for the eclipses to be visible in the lightcurve. As the inclination moves away from 90° the eclipses will become shallower, an they will move away from a `u' shape and become more `v'-shaped.
			
ecosw, esinw - The two parameters ecosw and esinw are used to account for eccentricity, e, and longitude of periastron angle ω. A negative ecosw shifts the secondary eclipse to a lower phase, while a positive value move it to a larger phase. esinw alters the relative widths of the eclipses; a positive value indicates the secondary is wider than the primary, while a negative value indicate the primary is wider than the secondary.
			 
Primary/secondary linear limb-darkening coefficient, up,s - Limb-darkening is the term used to describe how stars appear brighter at their centres in comparison to the limbs. Temperature increases with depth as you go towards the centre of a star. The limbs appear dimmer because the line-of-sight does not though these higher temperature regions. Mainly affects the shape of the contact points.
			
Primary/secondary gravity-darkening exponent, βp,s - For bolometric flux, the gravity-darkening exponent, βbol, determines how the local flux, F, or temperature, T, of a star is scaled in proportion to the surface gravity, g, through F α T4 α gβ  (Claret 1998, [1998A&AS..131..395C]). Affects the out-of-eclipse regions of a lightcurve.
			
Primary/secondary reflection coefficient - A brightening effect caused by radiant energy from one star heating the side of the companion star (or vice versa) and increasing the temperature in that region of the star. Can affect the lightcurve in two ways. Firstly, the brightness observed either side of the eclipses is increased, when the more luminous region is visible. Secondly, a sine wave variation in brightness across the lightcurve, meaning the is a substantial magnitude difference in the out-of-eclipse brightness at the start of each eclipse.

Mass ratio, q - Defined as the ratio between the masses of the two stars in the binary, M2/M1. Mainly affects the out-of-eclipse regions, and is only a small effect.
			
Third-light, l3 - The amount of light from sources outside the binary, e.g. a nearby background/foreground star, or a additional stars in the system. It will dilute the eclipses, in a similar way to inclination.
			
Light scale factor, l0 - This is a zero-point parameter that adjusts the position of the normalisation. It does not affect the overall shape of the lightcurve.

Phase correction, φc - Allows the observations to be shifted along the x-axis, in phase, to allow the deeper primary eclipse to sit at a phase of zero.
	# pyebop
