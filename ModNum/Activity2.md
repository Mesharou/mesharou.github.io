# Run an idealized ocean basin II

## Barotropic vorticity equation

We will revisit the results of  [Stommel](https://www.jgula.fr/ModNum/Stommel48.pdf) and [Munk] (https://www.jgula.fr/ModNum/Munk50.pdf) by diagnosing the gyre dynamics through the use of the barotropic vorticity equation. 

The barotropic vorticity equation can be computed online using the option ```DIAGNOSTICS_VRT```. The different terms are presented in [diagnostics\_barotropic\_vorticity.pdf](http://jgula.fr/ModNum/diagnostics_barotropic_vorticity.pdf)

-

## Westward intensification of gyres [Stommel, 1948]
 
   * **Run the BASIN test case from Activity 1 with additional diagnostics and in the Stommel configuration:**
   
       * Copy the BASIN test case from Activity 1 and create a new test case: (for example case2)
       * Add diagnostics by including in the file```cppdefs.h``` the following line (in the section corresponding to your BASIN test case):
      
        ```
        # define DIAGNOSTICS_VRT
        ```
       * Modify the ```croco.in``` to define outputs for the diagnostics (take this one [croco.in](https://www.jgula.fr/ModNum/croco.in.Basin))

      
   	   * Remove the non-linear terms (advection) by undefining the following key in the ```cppdefs.h```:
      
	     ```
	    # undef UV_ADV
	     ```
      
     * Re-compile and rerun the simulation

  * **Plot and analyze the barotropic vorticity budget:**
      * Using your preferred language (python , matlab, julia, etc.) plot together the different terms of the barotropic vorticity budget averaged over the last 2 years of the simulation for the BASIN test case from Activity 1 [available in the ```basin_diags_vrt_avg.nc```  file].
      * What is the first order balance over the interior of the gyres?
      * What about the western boundary current?
   
 ## Impact of the Beta-effect:

  	* Check what happens if you remove the latitudinal variation of the Coriolis parameter (beta-effect) from the previous case, to test the theory of Stommel. To change the value of beta, you need to copy and edit the file ```ana_grid.F``` :
 
		```
		# if defined BASIN
		                     depth=5000.
		                     f0=1.E-4 
		                     beta=0
		```
        
  	* Plot the different terms of the barotropic vorticity budget averaged over the last 2 years of the simulation. Compare them with the previous one.
    


 ## Viscous boundary layer (Munk, 1950)
	  * Put back the beta-effect. But use now a weaker drag and no-slip lateral conditions (in the ```croco.in```)
	  
		
		```
		bottom_drag:     RDRG(m/s),      RDRG2, Zob [m],  Cdb_min, Cdb_max
		                 3.e-4             0.      0.         0.      0.
		gamma2:
		                 -1.
		```
        
   * Plot the different terms of the barotropic vorticity budget averaged over the last 2 years of the simulation. Compare them with the previous one.
	     
 ## Non-linear effects
   * Starting form the previous case, put back the non-linear terms (advection) by defining the following key in the ```cppdefs.h```
   
     ```
     # define UV_ADV
     ```
   * Plot the different terms of the barotropic vorticity budget averaged over the last 2 years of the simulation. Compare them with the previous one.

   

 ## Make it more turbulent
   * Decrease the explicit dissipation in the ```croco.in```

        ```
lateral_visc:   VISC2    [m^2/sec ]
                 100.  0.
tracer_diff2: TNU2       [m^2/sec]
                 100. 0.
        ```
       
   * Edit the file ```param.h``` and increase the number of points:

   
        ```
#if defined BASIN
      parameter (LLm0=120,   MMm0=100,   N=10)
        ```
   * Estimate the largest theoretical barotropic and baroclinic time-steps for the default case [BASIN from Activity 2]. You can find more information [here]( https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.numerics.timestepping.html)    
   * Find the largest possible barotropic and baroclinic time-steps by running the code.
   * plot the different terms of the barotropic vorticity budget averaged over the last 2 years of the simulation

   
## Additional question

* You can similarly check the kinetic energy budget for the gyre by using the option ```DIAGNOSTICS_KE```
