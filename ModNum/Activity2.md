# Run an idealized ocean basin II
Credits: Jonathan Gula (jonathan.gula@univ-brest.fr)

## Barotropic vorticity equation

  * We will diagnose the barotropic vorticity equation from the simulations. See
    * https://journals.ametsoc.org/view/journals/phoc/31/10/1520-0485_2001_031_2871_wwbcir_2.0.co_2.xml?tab_body=pdf
    * http://jgula.fr/ModNum/diagnostics_croco.pdf

  * Rerun the BASIN test case from Activity 1 with additional diagnostics
      * Add diagnostics in the ```cppdefs.h```
        ```
        # define DIAGNOSTICS_VRT
        ```
      * Modify the ```croco.in``` (you can take this one https://www.jgula.fr/ModNum/croco.in.Basin)
      * Re-compile and rerun the simulation

  * <ins>**QUESTIONS**</ins>:
      * Using your preferred language (python , matlab, julia, etc.) plot together the different terms of the barotropic vorticity budget averaged over the last 5 years of the simulation for the BASIN test case from Activity 1 [available in the ```basin_diags_vrt_avg.nc```  file].
      * What is the first order balance over the interior of the gyres?
      * What about the western boundary current?
   
 ## Westward intensification of gyres (Stommel, 1948)
  * Copy the BASIN test case from Activity 1 and create a new test case: (for example case2)
  * Check what happens if you remove the latitudinal variation of the Coriolis parameter (beta-effect), to test the theory of Stommel. To change the value of beta, you need to copy and edit the file ```ana_grid.F``` :
    ![Alt text](https://github.com/quentinjamet/Tuto/blob/main/Figure/basin_coriolis.png "a title")
  * Plot the different terms of the barotropic vorticity budget averaged over the last 5 years of the simulation. Compare them with the previous one.
    

 ## Viscous boundary layer (Munk, 1950)
   * Use a weaker drag and free-slip lateral conditions (in the ```croco.in```)
     ![Alt text](https://github.com/quentinjamet/Tuto/blob/main/Figure/bottom_drag.png "a title")
   * Plot the different terms of the barotropic vorticity budget averaged over the last 5 years of the simulation. Compare them with the previous one.
     

 ## Non-linear effects
   * Check the impact of the non-linear terms (advection) by removing advection in the ```cppdefs.h```
     ```
     # undef UV_ADV
     ```
   * Plot the different terms of the barotropic vorticity budget averaged over the last 5 years of the simulation. Compare them with the previous one.

 ## Make it more turbulent
   * Decrease the explicit dissipation in the ```croco.in```
     ![Alt text](https://github.com/quentinjamet/Tuto/blob/main/Figure/lateral_dissip.png "a title")
   * Edit the file param.h and increase the number of points:
     ![Alt text](https://github.com/quentinjamet/Tuto/blob/main/Figure/resolution.png "a title")
   * Find the largest possible barotropic and baroclinic time-steps
   * plot the different terms of the barotropic vorticity budget averaged over the last 5 years of the simulation

