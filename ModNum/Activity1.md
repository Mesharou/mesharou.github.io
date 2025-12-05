# Activity 1 – Run an idealized ocean gyre 

## Install CROCO
  * Open a terminal (e.g. ```Ctrl+Alt+T```)
  * Create a directory (e.g. ```mkdir ~/ModNum/```) for the project, then go there (e.g. ```cd ~/ModNum/```)
  * ```git clone``` the CROCO code from the official website:

    ```
    git clone https://gitlab.inria.fr/croco-ocean/croco.git
    ```

  * You should see the following files (using ```ls ./croco/```), and the source code (i.e. the FORTRAN files ```*.F```) will be in the folder ```./croco/OCEAN/```

    ![Alt text](./Figure/CROCO_content.png "a title")

## Compile
  
  * Create a folder where you will run the model
	
	```
    mkdir ./case1
    [or mkdir -p ~/ModNum/case1]
	```

  * We need to edit the following files: ```jobcomp```, ```cppdefs.h```, ```param.h```, ```croco.in``` so copy them into the folder you just created:
  * 
    ```
    cp ~/ModNum/croco/OCEAN/jobcomp ~/ModNum/case1/
    cp ~/ModNum/croco/OCEAN/cppdefs.h ~/ModNum/case1/
    cp ~/ModNum/croco/OCEAN/param.h ~/ModNum/case1/
    cp ~/ModNum/croco/TEST_CASES/croco.in.Basin ~/ModNum/case1/croco.in
    ```

 * The model needs a fortran compiler and compatible netcdf libraries
    * With Linux on IUEM computers:
        * First, unload everything:
          ```
          module purge
          ```
        * Then load intel compilers and netcdf library:
          ```
          module load intel/12.1 netcdf/c-4.4.1.1-intel12 netcdf/fortran-4.4.4-intel12
          ```
        * You can verify that modules are loaded with
          ```
          module list
          ```
    * With MacOS:
        * You need to have gcc and netcdf installed (using homebrew or macports):
          ```
          brew install gcc
          brew install netcdf
          ```
        * You may need to edit the jobcomp to specify where to find the netcdf libraries (*nf-config doesn’t work with homebrew version of netcdf-fortran*)
          
          ```
          NETCDFLIB="-L/usr/local/lib -lnetcdf -lnetcdff"
          NETCDFINC="-I/usr/local/include"
          #NETCDFLIB=$(nf-config --flibs)
          #NETCDFINC=-I$(nf-config --includedir)
          ```
 * Go to your case and edit the ```jobcomp``` to specify the location of the source code:
```
SOURCE=~/ModNum/croco/OCEAN
```

 * Edit the ```cppdefs.h``` and choose the predefined test case $\color{red}{Basin}$:

	```
	#define BASIN
	...
	#undef REGIONAL
	```
 * Compile the code:
 
   ```
   ./jobcomp
   ```

## Run the model
 * Run the model:
   ```
       ./croco croco.in
       [or ./croco croco.in &> basin.out &]
   ```
* [If you run the model in backgrounf mode (i.e. ```./croco croco.in &> basin.out &```), you can stop it with the unix command  ```kill -9 PID```, where PID is the process ID associated with croco, which you can get with the ```top``` unix command.]
   
 * Look at model variables:
   ```
   ncdump -h basin_his.nc
   ```
 
 * Look at the output:
   ```
   ncview basin_his.nc
   ```

## Modify the namelist (croco.in)
See https://croco-ocean.gitlabpages.inria.fr/croco_doc/tutos/tutos.08.run.html
 * Modify the namelist (```croco.in```) to:
    * make the model run for 20 years (you can approximate 1 year = 360 days)
    * output of history files every 30 days
    * outputs of variables averaged over 5 years
      
 * Note that you also need to modify the ```cppdefs.h``` to output averaged files:
   ![Alt text](./Figure/Basin_averages.png "a title")

## Run the model using openMP
 * We need to edit the file ```cppdefs.h```
 * Find the part of the file where the BASIN case is defined ( ```#elif defined BASIN``` )
   and change ```# undef OPENMP``` into ```# define OPENMP```
 * Edit ```param.h``` to choose the number of processors (NPP) you want to use
   ![Alt text](./Figure/openmp.png "a title")
 * Recompile the code:
   ```
   ./jobcomp
   ```
 * Define the number of threads to use for your computer (the number should be the same than NPP
defined in the ```param.h``` file). Type in your terminal:
   * in bash:
     ```
     export OMP_NUM_THREADS=8
     ```
    * in csh/tcsh:
      ```
      setenv OMP_NUM_THREADS 8
      ```
 * run the model:
   ```
   ./croco croco.in
   ```

## QUESTIONS
 * About the configuration: 
   * What is the forcing?
   * What is the Coriolis parameter?
   * What is the boundary condition?
   * What is the bottom condition?
   * What is the accuracy of the advection scheme?
 * About the circulation:
    * Is it similar to Stommel’s gyre (http://ido.at.fcen.uba.ar/index_archivos/Stommel_1948.pdf)?
