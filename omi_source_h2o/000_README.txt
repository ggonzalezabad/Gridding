Gridding software (questions to ggonzalezabad@cfa.harvard.edu)
1. Compiling
      It uses ifort and icc. Other configurations have not been 
      tested but please feel free to do so and let me know the results.
      Basic options:
   	 - Compile tool >gmake target=omi_sao_avg
	 - Clean tool   >gmake clean

2. Inputs
      - Self explained control file: OMI_average_control.inp
      - Orbits to be considered in the gridding: test_list.txt 
        This filename is specified inside OMI_average_control.inp

3. Output
      - Self explained gridded data to OMHCHO_test.he5 (filename
        specified inside OMI_average_control.inp)

4. Running
   >./omi_sao_avg OMI_average_control.inp

5. I'm keeping this software under git at /data/tempo2/ggonzale/Gridding
