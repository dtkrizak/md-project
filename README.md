# MD-project

This project is from an Advanced Process Computations course. The basis for the project is given in the file "project guidelines.pdf". 

The fortran file "Krizak_MD_Project.f90" has parameters to set in the file for number of molecule, temperature, sigma, and the cutoff. The parameters have been set so that they apply to any inert gas, such as Argon. The program will first initialize the position and velocity arrays, followed by setting the root mean square velocity, and finally running the simulation. The simulation is run by moving the molecules slightly using the velocities and timestep, followed by recalculating the forces exerted on the molecules, then recalculating the velocities, and finally repeating the procedure until the alloted time is reached.

The fortran program was run and output data was collected and analyzed in the file "Krizak MD Project Data.xlsx". Graphs are available to visualize the velocites which are shown to provide the expected Gaussian distribution.