This solver is based on the OpenFOAM solver fsiFoam, as modified by Hendrik Reese.
Its purpose is to model the first oscillation cycles of a single laser-induced cavitation bubble near a elastic boundary.

To run a simulation, one must first properly install foam-extend-4.0, as well as the fluidStructure interaction package: 
%link fe40
%link fsi

In the following instructions it is assumed that OpenFOAM was installed on an Ubuntu operating system.
To execute OpenFOAM commands, one must first enter the OpenFOAM environment by 
executing the alias command for the installed OpenFOAM version (usually 'fe40'),
which has to have been sourced after OpenFOAM installation.

1. Integrate the folder "src" in the FluidStructureInteraction installation folder
2. Open a terminal in FluidStructureInteraction and compile the new solver via './Allwmake'.
3. Extract the folder "exampleCases" in any location.
4. Run the simulation by opening a terminal in the folder "fluid" in any simulation folder and executing './Allrun' for serial computation or './AllrunPar' for parallel computation.
5. If executed in parallel computation, reconstruct the simulation data using './reconstruct'.
6. To view the simulation results, open the file fluid.foam using ParaView.

The progress of the simulation may be monitored by running 'tail -f log.CavBubbleFsiFoam' or 'tail -f info.csv' in the folder "fluid" in the simulation folder.
Simulation settings and parameters may be altered in the text files within the simulation folder (e.g. bubble initial conditions in "constant/transportProperties" or simulation end time and field output interval in "system/controlDict").
The simulation geometry may be altered by changing the file "constant/polyMesh/blockMeshDict".
