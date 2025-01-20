This solver is based on the OpenFOAM solver fsiFoam, as modified by Hendrik Reese.
Its purpose is to model the first oscillation cycles of a single laser-induced cavitation bubble near a elastic boundary.
The development of this solver was funded by the German Research Foundation (Deutsche Forschungsgemeinschaft, DFG) under contract OH 75/4-1.

To run a simulation, one must first properly install foam-extend-4.0

https://sourceforge.net/projects/foam-extend/
Before compiling OpenFOAM, uncomment the lines 382-392 in the file src/finiteVolume/finiteVolume/fvSchemes/fvSchemes.C.

In the following instructions it is assumed that OpenFOAM was installed on an Ubuntu operating system.
To execute OpenFOAM commands, one must first enter the OpenFOAM environment by 
executing the alias command for the installed OpenFOAM version (usually 'fe40'),
which has to have been sourced after OpenFOAM installation.

1. Download the FluidSolidInteraction package: https://openfoamwiki.net/index.php/Extend-bazaar/Toolkits/Fluid-structure_interaction
2. In the OpenFOAM installation folder, along with a folder called "foam-extend-4.0", there should be another folder named after the PC user name followed by "-4.0". Unpack the package in that folder.
3. Integrate the folder "src" given here in the FluidSolidInteraction folder.
4. Before compiling the package, execute the following command:
sed -i -e 's=\(using namespace\)=#include <vector>\n\1=' fluidSolidInteraction/fluidSolvers/finiteVolume/RBFMeshMotionSolver/RBFMeshMotionSolver.C
5. Open a terminal in "FluidStructureInteraction/src" and compile it including the new solver via './Allwmake'.
6. Extract the folder "exampleCases" in any location.
7. Run the simulation by opening a terminal in the folder "fluid" in any simulation folder and executing './Allrun' for serial computation or './AllrunPar' for parallel computation.
8. If executed in parallel computation, reconstruct the simulation data using './reconstruct'.
9. To view the simulation results, open the file fluid.foam using ParaView.

The progress of the simulation may be monitored by running 'tail -f log.CavBubbleFsiFoam' or 'tail -f info.csv' in the folder "fluid" in the simulation folder.
Simulation settings and parameters may be altered in the text files within the simulation folder (e.g. bubble initial conditions in "constant/transportProperties" or simulation end time and field output interval in "system/controlDict").
The simulation geometry may be altered by changing the file "constant/polyMesh/blockMeshDict".
