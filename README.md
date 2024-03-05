# Overview
The repository contains Python files for the manuscript "Viscoelastic mechanics of tidally induced lake drainage in the Amery grounding zone", which is submitted to The Cryosphere. Below is a schematic of the numerical model.

![Image text](https://github.com/HwenZhang/TidalHydroFrac/blob/147148f5916b7197c94a07abe23951a49d448c2f/grounding_line_mesh_sensitivity/image/schematic.png)

The code is modified from the code for simulating viscous grounding line migration developed by [Stubblefield et al., 2021](http://dx.doi.org/10.1017/jfm.2021.394).

# Contents
The repository contains the essential files to simulate the tidal response in the reference case and make the figures in the manuscript.

### Scripts
* **create_mesh.py** is used to produce a mesh with a piecewise linear bottom profile, which will be used for simulations of long-term ice flow when initiating the mesh for tidal simulations.
* **figure_n.ipynb**(n=1~8) are used to make figures in the manuscript. Please note the code needs to be modified if used to visualise new results.
* **sigma_max.ipynb** computes $\sigma_{xx,max}$ on the surface of ice sheets.

### Meshes
Gmsh (https://gmsh.info/) is used for mesh generation in the code.
* **marine_DX12_Lngth20000_Slope2e_2** contains the mesh file for mesh initiation. 

* **tides_DX12_Lngth20000_Slope2e_2_U09.xml** is used to produce the reference case in the manuscript. The mesh file should be loaded with the corresponding initial state **w_init_DX12_L20000_Slope0_02_U09.h5**.

### Source files
The directory contains essential Python files of the program.
* **cases** constains the parameter files for different cases: **params_mesh.py** is for mesh initiation, and **params_tide.py** is for tidal-response simulations. 

* For details of other files, the reader is referred to our manuscript and [Stubblefield et al., 2021](http://dx.doi.org/10.1017/jfm.2021.394).

### Amery_data
The directory contains the water-free, 2-m WorldView-1 (https://earth.esa.int/eogateway/missions/worldview-1) DEM data for the supraglacial lake in the Amery grounding zone, where tidally induced lake drainage was observed.

# How to run
Here we present the tips to initiate the mesh and produce results for the reference case.

## Mesh initiation
Mesh initiation is essentially modelling the low-term marine ice sheet flow to get a steady state, which will be used as the initial condition for tidal-response simulation. The steps are:
* Replace the script **params.py** with **params_mesh.py** from **./source/cases/**.
* Rename **params_mesh.py** as **params.py** and set relevant parameters.
* Run the main file from the parent directory: `python3 ./source/main.py`.
* The procued mesh and initial conditions are all saved in **./results/{casename}**, where **{casename}** is specified in **params.py**.

## Reference case
The reference case represents the tidal response of a viscoelastic marine ice sheet with ocean tides whose amplitude is $1$m. The initial state and mesh are already obtained and attached. The basic workflow is
* Run the FEniCS Docker image. An example command for the author is
`docker run -ti --name fenics-grounding-line-test -w /home/fenics -v /Users/hwenzhang/fenicsproject/grounding_line_test:/home/fenics/shared -d -p 127.0.0.1:8888:8888  quay.io/fenicsproject/stable:latest`.
* Run the main file from the parent directory: `python3 ./source/main.py`.
* The results are saved in the directory with the chosen case name in **results**. The results include a subdirectory **field plot data** that contains the field of velocity, pressure and deviatoric stress,  and a subdirectory **line plot data** that contains the time series, i.e. grounding line positions and surface profile.
* How to postprocess the results is shown below.

## Postprocessing and Visulisation
Postprocessing includes calculating $\sigma_{xx,max}$ on the top surface, calculating stress intensity factor $K_1$ and setting up the model-based criterion. When results are obtained, compute $\sigma_{xx,max}$ with **sigma_max.ipynb**, the make plots using files **figure_n.ipynb**(n=1~8).

The figures in the manuscript are plotted using matplotlib (https://matplotlib.org/stable/), with the Python files given in **./scripts**. To run these codes, make sure that the relevant results are obtained and stored in **./results**, and that the paths in the files are modified.

## Other cases
For cases with a different bedslope angle or rheological parameters, please first run **create_mesh.py** as `python3 ./scripts/create_mesh.py` to create a .xdmf mesh file, and then initiate the mesh by following the above instructions. After getting the mesh and initial state, perform tidal-response simulations. Please adjust the solver parameters accordingly to maintain good solver performance.


# Reference
Stubblefield, A.G., Spiegelman, M., Creyts, T.T., 2021. Variational formulation of marine ice-sheet
and subglacial-lake grounding-line dynamics. Journal of Fluid Mechanics 919. doi:http://dx.doi.org/10.1017/jfm.2021.394.

Trusel, L.D., Pan, Z., Moussavi, M., 2022. Repeated tidally induced hydrofracture of a supraglacial
lake at the Amery ice shelf grounding zone. Geophysical Research Letters 49, e2021GL095661. doi:http://dx.doi.org/10.1029/2021GL095661.


# Acknowledgements
The authors thank Luke Trusel and Anton Fatula for help in the interpretation of their lake-drainage observations, and Ian Hewitt, Ching-Yao Lai, and the RIFT-O-MAT group for discussions on grounding line dynamics and model set-up. We also thank Dave May for valuable suggestions and codes for data visualisation. 

This research received funding from the European Research Council under Horizon 2020 research and innovation program grant agreement number 772255. For more details about the group, the readers are referred to [FOALAB](https://foalab.earth.ox.ac.uk/index.php).
