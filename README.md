# GRATE
1-D Simulation of Gravel Routing and Textural Evolution in Rivers

Copyright (C) 2012  Jon Tunnicliffe
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

GRATE has been in development since 2006 - it has had many incarnations for different applications in gravel bed sediment routing problems. The code provided here is built on the Qt Open Source framework; you will need Qt Creator and associated libraries to compile this code, but it should be relatively easy to transfer the core C++ elements to other compilers. While I have worked on development of this particular strand of the code, I acknowledge my collaborators Jeremy Walsh and Murray Hicks from NIWA Christchurch, who have been key developers and mentors in getting this project off the ground.

The code is inspired by Gary Parker's classic 1D sediment transport formulations (e.g. ACRONYM; see Gary's Morphodynamics Web Page: http://hydrolab.illinois.edu/people/parkerg/), though it employs many routines that have been developed in the intervening 25 years or so, including cross-section and floodplain specification, tributary additions, non-erodible sections, gravel abrasion, multiple lithologies, and optimisation of channel width. The abrasion and lithology routines are based on Cui et al. (e.g. TUGS, DREAM). The width optimisation model is based on algorithms developed by Brett Eaton (https://blogs.ubc.ca/beaton/working-code/). For calculating the water energy profile, one has the choice of employing a standard explicit backwater formulation based on Parker's model, or an implicit formulation that is based on Chaudhry (2008).

Application of the code can be found in Tunnicliffe and Church, 2015 (doi:10.1002/2014JF003370) and Tunnicliffe, Eaton, Fuller, Marden and Peacock (2018; forthcoming).

The code has been developed using the Qt open source IDE (https://www1.qt.io/download/). While the key algorithms are written in fairly generic C++, the 'main' header, 'mainwindow.cpp' and 'RwaveWin.ui' elements link operation to a GUI-based framework, allowing for easy model execution and observation of the evolving river profile. The latest version of the code (Nov 2017) has been compiled with Qt 5.9 (QtCreator 4.4). If you have the Qt framework installed, you can open 'GrateRip.pro' as a project, and compile and run from QtCreator.

The code consists essentially of four modules - MainWindow, RiverProfile, Hydro and Sed. MainWindow sets up the Qt graphical interface and handles user input, including start of the simulation. RiverProfile builds the data object that contains all of the cross-section and river profile information. It reads in the model setup information from a .dat file. The Hydro and Sed modules handle all of the routines that relate to hydraulics and sediment transport, respectively.

As development work proceeds, more documentation will be provided here.

## Building the code

### QtCreator (Windows, Linux, Mac)

Load either the `GrateRip.pro` or `CMakeLists.txt` as a project in QtCreator.

### Command line (Linux, Mac)

Using CMake:

```
mkdir build && cd build  # create a directory for the build and change to it
cmake .. -DBUILD_CLI=ON  # configure the build (we enable the CLI version since it is disabled by default)
make                     # build the code
ctest                    # test the code
```

This will build two executables; `GrateRip` is the GUI version, `GrateRipCLI` is the command line version. The same input files work with both.

The CLI/headless version has no dependency on Qt. To build without Qt, disable the GUI version and enable the CLI version at the configure step:

```
...
cmake .. -DBUILD_CLI=ON -DBUILD_GUI=OFF
...
```

This will build only the `GrateRipCLI` command line version.

Different compilers can be specified at the configure step by setting the `CXX` environment variable. For example, if you have the Intel C++ compiler (`icpc`) installed and wish to use it:

```
...
CXX=icpc cmake .. [OTHER_OPTIONS]
...
```

### CMake (Windows)

It should be possible to build the CMake version on Windows too.
