# GRATE
1-D Simulation of Gravel Routing and Textural Evolution in Rivers

Copyright (C) 2012  Jon Tunnicliffe
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

GRATE has been in development since 2006 - it has had many incarnations for different applications in gravel bed sediment routing problems. The code provided here is built on the Qt Open Source framework; you will need Qt Creator and associated libraries to compile this code, but it should be relatively easy to transfer the core C++ elements to other compilers. While I have worked on development of this particular strand of the code, I acknowledge my collaborators Jeremy Walsh and Murray Hicks from NIWA Christchurch, who have been key developers and mentors in getting this project off the ground.

The code is inspired by Gary Parker's classic 1D sediment transport formulations (e.g. ACRONYM; see Gary's Morphodynamics Web Page: http://hydrolab.illinois.edu/people/parkerg/), though it employs many routines that have been developed in the intervening 25 years or so, including cross-section and floodplain specification, tributary additions, non-erodible sections, gravel abrasion, multiple lithologies, and optimisation of channel width. The abrasion and lithology routines are based on Cui et al. (e.g. TUGS, DREAM). The width optimisation model is based on algorithms developed by Brett Eaton (https://blogs.ubc.ca/beaton/working-code/). For calculating the water energy profile, one has the choice of employing a standard explicit backwater formulation based on Parker's model, or an implicit formulation that is based on Chaudhry (2008).

Application of the code can be found in Tunnicliffe and Church, 2015 (doi:10.1002/2014JF003370) and Tunnicliffe, Eaton, Fuller, Marden and Peacock (2017; forthcoming).

The code consists essentially of four modules - MainWindow, RiverProfile, Hydro and Sed. MainWindow sets up the Qt graphical interface and handles user input, including start of the simulation. RiverProfile builds the data object that contains all of the cross-section and river profile information. It reads in the model setup information from a .dat file. The Hydro and Sed modules handle all of the routines that relate to hydraulics and sediment transport, respectively.

As development work proceeds, more documentation will be provided here.
