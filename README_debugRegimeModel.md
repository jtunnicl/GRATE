# Debug RegimeModel

Build the command line version (requires CMake):

# switch to branch
git checkout xml-input-merge-plot

# create a build directory
mkdir build
cd build

# compile
cmake .. -DBUILD_CLI=ON
make -j

# run for 4 steps (so goes in regimeModel once)
./GrateRipCLI 4

# there should now be the file plot_regimeModel.csv containing XS.width, XS.Qb_cap values

# plot using python script (requires Python with numpy, matplotlib)
python ../plot_function.py plot_regimeModel.csv

