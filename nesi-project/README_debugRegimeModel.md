# Debug RegimeModel

Build the command line version (requires CMake):

```
# clone and switch to branch
git clone https://github.com/jtunnicl/GRATE.git
cd GRATE
git checkout xml-input-merge-plot

# create a build directory
mkdir build
cd build

# compile
cmake .. -DDEBUG_REGIME_MODEL=ON
make -j

# run for 4 steps (so goes in regimeModel once)
./GrateRipCLI 4

# there should now be the file plot_regimeModel.csv containing XS.width, XS.Qb_cap values
# this can be loaded in Excel for example, or plotted with the Python script below

# plot using python script (requires Python with numpy, matplotlib)
python ../nesi-project/plot_function.py plot_regimeModel.csv
```
