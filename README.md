# ![gyros](https://img.freepik.com/free-vector/flat-design-nutritious-shawarma-illustration_23-2149009055.jpg?w=1380&t=st=1686171643~exp=1686172243~hmac=b6f974ccb3333dd79036ecacff75c857471d6d2689749324f9f493d924737da8)


Get light curves from NASA's *Kepler*, *K2*, and *TESS* space missions! Measure light curve properties, such as rotation, amplitude, and CDPP.

The main purpose of this package is to automate the downloading of *TESS* Full Frame Image Data and subsequent extraction of CPM lightcurves, but there are some nice functions to automate the downloading of *Kepler* and *K2* lightcurves as well. The *TESS* CPM code is containted in the `ffi_cpm.py` file, and the *Kepler* and *K2* code is containted in `lk_interface.py`.

## Dependencies
This package requires proper installation of the [`lightkurve`](https://docs.lightkurve.org/) and [`unpopular`](https://github.com/soichiro-hattori/unpopular). Make sure to visit those pages and install them properly. This package also makes use of `astropy`, which should come standard in an Anaconda environment. It could be good to make sure `astropy` is up-to-date.

## Download and Install
To use this package, clone the repository with your favorite method. I recommend opening the terminal and using `git clone https://github.com/jlbush23/gyros`. 

Once you've cloned it to your device, change directory into the module's directory with `cd gyros`.

Type `python setup.py install` into the terminal and hit enter to install `gyros` so that it may be called from any Python environment you want to use it in.

## Testing the FFI Download and CPM extraction

Run a test to make sure the *TESS* FFI download and CPM extraction are working properly. 

Before running the test, examine the `gyros_run.py` file in the `scripts` folder. This is an example runscript that can be used for any stellar population of interest.

While in the top `gyros/` directory in your terminal, run `python scripts/gyros_run.py`. 

This will download and extract CPM lightcurves for a sample of 3 stars from the MELANGE-4 stellar association. 

You will see the code begin to run and print its outputs. Let it run to completion.

A new folder named `test/` will appear locally in your `gyros/` directory. In it should be 3 subdirectories named after each star's TIC ID. Within those folders should be `.fits` files containing the CPM lightcurves for each sector of available *TESS* data, along with other information, including an estimate of the rotation period for that sector. Within `test/`, there should also be an output file with a column named `TESSCut_Avail` that is true if the FFI data was successfully downloaded.

## Using this code

Because the package is importable whenever you're in a Python environment, you can use the individual functions in any of the scripts in the `gyros/gyros/` directory.

To use the `gyros_run.py` script to automate the downloading and extraction of *TESS* CPM lightcurves for new populations of stars, you have two main options:
1) Edit the `gyros_run.py` script within the module each time you have a new population of stars, taking care to name the directories where you want your lightcurves and output files saved.
2) Manually create new directories in your computer's directory explorer, then copy the `gyros_run.py` file into that directory, and edit it there. 

Option 2 may be your best bet for lowest maintenance.


