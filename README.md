# gyros
Get light curves from NASA's *Kepler*, *K2*, and *TESS* space missions! Measure light curve properties, such as rotation, amplitude, and CDPP.

The main purpose of this package is to automate the downloading of *TESS* Full Frame Image Data and subsequent extraction of CPM lightcurves, but there are some nice functions to automate the downloading of *Kepler* and *K2* lightcurves as well. The *TESS* CPM code is containted in the `ffi_cpm.py` file, and the *Kepler* and *K2* code is containted in `lk_interface.py`.

## Dependencies
This package requires proper installation of the [`lightkurve`](https://docs.lightkurve.org/) and [`unpopular`](https://github.com/soichiro-hattori/unpopular/]). Make sure to visit those pages and install them properly. This package also makes use of `astropy`, which should come standard in an Anaconda environment. It could be good to make sure `astropy` is up-to-date.

## Download and Install
To use this package, clone the repository with your favorite method. I recommend opening the terminal and using `git clone https://github.com/jlbush23/gyros`. 

Once you've cloned it to your device, change directory into the module's directory with `cd gyros`.

Type `python setup.py install` into the terminal and hit enter to install `gyros` so that it may be called from any Python environment you want to use it in.

## Testing the FFI Download and CPM extraction

Run a test to make sure the *TESS* FFI download and CPM extraction are working properly. 

Before running the test, examine the `gyros_run.py` file in the `scripts` folder. This is an example runscript that can be used for any stellar population of interest.

While in the top `'gyros'` directory in your terminal, run `python scripts/gyros_run.py`. 

This will download and extract CPM lightcurves for a sample of 10 stars from the MELANGE-4 stellar association. 

You will see the code begin to run and print its outputs. Let it run to completion.

A new folder named `test` will appear locally in your `gyros` directory. In it should be 10 subdirectories named after each star's TIC ID. Within those folders should be `.fits` files containing the CPM lightcurves for each sector of available *TESS* data, along with other information, including an estimate of the rotation period for that sector.


