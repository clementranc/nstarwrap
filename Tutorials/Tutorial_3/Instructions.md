# Perform a 2-star fit using a fully optimized approach

## Introduction

## Step 1: save image in a binary file for fast and multiple access

As shown previously in the [Tutorial 2](https://github.com/clementranc/nstarwrap/blob/main/Tutorials/Tutorial_2/Notebook.ipynb), creating a local binary file of the image gives a lot of flexibility and improve the running time. The only drawback of this approach is the large binary file (typically 2 GB) that is created in your disk. This file can be seen as a temporary file that you can delete at the end of your modeling session if you wish.

Consequently, the first step is to create this binary file by running the `ATTACH` command with the option `pickle=True`. It is exactly what the script called `Step1.py` does. Before running it, make sure that the directory includes:
- an image, here image.fits
- a PSF model produced by DAOPHOT, here image.psf
- a list of stars of interest produced by DAOPHOT Group command, here image.grp.

The script is very simple (file `Step1.py`):
```python
import nstarwrap.daotools as daotools

if __name__ == "__main__":
    # We run nstarwrap to create a copy of the FITS file, 
    # optimized for (multiprocesses and) MCMCs.
    image_file_name = 'image.fits'
    nstar_wrapper = daotools.NstarPythonWrapper(image=image_file_name, pickle=True)
    ncol, nrow = nstar_wrapper.attach()
```

You can run the script from your terminal:
```
$ python Step1.py
```
Ouput:
```
/path-on-your-computer/pywrapper.so

     image...


                                      Picture size:   2290  2290
```
The file `image.fits.pkl` should be created in a couple of seconds. You are done with this step!

## Step 2: perform a 2-star MCMC fit using multiprocessing

Step 2 consists of running an MCMC, potentially using all the CPUs that are available in you computer. The `Step2.py` script is subtantially more complex than the tutorials showed until now, but it is this script that must be used to *efficiently* use DAOPHOT-II NSTAR for your research project. Run this script from the terminal with the command:
```
$ python Step1.py
```

Once the run has finished, I suggest you open a Jupyter Notebook to inspect the results, decide if you want to continue the MCMC or change some parameters etc. The [Notebook of Tuto 3]() is showing how to do that and plot the basic correlation plots. You might prefer to use a regular python script: an example is given in the file `scatter_plot.py`. You can run it with the command:
```
$ python scatter_plot.py
```
A plot `correlations.png` should have been created, similar the one below:



## Step 3: monitor and analyse the results









