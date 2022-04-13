import nstarwrap.daotools as daotools

if __name__ == "__main__":
    # We run nstarwrap to create a copy of the FITS file, 
    # optimized for (multiprocesses and) MCMCs.
    image_file_name = 'image.fits'
    nstar_wrapper = daotools.NstarPythonWrapper(image=image_file_name, pickle=True)
    ncol, nrow = nstar_wrapper.attach()
