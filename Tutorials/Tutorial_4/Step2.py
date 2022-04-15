#
#                    Clément Ranc, 14 April 2022
#
import emcee
import glob
import nstarwrap.daotools as daotools
import multiprocessing as mproc
import numpy as np
import os
import pickle
import re
from scipy.optimize import least_squares
import sys
import time

os.environ["OMP_NUM_THREADS"] = "1"  # Required by EMCEE documentation

# ------------------------------------------------------------------------------
# --------- UGLY PART WITH DAOPHOT INITIALISATION AND GLOBAL VARIABLES ---------
# ------------------------------------------------------------------------------
# This part is necessary because the object managing the back and forth with
# fortran is not pickalizable.

# ------ USER GLOBAL VARIABLES ------
# DAOPHOT
image_file_name = 'image.fits.pkl'
psf_model_file_name = 'image.psf'
group_file_name = 'image_4stars.grp'

# DAOPHOT.OPT
watch = 0.0
fitrad = 9.8
e1 = 0.75
e2 = 5.0

# Size of the box
box_xmin, box_xmax, box_ymin, box_ymax = 1080, 1133, 1157, 1203

# Initial positions if new run (will be used to initialize MCMC chains)
x1fg, y1fg = 1112.536, 1183.238
x2fg, y2fg = 1112.368, 1189.931
x3fg, y3fg = 1093.820, 1171.420
x4fg, y4fg = 1108.780, 1167.000
flux_ratio12_fg = 0.657333700306
flux_ratio13_fg = 0.173373372520
flux_ratio14_fg = 0.058225386230

err_factor = 0.122580628

# --------- IMPORT SHARED FORTRAN LIB ---------
nstar_wrapper = daotools.NstarPythonWrapper(image=image_file_name,
    watch_daophot=watch, fitrad_daophot=fitrad, e1_daophot=e1,
    e2_daophot=e2, psf=psf_model_file_name, group=group_file_name,
    fitting_box=[box_xmin, box_xmax, box_ymin, box_ymax], err_factor=err_factor)
ncol, nrow = nstar_wrapper.attach()

# ------------------------------------------------------------------------------
# ------------------ GENERAL PURPOSE FUNCTIONS ---------------------------------
# ------------------------------------------------------------------------------

err_factor_sqrt = np.sqrt(err_factor)

def run_nstar(x):
    '''Run NSTAR to get the chi2 and the total flux for specific centroid positions.
    '''
    chi2, total_flux = nstar_wrapper.nstar_goodness_of_fit(x)
    return chi2, total_flux

def cost_function_res(x):
    chi2, _ = nstar_wrapper.nstar_goodness_of_fit(x)
    return nstar_wrapper.residuals / err_factor_sqrt

def log_prior(x):
    if len(x) == 5:
        if (box_xmin < x[0] < box_xmax)\
            and (box_xmin < x[2] < box_xmax)\
            and (box_ymin < x[1] < box_ymax)\
            and (box_ymin < x[3] < box_ymax)\
            and (0 < x[4]):
            return 0.0

    if len(x) == 8:
        if (box_xmin < x[0] < box_xmax)\
            and (box_xmin < x[2] < box_xmax)\
            and (box_xmin < x[4] < box_xmax)\
            and (box_ymin < x[1] < box_ymax)\
            and (box_ymin < x[3] < box_ymax)\
            and (box_ymin < x[5] < box_ymax)\
            and (0 < x[6])\
            and (0 < x[7] < 1):
            return 0.0

    if len(x) == 11:
        if (box_xmin < x[0] < box_xmax)\
            and (box_xmin < x[2] < box_xmax)\
            and (box_xmin < x[4] < box_xmax)\
            and (box_xmin < x[6] < box_xmax)\
            and (box_ymin < x[1] < box_ymax)\
            and (box_ymin < x[3] < box_ymax)\
            and (box_ymin < x[5] < box_ymax)\
            and (box_ymin < x[7] < box_ymax)\
            and (0 < x[8])\
            and (0 < x[9] < 1)\
            and (0 < x[10] < 1):
            return 0.0

    return -np.inf

def log_prob(x):

    if len(x) == 5:
        fratio = x[4]
    if len(x) == 8:
        fratio = x[6]
    elif len(x) == 11:
        fratio = x[8]

    if fratio > 1:
        fratio = 1.0 / fratio
        x1, x2, y1, y2 = x[2], x[0], x[3], x[1]
        x[0] = x1
        x[1] = y1
        x[2] = x2
        x[3] = y2
    lp = log_prior(x)
    chi2, z0 = run_nstar(x)
    return -0.5 * chi2 + lp, z0

def check_filename(run_id):
    filename = f"{run_id}.h5"
    if os.path.exists(filename):
        print(f"File {filename} already exists.")
        value = input(f"Type a new run ID or enter to overwrite the file: ")
        if value == '':
            return run_id
        else:
            run_id = check_filename(value)
            return run_id
    else:
        return run_id

# ------------------------------------------------------------------------------
# ------------------------------------ MAIN ------------------------------------
# ------------------------------------------------------------------------------
if __name__ == "__main__":

    # User parameters
    flag_perform_levemberg = not True
    flag_perform_mcmc = True
    flag_continue = not True  # Do we continue a previous MCMC run?
    run_id = "4-star_fit"
    max_n = 10000  # Length of each each
    ndim = 11  # Number of variables
    nwalkers = 22  # Number of chains (should be at least 2 x ndim)
    verbose = 1  # 0=code does not ask anything.

    # Random seed for reproductibility and debugging
    SEED = 42

    # Initialize the random generator
    rng = np.random.default_rng(SEED)

    # Log file
    file_log = open(f'{run_id}.txt', 'a')
    file_log.write(f'START OF {run_id}\n')
    file_log.close()

    # First find solution with a least-square algorithm
    if flag_perform_levemberg:
        fnamelsq = f"{run_id}_lsq.pkl"

        if not os.path.exists(fnamelsq):
            p0 = np.array([x1fg, y1fg, x2fg, y2fg, x3fg, y3fg, x4fg, y4fg, 
                flux_ratio12_fg, flux_ratio13_fg, flux_ratio14_fg])

            if verbose > 0:
                print("Fit with Levenberg-Marquardt algorithm in progress...")
            start = time.time()
            lsq_result = least_squares(cost_function_res, p0, method='lm')
            end = time.time()
            run_time = end - start
            if verbose > 0:
                print(f"Run took {0:.1f} minute(s)".format(run_time))
            entries_to_remove = ('fun', 'grad', 'jac')
            for k in entries_to_remove:
                lsq_result.pop(k, None)
            save_file = open(fnamelsq, "wb")
            pickle.dump(lsq_result, save_file)
            save_file.close()
            if verbose > 0:
                print(f"Fit results saved in {fnamelsq}")
                print('Best fit:\n   ', lsq_result.x)
        else:
            save_file = open(fnamelsq, "rb")
            lsq_result = pickle.load(save_file)
            save_file.close()
            if verbose > 0:
                print("Fit with Levenberg-Marquardt algorithm requested.")
                print(f"File {fnamelsq} already exist. No run performed, but I load previous results.")
                print('Best fit:\n   ', lsq_result.x)

    if flag_perform_mcmc:
        filename = f'{run_id}.h5'
        if flag_continue & ~os.path.exists(filename):
            if verbose > 0:
                print(f"Previous run not found. Starting a new one with ID: {run_id}")
            flag_continue = False

        if not flag_continue & (verbose > 0):
            run_id = check_filename(run_id)
            filename = f"{run_id}.h5"

        # MCMC Initialization
        backend = emcee.backends.HDFBackend(filename)
        if flag_continue:
            if verbose > 0:
                print("Initial size: {0}".format(backend.iteration))
            p0 = None
            p0i = backend.get_chain(discard=0, flat=False)[-1][0]
            # If continue a previous run, ignore nwalkers and ndim
            nwalkers, ndim = backend.shape
        else:
            backend.reset(nwalkers, ndim)

            # Initial position of walkers inside fitting box
            star1 = x1fg + 1j * y1fg
            star2 = x2fg + 1j * y2fg
            star3 = x3fg + 1j * y3fg
            star4 = x4fg + 1j * y4fg

            mu, sigma = 0, 1
            inside_box = [True]
            while(any(inside_box)):
                p1 = star1 + rng.normal(mu, sigma, nwalkers) * np.exp(1j * rng.uniform(0, 2*np.pi, nwalkers))
                p2 = star2 + rng.normal(mu, sigma, nwalkers) * np.exp(1j * rng.uniform(0, 2*np.pi, nwalkers))
                p3 = star3 + rng.normal(mu, sigma, nwalkers) * np.exp(1j * rng.uniform(0, 2*np.pi, nwalkers))
                p4 = star4 + rng.normal(mu, sigma, nwalkers) * np.exp(1j * rng.uniform(0, 2*np.pi, nwalkers))

                p0 = np.array([
                    p1.real,
                    p1.imag,
                    p2.real,
                    p2.imag,
                    p3.real,
                    p3.imag,
                    p4.real,
                    p4.imag,
                    flux_ratio12_fg + rng.normal(0, 0.05, nwalkers),
                    flux_ratio13_fg + rng.normal(0, 0.05, nwalkers),
                    flux_ratio14_fg + rng.normal(0, 0.001, nwalkers)
                    ]).T

                p0[0] = np.array([x1fg, y1fg, x2fg, y2fg, x3fg, y3fg, x4fg, y4fg, 
                    flux_ratio12_fg, flux_ratio13_fg, flux_ratio14_fg])
                p0i = p0[0]

                inside_box = np.array([log_prior(a) for a in p0]) < 0


        # Run nstar once to evaluate chi-square of initial position and
        # get the final fitting box and degrees of freedom.
        chi2_0, total_flux_0 = run_nstar(p0i)
        txt = f'Degrees of freedom: {nstar_wrapper.dof:d}\n'
        txt += f'Zeropoint Mag: {nstar_wrapper.zpmag.value:.4f}\n'
        box_xmin, box_xmax = nstar_wrapper.opt['box_xmin'].value, nstar_wrapper.opt['box_xmax'].value
        box_ymin, box_ymax = nstar_wrapper.opt['box_ymin'].value, nstar_wrapper.opt['box_ymax'].value
        box_center = [0.5 * (box_xmax + box_xmin), 0.5 * (box_ymin + box_ymax)]
        box_size = [box_xmax - box_xmin, box_ymax - box_ymin]
        txt += f'Box center: {box_center[0]}   {box_center[1]}\n'
        txt += f'Fitting box size: {box_size[0]}   {box_size[1]} pixels\n'
        if verbose > 0: print(txt)
        file_log = open(f'{run_id}.txt', 'a')
        file_log.write(txt)
        file_log.close()

        # Set up the backend
        filename = f"{run_id}.h5"

        # Define additional non-fitted parameters we want to record
        nstar_params_name_light = ["z0"]
        dtype = [(a, np.float64) for a in nstar_params_name_light]

        # Multiprocessing
        txt = f'Number of available CPU: {mproc.cpu_count()}'
        ncpu = mproc.cpu_count()
        # ncpu = 1

        if ncpu > 1:
            if verbose > 0:
                print(f'Initialization of MCMC on {ncpu} CPU.')
            with mproc.Pool(ncpu) as pool:
                # Initialization of the sampler with all the previous parameters
                sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, blobs_dtype=dtype, backend=backend, pool=pool)

                # Start the MCMC run
                start = time.time()
                sampler.run_mcmc(p0, max_n, progress=True);
                end = time.time()
                run_time = end - start
        else:
            # Initialization of the sampler with all the previous parameters
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, blobs_dtype=dtype, backend=backend)

            # Start the MCMC run
            start = time.time()
            sampler.run_mcmc(p0, max_n, progress=True);
            end = time.time()
            run_time = end - start

        fortfiles = [f for f in os.listdir('./') if re.search(r'fort\.(0|[1-9][0-9]*)$', f)]
        if len(fortfiles) > 0:
            [os.remove(a) for a in fortfiles]

        # Display some basic information
        if verbose > 0:       
            print("Run took {0:.1f} minute(s)".format(run_time/60))
            print("Final size: {0}".format(backend.iteration))
            print(f'Degrees of freedom: {nstar_wrapper.dof:d}')
            print(f'Zeropoint Mag: {nstar_wrapper.zpmag.value:.4f}')




