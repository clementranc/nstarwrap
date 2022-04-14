#
# Example of script to inspect the correlations between the parameters
#                    Clément Ranc, 14 April 2022
#
import emcee
import moana
import numpy as np
import pandas as pd

if __name__ == "__main__":
    # Constants
    OSIRIS_PIX_SCALE = 9.952  # mas per pixel

    # Select your run
    run_id = '4-star_fit'

    # Select how many MCMC steps you want to throw away (burn-in)
    burnin = 10

    # Load MCMC output
    reader = emcee.backends.HDFBackend(f"{run_id}.h5")
    print(f'Autocorrelations')
    print(f'    Shape: {reader.get_chain(discard=0, flat=False).shape}')
    print(f'   {reader.get_autocorr_time(quiet=True)}')
    samples = reader.get_chain(discard=burnin, flat=True)
    log_prob_samples = reader.get_log_prob(discard=burnin, flat=True)
    nstar_out = reader.get_blobs(discard=burnin, flat=True)

    # Compute missing quantities and miscellaneous things
    colnames = np.array(['x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 
        'flux_ratio', 'flux_ratio13', 'flux_ratio14'])
    mcmc_out = pd.DataFrame(samples, columns=colnames)
    mcmc_out['FTOTAL'] = nstar_out['z0']

    # Don't forget to take into account when flux ratio > 1.0
    # TODO: remove that!
    mask = mcmc_out['flux_ratio'] > 1
    mcmc_out.loc[mask, 'flux_ratio'] = 1 / mcmc_out.loc[mask, 'flux_ratio']
    x1 = mcmc_out.loc[mask, 'x2']
    y1 = mcmc_out.loc[mask, 'y2']
    x2 = mcmc_out.loc[mask, 'x1']
    y2 = mcmc_out.loc[mask, 'y1']
    mcmc_out.loc[mask, 'x1'] = x1
    mcmc_out.loc[mask, 'y1'] = y1
    mcmc_out.loc[mask, 'x2'] = x2
    mcmc_out.loc[mask, 'y2'] = y2
    mcmc_out['chi2'] = - 2.0 * log_prob_samples
    mcmc_out['dchi2'] = mcmc_out['chi2'] - np.min(mcmc_out['chi2'])
    print("\nBest-fit:")
    print(mcmc_out[mcmc_out['chi2']==np.min(mcmc_out['chi2'])].drop_duplicates(inplace=False))

    # Load info from the log file
    with open(f"{run_id}.txt", 'r') as file: logfile = file.read()
    dof = [int(a.split(':')[-1].strip()) for a in logfile.split('\n') 
           if a!='' if 'egrees of freedom:' in a][-1]
    print(f"Reduced chi-square: {np.min(mcmc_out['chi2']) / dof:.6f}")

    # Plot scatter plot
    labels = ['x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 
        'flux_ratio', 'flux_ratio13', 'flux_ratio14']
    label_display = dict({
        'x1': r'x_1', 'x2': r'x_2', 'x3': r'x_3', 'x4': r'x_4', 
        'y1': r'y_1', 'y2': r'y_2', 'y3': r'y_3', 'y4': r'y_4',   
        'flux_ratio': r'\varphi_l/\varphi_s',
        'flux_ratio13': r'\varphi_3/\varphi_s',
        'flux_ratio14': r'\varphi_4/\varphi_s'})

    mask = mcmc_out['dchi2'] < 40

    posterior = moana.SampledPosterior(mcmc_out[mask], labels=labels)
    fig, ax = posterior.corner_plot(credible_intervals=False, diagonal='chi2', 
        rotation=45, show_samples=True, labels=label_display);

    bf = mcmc_out[mcmc_out['chi2']==np.min(mcmc_out['chi2'])][labels].values[0]
    [ax[j][i].scatter(bf[i], bf[j], s=20, lw=0.8, marker='+', c='w') 
     for i in range(len(labels)) for j in range(len(labels)) if i<j]

    fig.subplots_adjust(hspace=0.15, wspace=0.15)

    # Save the figure
    fig.savefig('correlation.png', bbox_inches='tight', dpi=300, pad_inches=0.1)
