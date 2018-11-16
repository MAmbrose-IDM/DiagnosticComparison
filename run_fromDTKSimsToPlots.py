import pickle
import os
import numpy as np

from formatDataFromDTKSims import formatDataFromDTKSims, formatPrevDataFromDTKSims
from calculateProbsLikeFromDTKSims import calculate_prob_positive_sample, calculate_likelihood_given_observation
from createPlotPanelsFromDTKSims import plot_panel_prob_num_pos_circulation, plot_panel_prob_num_pos_no_circulation,\
    plot_panel_prob_positive_sample_circulation, plot_panel_prob_positive_sample_no_circulation, \
    plot_panel_likelihood_circulation, plot_panel_likelihood_no_circulation, \
    plot_panel_prevalence, plot_panel_prevalence_with_without_HS

#
# Parameters
#
os.chdir('C:\\Users\\mambrose\\OneDrive - IDMOD\\MalariaDiagnostics\\Serology')

# simResults_filename = 'C:/Users/mambrose/Dropbox (IDM)/Malaria Team Folder/projects/serological_diagnostic/simulation_data/v3_noHS_all_data.csv'
# prev_filename = 'C:/Users/mambrose/Dropbox (IDM)/Malaria Team Folder/projects/serological_diagnostic/simulation_data/v3_noHS_prevalence_data.csv'
# filename_suffix = 'v3'
simResults_filename = 'C:/Users/mambrose/Dropbox (IDM)/Malaria Team Folder/projects/serological_diagnostic/simulation_data/v3_withHS_all_data.csv'
prev_filename = 'C:/Users/mambrose/Dropbox (IDM)/Malaria Team Folder/projects/serological_diagnostic/simulation_data/v3_withHS_prevalence_data.csv'
filename_suffix = 'v3_withHS'


# plot with and without HS?
plot_with_without_hs = True
filename_suffix1 = 'v3'
filename_suffix2 = 'v3_withHS'

need_to_generate_pickle_files = True
need_to_generate_plots = True

# maximum individuals that could test positive in a node (should be something like population size)
max_pop_size = 500
false_pos_rates = [0.05, 0.08, 0.03]
# test_col_names = ['RDTTestStatus', 'hsRDTTestStatus', 'SerologyTestStatus']
# test_pos_names = ['RDTPositive', 'hsRDTPositive', 'SeroPositive']
test_col_names = ['RDTPositive', 'hsRDTPositive', 'SerologyTestStatus']
test_pos_names = [1, 1, 'SeroPositive']
test_names = ['RDT', 'hsRDT', 'Serology']
# Additional parameters needed to look at the likelihood of circulation (or no circulation) under different assumptions
p_circulation = [0.1, 0.3, 0.5]
ss_values = [int(round(300 * y)) for y in [0.1, 0.3, 0.5]]


#
# save reformatted data files as pickle files
#
if need_to_generate_pickle_files == True:

    print('creating prevalence files...')
    formatPrevDataFromDTKSims(prev_filename, filename_suffix)

    print('creating probability density files for number of individuals positive in the population...')

    output = formatDataFromDTKSims(simResults_filename=simResults_filename, false_pos_rates=false_pos_rates,
                                   test_col_names=test_col_names, test_pos_names=test_pos_names,
                                   max_pop_size=max_pop_size, filename_suffix=filename_suffix)

    all_sampling_dates = output[0]
    all_LHs = output[1]
    #
    # calculate the probabilities of at least one positive sample given circulation or given no circulation
    #
    print('calculating probabilities that at least one sample will be positive...')
    calculate_prob_positive_sample(filename_suffix=filename_suffix, all_sampling_dates=all_sampling_dates,
                                   all_LHs=all_LHs, test_names=test_names)


    #
    # calculate the likelihood of circulation or of no circulation given observed sample results
    #
    print('calculating likelihood of circulation or no circulation given observations...')
    calculate_likelihood_given_observation(filename_suffix=filename_suffix, all_sampling_dates=all_sampling_dates,
                                           all_LHs=all_LHs, test_names=test_names,
                                           p_circulation=p_circulation, ss_values=ss_values)

else:
    all_sampling_dates = [45, 190, 300, 360]
    all_LHs = [0.6, 0.7, 0.8, 1.0]

sampling_date_names = ['Feb. 14', 'Jul. 9', 'Oct. 27', 'Dec. 26']
# calculate the average prevalence (equally weighted across all days of the year) for each LH scenario
print('calculating average prevalences...')
ave_prevs = [None] * len(all_LHs)
for lh in range(len(all_LHs)):
    with open("simOutputs_DTK/prevalenceData_xLH%i_%s.p" % (round(all_LHs[lh] * 100),
                                                            filename_suffix), "rb") as f:
        prev_list = pickle.load(f)
    ave_prevs[lh] = np.mean(prev_list[0])

# Plots
if need_to_generate_plots == True:
    print('creating plot panels...')

    # To greatly reduce the difficulty of this problem, we assume the population size is constant at the mean value across
    #     all DTK simulations. This assumption may possibly give rise to some biases and will ideally be explored through
    #     a sensitivity analysis in the future (something as easy as using the max and min values and seeing whether
    #     results change substantially).
    # load the population sizes across all simulations and take average
    with open("simOutputs_DTK/pop_size_sim_all_%s.p" % filename_suffix, "rb") as f:
        pop_size_sim = pickle.load(f)
    pop_size = int(round(np.mean(pop_size_sim)))

    #
    # Prevalence of infection through the year
    #
    plot_panel_prevalence(filename_suffix=filename_suffix, all_sampling_dates=all_sampling_dates, all_LHs=all_LHs,
                                        sampling_date_names=sampling_date_names, ave_prevs=ave_prevs)

    if plot_with_without_hs:
        plot_panel_prevalence_with_without_HS(filename_suffix1=filename_suffix1, filename_suffix2=filename_suffix2,
                                              all_sampling_dates=all_sampling_dates, all_LHs=all_LHs,
                                              sampling_date_names=sampling_date_names, ave_prevs=ave_prevs)


    #
    # Probability distribution for number of individuals positive given circulation
    #
    plot_panel_prob_num_pos_circulation(pop_size=pop_size, test_names=test_names, filename_suffix=filename_suffix,
                                        all_sampling_dates=all_sampling_dates, all_LHs=all_LHs,
                                        sampling_date_names=sampling_date_names, ave_prevs=ave_prevs)


    #
    # Probability distribution for number of individuals positive given no circulation
    #
    plot_panel_prob_num_pos_no_circulation(pop_size=pop_size, test_names=test_names, filename_suffix=filename_suffix,
                                           all_sampling_dates=all_sampling_dates, all_LHs=all_LHs)


    #
    # Probability of at least one positive sample given circulation
    #
    plot_panel_prob_positive_sample_circulation(pop_size=pop_size, test_names=test_names,
                                                filename_suffix=filename_suffix, all_sampling_dates=all_sampling_dates,
                                                all_LHs=all_LHs, sampling_date_names=sampling_date_names,
                                                ave_prevs=ave_prevs)


    #
    # Probability of at least one positive sample given no circulation
    #
    plot_panel_prob_positive_sample_no_circulation(pop_size=pop_size, test_names=test_names,
                                                   filename_suffix=filename_suffix, all_sampling_dates=all_sampling_dates,
                                                   all_LHs=all_LHs)


    #
    # Likelihood of circulation given observation
    #
    plot_panel_likelihood_circulation(pop_size=pop_size, test_names=test_names, filename_suffix=filename_suffix,
                                      all_sampling_dates=all_sampling_dates, all_LHs=all_LHs, ss_values=ss_values,
                                      p_circulation=p_circulation, sampling_date_names=sampling_date_names,
                                      ave_prevs=ave_prevs)


    #
    # Likelihood of no circulation given observation
    #
    plot_panel_likelihood_no_circulation(pop_size=pop_size, test_names=test_names, filename_suffix=filename_suffix,
                                      all_sampling_dates=all_sampling_dates, all_LHs=all_LHs, ss_values=ss_values,
                                      p_circulation=p_circulation, sampling_date_names=sampling_date_names,
                                      ave_prevs=ave_prevs)












