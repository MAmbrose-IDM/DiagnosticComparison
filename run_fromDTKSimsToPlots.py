from formatDataFromDTKSims import formatDataFromDTKSims

from calculateProbsLikeFromDTKSims import calculate_prob_positive_sample, calculate_likelihood_given_observation

import os
os.chdir('C:\\Users\\mambrose\\OneDrive - IDMOD\\MalariaDiagnostics\\Serology')

#
# Parameters
#
simResults_filename = 'C:/Users/mambrose/Dropbox (IDM)/Malaria Team Folder/projects/serological_diagnostic/simulation_data/test_all_data.csv'
filename_suffix = 'test1'

# maximum individuals that could test positive in a node (should be something like population size)
max_pop_size = 500
false_pos_rates = [0.001, 0.002, 0.001]
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








