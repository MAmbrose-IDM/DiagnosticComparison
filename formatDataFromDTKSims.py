import pandas as pd
import numpy as np
import pickle

from scipy.stats import norm, binom

# import os
# os.chdir('C:\\Users\\mambrose\\OneDrive - IDMOD\\MalariaDiagnostics\\Serology')

# Read data from file 'filename.csv'
# (in the same directory that your python process is based)
# Control delimiters, rows, column names with read_csv (see later)
allSimResults = pd.read_csv('C:/Users/mambrose/Dropbox (IDM)/Malaria Team Folder/projects/serological_diagnostic/simulation_data/v1_no_HS_all_data.csv')

# maximum individuals that could test positive in a node (should be something like population size)
max_pop_size = 500
false_pos_rates = [0.001, 0.002, 0.001]

# current scenarios are various combinations of sampling date and larval habitat (proxy for intensity) - there are
#  four different values for each, giving a total of 16 scenarios to explore
all_sampling_dates = allSimResults.survey_date.unique()
all_LHs = allSimResults.x_Temporary_Larval_Habitat.unique()
all_node_ids = allSimResults.node_id.unique()
all_sim_nums = allSimResults.Run_Number.unique()
test_col_names = ['RDTTestStatus', 'hsRDTTestStatus', 'SerologyTestStatus']
test_pos_names = ['RDTPositive', 'hsRDTPositive', 'SeroPositive']

# how many individuals in each simulation/realization?
pop_size_sim = [None] * (len(all_node_ids) * len(all_sim_nums) * len(all_sampling_dates) * len(all_LHs))
index = 0

for s1 in range(len(all_sampling_dates)):
    s1_cur_rows = (allSimResults['survey_date'] == all_sampling_dates[s1])
    s1_cur_dataframe = allSimResults[s1_cur_rows]
    for s2 in range(len(all_LHs)):
        s1s2_cur_dataframe = s1_cur_dataframe[(s1_cur_dataframe['x_Temporary_Larval_Habitat'] == all_LHs[s2])]

        # count the number of times we get a certain number of positive individuals across all realization of this
        #    scenario (for now, count both nodes independently for more power)
        freq_pos_counts_circulation = [([0] * max_pop_size) for y in range(len(test_col_names))]
        freq_pos_counts_no_circulation = [([0] * max_pop_size) for y in range(len(test_col_names))]

        for nn in all_node_ids:
            s1s2nn_cur_dataframe = s1s2_cur_dataframe[(s1s2_cur_dataframe['node_id'] == nn)]
            for sim in range(len(all_sim_nums)):
                cur_rows = s1s2nn_cur_dataframe['Run_Number'] == sim
                cur_dataframe = s1s2nn_cur_dataframe[cur_rows]

                # save the number of individuals in this village at time of sampling
                pop_size_sim[index] = sum(cur_rows)

                for test in range(len(test_col_names)):
                    # count the number of positive results for each test within this simulated dataset
                    num_pos_cur = cur_dataframe.loc[:, test_col_names[test]] == test_pos_names[test]
                    freq_pos_counts_circulation[test][sum(num_pos_cur)] += 1

                    # false positives if no circulation
                    num_false_pos_cur = np.random.binomial(n=pop_size_sim[index], p=false_pos_rates[test])
                    freq_pos_counts_no_circulation[test][num_false_pos_cur] += 1
                index += 1

        num_sims = sum(freq_pos_counts_circulation[0])
        prob_num_pos_circulation = [[freq_pos_counts_circulation[x][y] / num_sims for y in range(len(freq_pos_counts_circulation[x]))]
                                    for x in range(len(freq_pos_counts_circulation))]
        prob_num_pos_no_circulation = [[freq_pos_counts_no_circulation[x][y] / num_sims for y in range(len(freq_pos_counts_no_circulation[x]))]
                                       for x in range(len(freq_pos_counts_no_circulation))]
        # save the counts for this scenario as a pickle file
        with open("simOutputs_DTK/prob_num_pos_circulation_samplingDate%i_xLH%i.p" % (all_sampling_dates[s1],
                                                                                      round(all_LHs[s2] * 100)), "wb") as f:
            pickle.dump(prob_num_pos_circulation, f)

        with open("simOutputs_DTK/prob_num_pos_no_circulation_samplingDate%i_xLH%i.p" % (all_sampling_dates[s1],
                                                                                         round(all_LHs[s2] * 100)), "wb") as f:
            pickle.dump(prob_num_pos_no_circulation, f)

# save the population sizes across all simulations
with open("simOutputs_DTK/pop_size_sim_all.p", "wb") as f:
    pickle.dump(pop_size_sim, f)




# The probability of getting a certain number of positive individuals in the population for each scenario was obtained
#    from a limited number of simulations. Although it likely approximates the true probability distribution, rarer
#    events may not have been observed. As such, we could potentially fit a probability distribution to the observed
#    results if we are willing to specify the form of this distribution. For instance, below we describe the
#    probabilities using a binomial distribution
# # load the population sizes across all simulations and take average
# with open("simOutputs_DTK/pop_size_sim_all.p", "rb") as f:
#     pop_size_sim = pickle.load(f)
# pop_size = int(round(np.mean(pop_size_sim)))
#
# for s1 in range(len(all_sampling_dates)):
#     for s2 in range(len(all_LHs)):
#         # Given circulation
#         # load file describing the probability of each possible number of positive individuals for this scenario
#         with open("simOutputs_DTK/prob_num_pos_circulation_samplingDate%i_xLH%i.p" % (all_sampling_dates[s1],
#                                                                                       round(all_LHs[s2] * 100)), "rb") as f:
#             prob_num_pos_circulation = pickle.load(f)
#
#         # fit a binomial distribution to these values
#         data_0_circulation = [[y] * int(np.round(prob_num_pos_circulation[0][y] * 1000)) for y in
#                   range(len(prob_num_pos_circulation[0]))]
#         data_circulation = [val for sublist in data_0_circulation for val in sublist]
#
#         # Fit a binomial distribution to the data:
#         n_b, p_b = pop_size, np.mean(data_circulation) / pop_size
#
#         # Given no circulation
#         # load file describing the probability of each possible number of positive individuals for this scenario
#         with open("simOutputs_DTK/prob_num_pos_no_circulation_samplingDate%i_xLH%i.p" % (all_sampling_dates[s1],
#                                                                                          round(all_LHs[s2] * 100)), "rb") as f:
#             prob_num_pos_no_circulation = pickle.load(f)
#
#         # fit a normal distribution to these values
#         data_0_no_circulation = [[y] * int(np.round(prob_num_pos_circulation[0][y] * 1000)) for y in
#                   range(len(prob_num_pos_circulation[0]))]
#         data_no_circulation = [val for sublist in data_0_no_circulation for val in sublist]
#
#         # Fit a binomial distribution to the data:
#         n_b, p_b = pop_size, np.mean(data_no_circulation) / pop_size
#
# # could then calculate and save probabilities for different values of s_n
#

# # some plots to explore how well this would work...
# import numpy as np
# from scipy.stats import norm
# import matplotlib.pyplot as plt
#
#
# data_0 = [[y]*int(np.round(prob_num_pos_no_circulation[0][y]*1000)) for y in range(len(prob_num_pos_no_circulation[0]))]
# data = [val for sublist in data_0 for val in sublist]
#
# # Fit a normal distribution to the data:
# mu, std = norm.fit(data)
# n_b, p_b = 271, np.mean(data)/271
#
# # Plot the histogram.
# plt.hist(data, bins=3, density=True, alpha=0.6, color='g')
#
# # Plot the PDF.
# xmin, xmax = plt.xlim()
# x = np.linspace(xmin, xmax, 100)
# p_norm = norm.pdf(x, mu, std)
# p_binom = binom.pmf(x, n_b, p_b)
# plt.plot(x, p_norm, 'k', linewidth=2)
# # plt.plot(x, p_binom, 'k', linewidth=2, color='red')
# x_binom = range(int(np.round(xmin)), int(np.round(xmax)))
# plt.plot(x_binom, binom.pmf(x_binom, n_b, p_b), 'bo', ms=8, label='binom pmf')
# title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
# plt.title(title)
#
# plt.show()







