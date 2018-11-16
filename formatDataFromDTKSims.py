import pandas as pd
import numpy as np
import pickle

from scipy.stats import norm, binom

def formatDataFromDTKSims(simResults_filename, false_pos_rates, test_col_names, test_pos_names, max_pop_size,
                          filename_suffix):
    """
    Read in output from many DTK simulations, all stored in a csv file. The input file must have columns for the
        test results, node (village), day of sampling, LH multiplier, and realization number.
        Save a .p pickle file for each sampling date and LH multiplier scenario. The pickle file contains a nested list,
        where the outer list gives the test and ith element of the inner list gives the probability that i individuals
        would have a positive test if surveyed.
    :param simResults_filename: path and name of file where DTK output was stored
    :param false_pos_rates: for each of the tests, give the probability a test will be positive even if the
        individual has had no exposure
    :param test_col_names: names of the columns for each of the diagnostic tests in the DTK simResults file
    :param test_pos_names: values used to indicate a positive results for each of the tests in the DTK simResults file
    :param max_pop_size: maximum population size we will encounter across all simulations in this input file
    :param filename_suffix: string to add to end of newly created files to identify simulation batch
    :return:
    """
    allSimResults = pd.read_csv(simResults_filename)

    # current scenarios are various combinations of sampling date and larval habitat (proxy for intensity) - there are
    #  four different values for each, giving a total of 16 scenarios to explore
    all_sampling_dates = allSimResults.survey_date.unique()
    all_LHs = allSimResults.x_Temporary_Larval_Habitat.unique()
    all_node_ids = allSimResults.node_id.unique()
    all_sim_nums = allSimResults.Run_Number.unique()

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

                        # add in false positives among those with a negative test
                        num_neg_cur = pop_size_sim[index] - sum(num_pos_cur)
                        num_false_pos = np.random.binomial(n=num_neg_cur, p=false_pos_rates[test])

                        # total number of positive individuals = true positives plus false positives
                        freq_pos_counts_circulation[test][np.sum(num_pos_cur) + num_false_pos] += 1

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
            with open("simOutputs_DTK/prob_num_pos_circulation_samplingDate%i_xLH%i_%s.p" % (all_sampling_dates[s1],
                                                                                             round(all_LHs[s2] * 100),
                                                                                             filename_suffix), "wb") as f:
                pickle.dump(prob_num_pos_circulation, f)

            with open("simOutputs_DTK/prob_num_pos_no_circulation_samplingDate%i_xLH%i_%s.p" % (all_sampling_dates[s1],
                                                                                                round(all_LHs[s2] * 100),
                                                                                                filename_suffix), "wb") as f:
                pickle.dump(prob_num_pos_no_circulation, f)

    # save the population sizes across all simulations
    with open("simOutputs_DTK/pop_size_sim_all_%s.p" % filename_suffix, "wb") as f:
        pickle.dump(pop_size_sim, f)

    # return sampling dates and LH values
    return [all_sampling_dates, all_LHs]


def formatPrevDataFromDTKSims(prev_filename, filename_suffix):
    """
    Read in output from many DTK simulations, all stored in a csv file. The input file must have a column for the
        malaria prevalence at each day.
        Save a .p pickle file for each larval habitat scenario containing three lists.
        The first list contains the mean prevalence for each day, the second list contains the 95%CI max
        (97.5th percentile), and the third list contains the 95%CI min (the 2.5th percentile).
    :param prev_filename: path and name of file where DTK output was stored
    :param filename_suffix: string to add to end of newly created files to identify simulation batch
    """
    allSimResults = pd.read_csv(prev_filename)

    # four different larval habitat multipliers (proxy for intensity) - there are
    #  four different values for each, giving a total of 16 scenarios to explore
    all_LHs = allSimResults.x_Temporary_Larval_Habitat.unique()
    all_node_ids = allSimResults.node.unique()
    all_sim_nums = allSimResults.Run_Number.unique()


    for lh in range(len(all_LHs)):
        # store a list of the prevalences recorded on each day for this scenario
        #   outer list is day, inner list is the prevalences recorded on that day across all simulations
        # prev_list = [([None] * (len(all_node_ids) * len(all_sim_nums))) for y in np.max(allSimResults['time'])]
        prev_list = [[None] * (np.max(allSimResults['time']) + 1) for y in range(3)]

        for dd in range(np.max(allSimResults['time'])+1):
            lh_rows = allSimResults['x_Temporary_Larval_Habitat'] == all_LHs[lh]
            dd_rows = allSimResults['time'] == dd
            new_df = allSimResults[lh_rows & dd_rows]
            prev_values = list(new_df.Prevalence)
            prev_list[0][dd] = np.median(prev_values)
            prev_list[1][dd] = np.quantile(prev_values, q=0.975)
            prev_list[2][dd] = np.quantile(prev_values, q=0.025)

            # save the mean and 95CI values for this scenario as a pickle file
            with open("simOutputs_DTK/prevalenceData_xLH%i_%s.p" % (round(all_LHs[lh] * 100),
                                                                    filename_suffix), "wb") as f:
                pickle.dump(prev_list, f)




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







