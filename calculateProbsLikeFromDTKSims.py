# Create pickle files for:
# 1) the probabilities of detecting one or more positives in S samples given circulation (for each scenario)
# 2) the probabilities of detecting one or more positives in S samples given NO circulation (for each scenario)
# 3) the likelihood of no circulation given that S(+) out of S samples test positive (for each scenario, different
#     values of S, and different values of P(circulation)

import numpy as np
import pickle
from operator import mul
from functools import reduce

# To greatly reduce the difficulty of this problem, we assume the population size is constant at the mean value across
#     all DTK simulations. This assumption may possibly give rise to some biases and will ideally be explored through
#     a sensitivity analysis in the future (something as easy as using the max and min values and seeing whether
#     results change substantially).
# load the population sizes across all simulations and take average
with open("simOutputs_DTK/pop_size_sim_all.p", "rb") as f:
    pop_size_sim = pickle.load(f)
pop_size = int(round(np.mean(pop_size_sim)))

# Scenario names
all_sampling_dates = [45, 190, 300, 360]
all_LHs = [0.6, 0.7, 0.8, 1.0]

# maximum individuals that could test positive in a node (should be something like population size)
max_pop_size = 500
false_pos_rates = [0.001, 0.002, 0.001]

# Test names
test_names = ['RDT', 'hsRDT', 'Serology']

# Probabilities of detecting one or more positives in S samples given circulation and given no circulation
#   (repeat for each scenario)
for s1 in range(len(all_sampling_dates)):
    for s2 in range(len(all_LHs)):
        # load files describing the probability of each possible number of positive individuals for this scenario
        with open("simOutputs_DTK/prob_num_pos_circulation_samplingDate%i_xLH%i.p" % (all_sampling_dates[s1],
                                                                                      round(all_LHs[s2] * 100)), "rb") as f:
            freq_pos_counts_circulation = pickle.load(f)
        with open("simOutputs_DTK/prob_num_pos_no_circulation_samplingDate%i_xLH%i.p" % (all_sampling_dates[s1],
                                                                                         round(all_LHs[s2] * 100)), "rb") as f:
            freq_pos_counts_no_circulation = pickle.load(f)

        # create lists to store the probabilities for all values of S and all tests
        prob_pos_sample_given_circulation = [([0] * (pop_size+1)) for y in range(len(test_names))]
        prob_pos_sample_given_no_circulation = [([0] * (pop_size+1)) for y in range(len(test_names))]

        # iterate over values of S (from 0 to pop_size)
        for ss in range(pop_size+1):

            # Let p1 = P(S_p >= 1 | N_p=n_p ) = 1-P(S_p == 0 | N_p=n_p )
            # Let p2 = P(N_p == n_p | circulation)
            # Then P(S_p >= 1 | circulation) = sum from n_p=0 to n_p=pop_size of p1*p2
            # S_p = the total number of individuals among those sampled who test positive
            # N_p = the total number of individuals in the population who would have been positive if tested

            # iterate over values of n_p, calculating p1, p2, and adding the product to prob_at_least_one_positive
            prob_at_least_one_positive_circulation = [0]*len(test_names)
            prob_at_least_one_positive_no_circulation = [0]*len(test_names)
            for n_p in range(pop_size+1):
                # probability of S_p==0 observed positive samples given this value of n_p
                if n_p == 0:
                    p1 = 0  # no chance of getting a positive sample no matter what the value of ss
                elif ss == 0:
                    p1 = 0  # no chance of a positive sample if you don't collect any samples
                elif (pop_size - n_p) >= ss:
                    # calculate ((N-n_p) choose S) / (N choose S), where N=pop_size, S=number of samples collected (ss)
                    # (N-n_p)_choose_S = (N-n_p)!/S!/(N-n_p-S)! = (N-n_p)*...*(N-n_p-S+1)/S!
                    # N_choose_S = N!/S!/(N-S)! = N*...*(N-S+1)/S!
                    # (N - n_p)_choose_S / N_choose_S = (N-n_p)*...*(N-n_p-S+1)/(N*...*(N-S+1))
                    # note that range(x) goes up to (x-1), so we add one to the upper values
                    p1 = 1 - (reduce(mul, list(range((pop_size - n_p - ss + 1),
                                                     (pop_size - n_p + 1))))
                              / reduce(mul, list(range((pop_size - ss + 1), (pop_size + 1)))))

                else:
                    p1 = 1 - 0
                # iterate through tests, each of which will have a different value for p2
                for test in range(len(test_names)):
                    p2 = freq_pos_counts_circulation[test][n_p]
                    prob_at_least_one_positive_circulation[test] += p1*p2
                    p2_no_circulation = freq_pos_counts_no_circulation[test][n_p]
                    prob_at_least_one_positive_no_circulation[test] += p1*p2_no_circulation

            # store the results for this value of ss
            for test in range(len(test_names)):
                prob_pos_sample_given_circulation[test][ss] = prob_at_least_one_positive_circulation[test]
                prob_pos_sample_given_no_circulation[test][ss] = prob_at_least_one_positive_no_circulation[test]

        # save the nested list containing the relevant probabilities for this scenario
        with open("simOutputs_DTK/prob_pos_sample_given_circulation_samplingDate%i_xLH%i.p" % (all_sampling_dates[s1],
                                                                                               round(all_LHs[s2] * 100)), "wb") as f:
            pickle.dump(prob_pos_sample_given_circulation, f)

        with open("simOutputs_DTK/prob_pos_sample_given_no_circulation_samplingDate%i_xLH%i.p" % (all_sampling_dates[s1],
                                                                                               round(all_LHs[s2] * 100)), "wb") as f:
            pickle.dump(prob_pos_sample_given_no_circulation, f)





