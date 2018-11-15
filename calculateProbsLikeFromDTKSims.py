# Create pickle files for:
# 1) the probabilities of detecting one or more positives in S samples given circulation (for each scenario)
# 2) the probabilities of detecting one or more positives in S samples given NO circulation (for each scenario)
# 3) the likelihood of no circulation given that S(+) out of S samples test positive (for each scenario, different
#     values of S, and different values of P(circulation)

import numpy as np
import pickle
from operator import mul
from functools import reduce

# Quick key
# S_p = the total number of individuals among those sampled who test positive
# S_n = the total number of individuals among those sampled who test negative
# N_p = the total number of individuals in the population who would have been positive if tested
# N_n = the total number of individuals in the population who would have been negative if tested

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
            prob_num_pos_circulation = pickle.load(f)
        with open("simOutputs_DTK/prob_num_pos_no_circulation_samplingDate%i_xLH%i.p" % (all_sampling_dates[s1],
                                                                                         round(all_LHs[s2] * 100)), "rb") as f:
            prob_num_pos_no_circulation = pickle.load(f)

        # create lists to store the probabilities for all values of S and all tests
        prob_pos_sample_given_circulation = [([0] * (pop_size+1)) for y in range(len(test_names))]
        prob_pos_sample_given_no_circulation = [([0] * (pop_size+1)) for y in range(len(test_names))]

        # iterate over values of S (from 0 to pop_size)
        for ss in range(pop_size+1):

            # Let p1 = P(S_p >= 1 | N_p=n_p ) = 1-P(S_p == 0 | N_p=n_p )
            # Let p2 = P(N_p == n_p | circulation)
            # Then P(S_p >= 1 | circulation) = sum from n_p=0 to n_p=pop_size of p1*p2

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
                    p2 = prob_num_pos_circulation[test][n_p]
                    prob_at_least_one_positive_circulation[test] += p1*p2
                    p2_no_circulation = prob_num_pos_no_circulation[test][n_p]
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


# Calculate likelihood of circulation given s_p positive samples out of ss total samples
p_circulation = [0.1, 0.3, 0.5]
ss_values = [int(round(pop_size * y)) for y in [0.1, 0.3, 0.5]]


# Probabilities of detecting one or more positives in S samples given circulation and given no circulation
#   (repeat for each scenario)
for s1 in range(len(all_sampling_dates)):
    for s2 in range(len(all_LHs)):
        # load files describing the probability of each possible number of positive individuals for this scenario
        with open("simOutputs_DTK/prob_num_pos_circulation_samplingDate%i_xLH%i.p" % (all_sampling_dates[s1],
                                                                                      round(all_LHs[s2] * 100)), "rb") as f:
            prob_num_pos_circulation = pickle.load(f)
        with open("simOutputs_DTK/prob_num_pos_no_circulation_samplingDate%i_xLH%i.p" % (all_sampling_dates[s1],
                                                                                         round(all_LHs[s2] * 100)), "rb") as f:
            prob_num_pos_no_circulation = pickle.load(f)

        # iterate over values of S (from the ss_values list)
        for ss_index, ss in enumerate(ss_values):

            # create lists to store the likelihoods for all values of s_n and all tests. For calculation efficiency,
            #    we have this nested list nested within an outer list for each of the p_circulation values, though the
            #    output for different p_circulation values will be saved in separate files
            #    Indexing will be as follows: [prob_circulation_index][test_index][s_n]
            like_circulation_given_s_n = [[([0] * (ss + 1)) for y in range(len(test_names))] for z in range(len(p_circulation))]
            like_no_circulation_given_s_n = [[([0] * (ss + 1)) for y in range(len(test_names))] for z in range(len(p_circulation))]
            # store sum (as we iterate over values of n_p) for each component
            circulation_sum_component = [[([0] * (ss + 1)) for y in range(len(test_names))] for z in range(len(p_circulation))]
            no_circulation_sum_component = [[([0] * (ss + 1)) for y in range(len(test_names))] for z in range(len(p_circulation))]

            # iterate over possible values of s_n (must be less than or equal to ss
            for s_n in range(ss + 1):
                # Calculate L(no circulation | S_n==s_n) ...and repeat for each value of p_circulation and for each test:
                #  = P(S_n==s_n | no circulation) * P(no circulation) / (P(S_n==s_n | no circulation) * P(no circulation)
                #                                                        + P(S_n==s_n | circulation) * P(circulation))

                # Let p1 = P(S_n == s_n | N_p=n_p )
                # Let p2_circulation = P(N_p == n_p | circulation)
                # Let p2_no_circulation = P(N_p == n_p | no circulation)
                # Then, for a given set of values for ss, s_n, prob_circulation,
                # L(no circulation | S_n==s_n) = sum({n_p=0->pop_size} of (p1 * p2_no_circulation * prob_no_circulation))
                #     / (  sum({n_p=0->pop_size} of (p1 * p2_no_circulation * prob_no_circulation))
                #        + sum({n_p=0->pop_size} of (p1 * p2_circulation * prob_circulation))    )

                # iterate over values of n_p, calculating p1, p2, p2_no_circulation, and adding the product to
                #     the relevant list positions
                for n_p in range(pop_size+1):
                    # probability of S_n==s_n observed positive samples given this value of n_p
                    if (pop_size-n_p) < s_n:
                        p1 = 0  # no chance of sampling more negatives than there are in the entire population
                    elif n_p < (ss-s_n):
                        p1 = 0  # no chance of sampling more positives than there are in the entire population
                    else:
                        # calculate ((N-n_p) choose (s_n)) * ((n_p) choose (ss-s_n)) / (N choose ss), where N=pop_size
                        #   = ( (N-n_p-s_n+1)*...*(N-n_p) ) * ( (n_p-(ss-s_n)+1)*...*(n_p) ) * ( (ss-s_n+1)*...*(ss) ) /
                        #       ( (1)*...*(s_n) ) / ( (N-ss+1)*...*(N) )
                        # note that range(x) goes up to (x-1), so we add one to the upper values
                        if s_n == 0:
                            # p1 = (reduce(mul, list(range((n_p - ss + 1), (n_p + 1))))
                            #           / reduce(mul, list(range((pop_size - ss + 1), (pop_size + 1)))))
                            p1_c_log = ( (sum([np.log(y) for y in range((n_p - (ss-s_n) + 1), (n_p + 1))]))
                                         - (sum([np.log(y) for y in range((pop_size - ss + 1), (pop_size + 1))])) )
                            p1 = np.exp(p1_c_log)
                        elif s_n == ss:
                            # p1 = (reduce(mul, list(range((pop_size - n_p - ss + 1), (pop_size - n_p + 1))))
                            #           / reduce(mul, list(range((pop_size - ss + 1), (pop_size + 1)))))
                            p1_c_log = ( (sum([np.log(y) for y in range((pop_size - n_p - s_n + 1), (pop_size - n_p + 1))]))
                                         - (sum([np.log(y) for y in range((pop_size - ss + 1), (pop_size + 1))])) )
                            p1 = np.exp(p1_c_log)
                        else:
                            # p1 = (reduce(mul, list(range((pop_size - n_p - s_n + 1), (pop_size - n_p + 1))))
                            #           * reduce(mul, list(range((n_p - (ss-s_n) + 1), (n_p + 1))))
                            #           * reduce(mul, list(range(((ss-s_n) + 1), (ss + 1))))
                            #           / reduce(mul, list(range(1, (s_n + 1))))
                            #           / reduce(mul, list(range((pop_size - ss + 1), (pop_size + 1)))))

                            # re-do, now using logs to avoid overflow errors
                            # the log of the complement of p1
                            p1_c_log = ( (sum([np.log(y) for y in range((pop_size - n_p - s_n + 1), (pop_size - n_p + 1))]))
                                         + (sum([np.log(y) for y in range((n_p - (ss-s_n) + 1), (n_p + 1))]))
                                         + (sum([np.log(y) for y in range(((ss-s_n) + 1), (ss + 1))]))
                                         - (sum([np.log(y) for y in range(1, (s_n + 1))]))
                                         - (sum([np.log(y) for y in range((pop_size - ss + 1), (pop_size + 1))])) )
                            p1 = np.exp(p1_c_log)

                    for p_c in range(len(p_circulation)):
                        for test in range(len(test_names)):
                            p2_circulation = prob_num_pos_circulation[test][n_p]
                            p2_no_circulation = prob_num_pos_no_circulation[test][n_p]
                            prob_circulation = p_circulation[p_c]
                            circulation_sum_component[p_c][test][s_n] += (p1 * p2_circulation * prob_circulation)
                            no_circulation_sum_component[p_c][test][s_n] += (p1 * p2_no_circulation * (1-prob_circulation))

            # after accumulating all the sums, calculate the final likelihood for each value of s_n, test, prob_circulation
            for p_c in range(len(p_circulation)):
                for test in range(len(test_names)):
                    for s_n in range(ss+1):
                        denom_value = ((circulation_sum_component[p_c][test][s_n] + no_circulation_sum_component[p_c][test][s_n]))
                        if denom_value > 0:
                            like_circulation_given_s_n[p_c][test][s_n] = (circulation_sum_component[p_c][test][s_n]
                                                                          / denom_value)
                            like_no_circulation_given_s_n[p_c][test][s_n] = 1 - like_circulation_given_s_n[p_c][test][s_n]
                        else:
                            like_circulation_given_s_n[p_c][test][s_n] = -9
                            like_no_circulation_given_s_n[p_c][test][s_n] = -9

                            # Save the results in pickle files
            for p_c in range(len(p_circulation)):
                # save the nested list containing the relevant probabilities for this scenario
                with open("simOutputs_DTK/lik_circulation_numSamples%i_PCirculation%i_samplingDate%i_xLH%i.p" % (ss, round(p_circulation[p_c] * 100), all_sampling_dates[s1],
                                                                                                  round(all_LHs[s2] * 100)), "wb") as f:
                    pickle.dump(like_circulation_given_s_n[p_c], f)

                with open("simOutputs_DTK/lik_no_circulation_numSamples%i_PCirculation%i_samplingDate%i_xLH%i.p" % (ss, round(p_circulation[p_c] * 100), all_sampling_dates[s1],
                                                                                                     round(all_LHs[s2] * 100)), "wb") as f:
                    pickle.dump(like_no_circulation_given_s_n[p_c], f)




