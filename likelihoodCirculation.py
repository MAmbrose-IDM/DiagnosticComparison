import numpy as np
from operator import mul
from functools import reduce


def likelihood_no_circulation_all_negative(pop_size, number_tested, number_pos_circulation, number_pos_no_circulation,
                              prob_circulation):
    """
    Calculate the likelihood that there is no circulation given that number_tested individuals are sampled
        and are all negative if we're dealing with a population where number_pos_circulation individuals would test positive
        out of a total population of pop_size individuals if there were circulation and number_pos_no_circulation
        would test positive in a population without circulation. Assume sampling without replacement.

    Arguments:
        :param pop_size - number of individuals in the sampled population
        :param number_tested - number of individuals tested
        :param number_pos_circulation - number of individuals who would have a positive diagnostic test if circulation were
            occurring
        :param number_pos_no_circulation - number of individuals who would have a positive diagnostic test if circulation
            were not occurring
        :param prob_circulation - probability circulation is occurring in any given sampled village
    """
    prob_all_negative_given_circulation = ((reduce(mul,
                                                   list(range((pop_size - number_pos_circulation - number_tested + 1),
                                                              (pop_size - number_pos_circulation + 1))))
                                            / reduce(mul, list(range((pop_size - number_tested + 1), (pop_size + 1))))))

    prob_all_negative_given_no_circulation = ((reduce(mul,
                                                      list(range((pop_size - number_pos_no_circulation - number_tested + 1),
                                                                 (pop_size - number_pos_no_circulation + 1))))
                                               / reduce(mul, list(range((pop_size - number_tested + 1), (pop_size + 1))))))

    likelihood_no_circulation = (prob_all_negative_given_no_circulation * (1 - prob_circulation)) \
        / ((prob_all_negative_given_no_circulation * (1 - prob_circulation))
           + (prob_all_negative_given_circulation * prob_circulation))

    return likelihood_no_circulation


def calculate_likelihood_no_circulation_each_sim(sim_output_list, pop_size, number_tested, false_pos_prob,
                                              prob_circulation):
    """

    :param sim_output_list: number of individuals whose test would be positive if they were sampled on the
        surveillance day if circulation were occurring.
        It is a list of list of lists. This is a nested list five levels deep. The outermost level
        is the surveillance scenario, the next should be set to one to give number of tests positive, the next is the
        diagnostic test, the next is the transmission intensity and the inner-most list has all results from different
        simulation runs.
    :param pop_size:  number of individuals in surveilled population
    :param number_tested - number of individuals tested
    :param false_pos_prob: for each diagnostic test, what is the probability an individual without exposure tests
        positive
    :param prob_circulation - probability circulation is occurring in any given sampled village
    :return: the number of samples that would need to be collected for each simulated population so that the probability
        of detecting a positive test would be at least confidence_level. Returned in the same format as
        sim_number_positive (with the same list levels).
    """
    # save results in likelihood_circulation_all_scenarios_sims
    likelihood_no_circulation_all_scenarios_sims = [[[[[None]*len(sim_output_list[0][1][0][0])
                                                    for m in range(len(sim_output_list[0][1][0]))]
                                                    for n in range(len(sim_output_list[0][1]))]
                                                    for o in range(len(sim_output_list[0]))]
                                                    for p in range(len(sim_output_list))]

    # surveillance scenario
    for i1 in range(len(sim_output_list)):
        # diagnostic test
        for i3 in range(len(sim_output_list[i1][1][i3])):
            # transmission intensity
            for i4 in range(len(sim_output_list[i1][1][i3][i4])):
                # simulation run
                for i5 in range(len(sim_output_list[i1][1][i3][i4][i5])):
                    likelihood_no_circulation_all_scenarios_sims[i1][1][i3][i4][i5] = \
                        likelihood_no_circulation_all_negative(pop_size=pop_size,
                                                               number_tested=number_tested,
                                                               number_pos_circulation=sim_output_list[i1][1][i3][i4][i5],
                                                               number_pos_no_circulation=np.random.binomial(n=pop_size,
                                                                                                            p=false_pos_prob[i3]),
                                                               prob_circulation=prob_circulation)

    return likelihood_no_circulation_all_scenarios_sims




# TODO: all the rest of this script... not finished as of 11/12/2018
###############
# calculate how many samples are needed to have a likelihood of no circulation less than a threshold if all samples are negative


def calculate_num_samples_needed_no_circulation(pop_size, number_pos_circulation, number_pos_no_circulation,
                                                prob_circulation, confidence_level):
    """

    Given information about the number of individuals who would test positive on a given surveillance date if
        circulation were occurring, the number that would test positive if circulation were not occurring, the
        total population size, and the probability a village has circulation,  return the number of individuals
        that would need to be sampled to result in a certain likelihood that a village with all negative samples
        is truly negative.  Assume sampling without replacement (so we're calculating the number of
        unique individuals that should be sampled).

    :param pop_size: number of individuals in surveilled population
    :param number_pos_circulation: number of individuals that would test positive for this diagnostic test if sampled
        when circulation was occurring
    :param number_pos_no_circulation: number of individuals that would test positive for this diagnostic test if sampled
        when circulation was not occurring
    :param prob_circulation: probability a randomly selected village will be positive
    :param confidence_level: if circulation is occuring, how certain do we want to be that our sampling design will
        detect it?
    :return: the number of samples that need to be collected from this population if we want to be confidence_level
        sure that obtaining all negative samples will mean no circulation
    """

    # If all individuals would test positive if circulation were occurring, can sample one individual and be certain
    if number_pos_circulation == pop_size:
        return 1
    # # Check whether no samples are positive - if so, no number of samples will be sufficient
    # if num_test_positive == 0:
    #     return math.inf

    # Find the minimum number of samples that would need to test negative so that likelihood of no circulation probability of detection is above confidence_level
    upper_value = pop_size

    # probability infection would be detected if upper.value individuals were sampled
    if upper_value == 0:
        current_prob = 0
    else:
        current_prob = prob_detected_without_replacement(pop_size=pop_size,
                                                         number_pos=num_test_positive,
                                                         number_tested=upper_value,
                                                         times_sampled=1
                                                         )

    # currently coded to search each value; if population sizes are large, a more efficient search algorithm would
    #   be preferable
    while current_prob > confidence_level:
        # decrease the number of samples taken and recalculate the new current_prob
        upper_value -= 1
        if upper_value == 0:
            current_prob = 0
        else:
            current_prob = prob_detected_without_replacement(pop_size=pop_size,
                                                             number_pos=num_test_positive,
                                                             number_tested=upper_value,
                                                             times_sampled=1
                                                             )


def calculate_num_samples_no_circulation_each_sim(sim_output_list, pop_size, number_tested, false_pos_prob,
                                              prob_circulation, confidence_level):
    """

    :param sim_output_list: number of individuals whose test would be positive if they were sampled on the
        surveillance day if circulation were occurring.
        It is a list of list of lists. This is a nested list five levels deep. The outermost level
        is the surveillance scenario, the next should be set to one to give number of tests positive, the next is the
        diagnostic test, the next is the transmission intensity and the inner-most list has all results from different
        simulation runs.
    :param pop_size:  number of individuals in surveilled population
    :param number_tested - number of individuals tested
    :param false_pos_prob: for each diagnostic test, what is the probability an individual without exposure tests
        positive
    :param prob_circulation - probability circulation is occurring in any given sampled village
    :return: the number of samples that would need to be collected for each simulated population so that the probability
        of detecting a positive test would be at least confidence_level. Returned in the same format as
        sim_number_positive (with the same list levels).
    """
    # save results in likelihood_circulation_all_scenarios_sims
    likelihood_circulation_all_scenarios_sims = [[[[[None]*len(sim_output_list[0][1][0][0])
                                                    for m in range(len(sim_output_list[0][1][0]))]
                                                   for n in range(len(sim_output_list[0][1]))]
                                                  for o in range(len(sim_output_list[0]))]
                                                 for p in range(len(sim_output_list))]

    # surveillance scenario
    for i1 in range(len(sim_output_list)):
        # diagnostic test
        for i3 in range(len(sim_output_list[i1][1][i3])):
            # transmission intensity
            for i4 in range(len(sim_output_list[i1][1][i3][i4])):
                # simulation run
                for i5 in range(len(sim_output_list[i1][1][i3][i4][i5])):
                    likelihood_circulation_all_scenarios_sims[i1][1][i3][i4][i5] = \
                        likelihood_circulation_all_negative(pop_size=pop_size,
                                                            number_tested=number_tested,
                                                            number_pos_circulation=sim_output_list[i1][1][i3][i4][i5],
                                                            number_pos_no_circulation=np.random.binomial(n=pop_size,
                                                                                                         p=false_pos_prob[i3]),
                                                            prob_circulation=prob_circulation)

    return likelihood_circulation_all_scenarios_sims



