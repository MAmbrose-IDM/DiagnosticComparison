# import numpy as np
import math
from operator import mul
from functools import reduce


def prob_detected_without_replacement(pop_size, number_pos, number_tested, times_sampled):
    """
    Calculate the probability that at least one individual will test positive if number_tested individuals are sampled
        during each surveillance period in a population where number_pos individuals would test positive out of a total
        population of pop_size individuals. Assume sampling without replacement.

    Arguments:
        pop_size - number of individuals in the sampled population
        number_pos - number of individuals who would have a positive diagnostic test
        number_tested - number of individuals tested
        times_sampled - number of times a year the population is surveyed (assume surveys are independent from one
             another)
    """
    prob_positive = 1 - ((reduce(mul, list(range((pop_size - number_pos - number_tested + 1),
                                                 (pop_size - number_pos + 1))))
                         / reduce(mul, list(range((pop_size - number_tested + 1), (pop_size + 1)))))
                         ** times_sampled)

    # Previous version used factorials
    # Note: the factorial version doesn't work when pop_size is large.
    # current_prob = 1 - ((math.factorial(pop_size - number_pos_samples[ii])
    #                    / math.factorial(pop_size - number_pos_samples[ii] - upper_value)
    #                    * math.factorial(pop_size - upper_value)
    #                    / math.factorial(pop_size)) ** times_sampled)
    return prob_positive


def calculate_num_samples_needed(pop_size=500, num_test_positive=25, confidence_level=0.95):
    """

    Given information about the number of individuals who would test positive on a given surveillance date and the
        total population size, return the number of individuals that would need to be sampled to result in a certain
        probability of detecting the infection. Assume sampling without replacement (so we're calculating the number of
        unique individuals that should be sampled).

    :param pop_size: number of individuals in surveilled population
    :param num_test_positive: number of individuals that would test positive for this diagnostic test if sampled
    :param confidence_level: if circulation is occuring, how certain do we want to be that our sampling design will
        detect it?
    :return: the number of samples that need to be collected from this population if we want to be confidence_level
        sure that we will detect positive samples if the disease is circulating
    """

    # Check whether all samples are positive - if so, can just sample one individual
    if num_test_positive == pop_size:
        return 1
    # Check whether no samples are positive - if so, no number of samples will be sufficient
    if num_test_positive == 0:
        return math.inf

    # Find the minimum number of samples so that probability of detection is above confidence_level
    upper_value = min((pop_size - num_test_positive), math.ceil(math.log(1 - confidence_level)
                                                                / math.log(1 - num_test_positive/pop_size)))
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

    return upper_value + 1  # upper.value is one below the threshold probability


def calculate_num_samples_each_sim(sim_number_positive, pop_size=500, confidence_level=0.95):
    """

    :param sim_number_positive: number of individuals whose test would be positive if they were sampled on the
        surveillance day. It is a list of list of lists. In this object,
        the outer-most level has one entry for each of the diagnostic tests; the middle level has one
        entry for each scenario on the number of cases in a year (from num_case_in_year_list); and the inner-most level
        has one entry for each simulation run.
    :param pop_size:  number of individuals in surveilled population
    :param confidence_level:  if circulation is occuring, how certain do we want to be that our sampling design will
        detect it?
    :return: the number of samples that would need to be collected for each simulated population so that the probability
        of detecting a positive test would be at least confidence_level. Returned in the same format as
        sim_number_positive (with the same list levels).
    """
    # save results in sim_samples_needed
    sim_samples_needed = [[[None]*len(sim_number_positive[0][0]) for n in range(len(sim_number_positive[0]))]
                             for m in range(len(sim_number_positive))]
    for i1 in range(len(sim_number_positive)):
        for i2 in range(len(sim_number_positive[i1])):
            for i3 in range(len(sim_number_positive[i1][i2])):
                sim_samples_needed[i1][i2][i3] = calculate_num_samples_needed(pop_size=pop_size,
                                                                              num_test_positive=sim_number_positive[i1][i2][i3],
                                                                              confidence_level=confidence_level)
    return sim_samples_needed

# # test function
# print(calculate_num_samples_needed())
