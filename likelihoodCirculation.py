# import numpy as np
from operator import mul
from functools import reduce


def likelihood_circulation_all_negative(pop_size, number_tested, number_pos_circulation, number_pos_no_circulation,
                              prob_circulation):
    """
    Calculate the probability that at least one individual will test positive if number_tested individuals are sampled
        during each surveillance period in a population where number_pos individuals would test positive out of a total
        population of pop_size individuals. Assume sampling without replacement.

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

    likelihood_circulation = (prob_all_negative_given_circulation * prob_circulation) \
                             / ((prob_all_negative_given_no_circulation * (1 - prob_circulation))
                                + (prob_all_negative_given_circulation * prob_circulation))

    return likelihood_circulation
