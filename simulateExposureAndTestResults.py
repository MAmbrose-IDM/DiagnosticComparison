# Create a dataset of the dates of most recent exposure for all individuals in a population, assuming a certain rate of
#  new cases and assuming that new infections are selected randomly from the population (no role of age, immunity
#  status, etc. Let num_case_in_year gives the number of new cases that occur each year. If num_case_in_year < 365,
#  select one individual randomly from the population to be infected every num_case_in_year/365 days. If
#  num_case_in_year>365, select an average of num_case_in_year/365 individuals every day. Selected individuals have
#  their date of most recent infection updated. After burn-in, record the number of days since the most recent
#  infection for all individuals in the population on 'surveillance' days.

# In current iteration, not at all mechanistic, simply assumes individuals' infections are distributed randomly
#  throughout the year and ignores immunity from previous exposure.

import numpy as np


def time_since_infection(pop_size=500, num_case_in_year=25, surveillance_days=None, initialize_years=5,
                         seasonal_scalar=None):
    """

    :param pop_size: number of individuals in surveilled population
    :param num_case_in_year: average number of malaria cases expected to occur in a year among all individuals in the
    population (assumed to be evenly distributed through time and assume individuals are selected to be infected at
    random)
    :param surveillance_days: list of days on which surveillance is carried out
    :param initialize_years: how many years should the simulations be run before we start collecting data
    :param seasonal_scalar: relative intensity of transmission on this day of the year.
    :return: list of lists. Outer list is for each surveillance day. For a given surveillance day, inner list contains
    the number of days since last infection for each individual in the population
    """

    if surveillance_days is None:
        surveillance_days = [365]
    if seasonal_scalar is None:
        seasonal_scalar = [1] * 365

    # list of lists giving the number of days since last infection for each individual (inner list) as observed on each
    #  surveillance day (outer list)
    sampled_days_since_infection = [[None]*pop_size for n in range(len(surveillance_days))]

    # list of current number of days since last infection (start at day 0, with everyone having 10000 days since last
    #  infection)
    day_of_last_infection = [-10000]*pop_size

    # total number of days to run will be 365 * (initialize_years + 1). Surveys happen in last year.
    total_days = 365 * (initialize_years + 1)
    survey_day_num = [sd + 365*initialize_years for sd in surveillance_days]

    # average number of new cases on each day (this part will be stochastic as we select how many individuals for
    #  a given day)
    # old version with constant intensity through time: ave_num_cases_each_day = num_case_in_year/365
    # new version has a different expected number for each day based on seasonality
    ave_num_cases_each_day = [num_case_in_year/365/pop_size * y for y in seasonal_scalar] * (initialize_years + 1)

    # run through days, updating current_days_since_infection where relevant. If it is a surveillance day, record
    #  current values
    for dd in range(1, (total_days + 1)):
        # draw number of cases to add on this day
        # old version with constant intensity through time: num_new_cases = np.random.binomial(n=pop_size,
        #             p=ave_num_cases_each_day/pop_size)
        num_new_cases = np.random.binomial(n=pop_size, p=ave_num_cases_each_day[dd-1])
        # check whether anyone was infected this day
        if num_new_cases > 0:
            # determine which individuals get infected
            for ii in range(num_new_cases):
                day_of_last_infection[np.random.randint(low=0, high=pop_size)] = dd
            # is dd one of the surveillance days?
        if survey_day_num.count(dd) > 0:
            survey_number = survey_day_num.index(dd)
            sampled_days_since_infection[survey_number] = [dd - y for y in day_of_last_infection]

    return sampled_days_since_infection


def get_test_results(days_since_infection, test_detection_probs=None, false_pos_prob=None):
    """
    assign test results to each individual at each surveillance point based on number of days since last infection and
     on sensitivity/specificity of each test

    :param days_since_infection: list giving the number of days since last infection for all individuals in the
    population (one of the inner lists returned from time_since_infection)
    :param test_detection_probs: for each test (outer list) what is probability of a positive result if it has been ii
    days since last infection, where ii=0 is the first index in the inner list
    :param false_pos_prob: probability an individual has a positive test independent of exposure history
    :return: list of lists: outer list is for each diagnostic test in test_detection_probs, inner list gives whether
    each individual was sampled to test positive (1) or negative (0) for that test
    """
    if test_detection_probs is None:
        test_detection_probs = [[0.0]*20 + [0.8]*40 + [0.2]*10 + [0.0]*20]
    if false_pos_prob is None:
        false_pos_prob = [0.0]*len(test_detection_probs)

    store_test_results = [[None]*len(days_since_infection) for n in range(len(test_detection_probs))]

    for tt in range(len(test_detection_probs)):
        cur_test_probs = test_detection_probs[tt]
        for ii in range(len(days_since_infection)):
            if days_since_infection[ii] >= len(cur_test_probs):
                cur_prob = 0.0
            else:
                cur_prob = cur_test_probs[days_since_infection[ii]]
            # could test positive from actual infection or from false positive
            prob_pos = 1.0 - (1.0-cur_prob) * (1.0-false_pos_prob[tt])
            store_test_results[tt][ii] = np.random.binomial(n=1, p=prob_pos)
    return store_test_results


def simulate_over_transmission_intensities(pop_size=500, num_case_in_year_list=None,
                                           initialize_years=5, seasonal_scalar=None, surveillance_days=None,
                                           test_detection_probs=None, false_pos_prob=None, num_sims_each=1):
    """

    Iterate over values in num_case_in_year_list (each representing a different transmission intensity). For each
    value, run num_sims_each simulations to get the time since last infection for each individual in the population.
    Then use test_detection_probs to establish whether each individual would test positive or negative at the
    surveillance time points indicated in surveillance_days.

    :param pop_size: number of individuals in surveilled population
    :param num_case_in_year_list: list containing values to use for simulations under different transmission
        intensities. Gives the average number of malaria cases expected to occur in a year among all individuals in the
        population (assumed to be evenly distributed through time and assume individuals are selected to be infected at
        random)
    :param initialize_years: how many years should the simulations be run before we start collecting data
    :param seasonal_scalar: relative intensity of transmission on this day of the year.
    :param surveillance_days: list giving day for survillance (currently, must send a list but only first element used)
    :param test_detection_probs: for each test (outer list) what is probability of a positive result if it has been ii
        days since last infection, where ii=0 is the first index in the inner list
    :param false_pos_prob: probability an individual has a positive test independent of exposure history
    :param num_sims_each: number of simulations to run for each of the transmission intensities (each of the values in
        num_case_in_year_list)
    :return: a list of two items. The first item gives the time since the last infection for individuals
        in the population. It is a list of lists of lists. The outer-most list is the transmission intensity scenario
        ( the value used from num_case_in_year_list); the middle list is the simulation; the inner list is the
        individuals in the population. The second item gives the number of individuals in the population that would
        test positive if sampled on the surveillance day. It is a list of list of lists. In this object,
        the outer-most level has one entry for each of the diagnostic tests; the middle level has one
        entry for each scenario on the number of cases in a year (from num_case_in_year_list); and the inner-most level
        has one entry for each simulation run.
    """
    if num_case_in_year_list is None:
        num_case_in_year_list = [5, 25, 100]
    if test_detection_probs is None:
        test_detection_probs = [[0.0]*20 + [0.8]*40 + [0.2]*10 + [0.0]*20]
    if false_pos_prob is None:
        false_pos_prob = [0.0]*len(test_detection_probs)
    if surveillance_days is None:
        surveillance_days = [365]
    if seasonal_scalar is None:
        seasonal_scalar = [1] * 365

    # Iterate through transmission scenarios, running num_sims_each simulations for each.
    # Create a list of list of lists variable to store the time since last infection for all
    #    1) transmission intensities, 2) simulations, and 3) individuals.
    sim_time_since_last_infection = [[[None]*pop_size for n in range(num_sims_each)]
                                     for m in range(len(num_case_in_year_list))]
    # Create a list of lists of lists variable to store the fraction of tests positive in the population for all
    #     1) diagnostic tests, 2) transmission intensities, and 3) simulations
    sim_number_positive = [[[None]*num_sims_each for n in range(len(num_case_in_year_list))]
                             for m in range(len(test_detection_probs))]
    for ii in range(len(num_case_in_year_list)):  # iterate through transmission intensities
        print('Simulating transmission intensity number: %i' % ii)
        for ss in range(num_sims_each):  # iterate through simulations within a transmission intensity

            sim_time_since_last_infection[ii][ss] \
                = (time_since_infection(pop_size=pop_size, num_case_in_year=num_case_in_year_list[ii],
                                        surveillance_days=surveillance_days, initialize_years=initialize_years,
                                        seasonal_scalar=seasonal_scalar))[0]
            sim_fraction_positive_all_tests \
                = get_test_results(days_since_infection=sim_time_since_last_infection[ii][ss],
                                   test_detection_probs=test_detection_probs, false_pos_prob=false_pos_prob)
            for dd in range(len(test_detection_probs)):
                sim_number_positive[dd][ii][ss] = np.sum(sim_fraction_positive_all_tests[dd])

    return [sim_time_since_last_infection, sim_number_positive]


# # test the functions
# output = time_since_infection(pop_size=500, num_case_in_year=50, surveillance_days=[365], initialize_years=5)
# print('number of infections within last 30 days of survey 1:')
# print(sum([y < 30 for y in output[0]]))
# # specify detection probabilities for each of the tests (index corresponds to number of days after infection)
# test_detection_probs = [([0]*2 + list(np.linspace(0.001, 0.7, 5)) + list(np.linspace(0.7, 0.5, 20))
#                          + list(np.linspace(0.5, 0, 20))),
#                         ([0]*2 + list(np.linspace(0.001, 0.9, 5)) + list(np.linspace(0.9, 0.7, 50))
#                          + list(np.linspace(0.7, 0, 30))),
#                         ([0]*2 + list(np.linspace(0.001, 0.9, 14)) + list(np.linspace(0.9, 0.5, 100))
#                          + list(np.linspace(0.5, 0, 120)))]
#
# output2 = get_test_results(days_since_infection=output[0], test_detection_probs=test_detection_probs)
# print('number of tests positive for each test type:')
# for ii in range(len(output2)):
#     print(sum(output2[ii]))
#
# print((simulate_over_transmission_intensities(num_case_in_year_list=[50, 100, 600],
#                                               test_detection_probs=test_detection_probs))[1][0][2][0])


