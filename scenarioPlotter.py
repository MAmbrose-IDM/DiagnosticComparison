# Goal of this script is to read in results from malaria simulations and use them to assess the sampling effort
#   required to detect malaria circulation in a village under different circumstances.
# Sampling effort can involve the number of individuals sampled as well as the number and timing of sampling trips.
# Transmission circumstances can involve transmission intensity, seasonality, new introductions from
#    outside the village, recent fade-out from village
# In addition, test sensitivity and behavior may vary by type and LOD of particular test

import numpy as np
import math
import matplotlib.pyplot as plt
# from operator import mul
# from functools import reduce
from matplotlib import cm

from calculateNumSamplesNeeded import prob_detected_without_replacement


def samples_needed_each_test_plotter(x_values, values_for_plot, pop_size, diagnostic_names, confidence_level,
                                     num_case_in_year_list, ax=None):
    """
    :param x_values: x-coordinates for lines that are to be plotted
    :param values_for_plot: points to corresponding to the y-coordinates of lines. Object is a list of lists of lists.
        The outer list contains plotting information for each of the diagnostic tests. The middle level contains three
        lists: one for the mean values, one for the upper 95%CI and one for the lower 95%CI bounds.
    :param pop_size: number of individuals in the population
    :param diagnostic_names: list of character strings naming the diagnostic tests plotted
    :param confidence_level: level of certainty that no positive samples means no circulation
    :param num_case_in_year_list: list containing values used for simulations under different transmission
        intensities. Gives the average number of malaria cases expected to occur in a year among all individuals in the
        population (assumed to be evenly distributed through time and assume individuals are selected to be infected at
        random)
    :param ax: axes to use for the current plot
    """
    if ax is None:
        ax = plt.gca()

    colormap0 = cm.tab10.colors
    color_order = [5, 9, 8, 7, 6, 2, 1, 0, 4, 8, 3]
    colormap = [colormap0[y] for y in color_order]

    for pp in range(len(values_for_plot)):
        ax.plot(x_values, values_for_plot[pp][0], 'k', color=colormap[pp], linewidth=3.0,  label=diagnostic_names[pp])
        ax.fill_between(x_values, values_for_plot[pp][1], values_for_plot[pp][2],
                        alpha=0.25, facecolor=colormap[pp])

    # ax = plt.subplot(111)
    # plot the color for each circulation intensity along the x-axis
    # get colors
    start = 0.0
    stop = 1.0
    number_of_lines = len(num_case_in_year_list)
    cm_subsection = np.linspace(start, stop, number_of_lines)

    colors = [cm.jet(x) for x in cm_subsection]
    ax.scatter(num_case_in_year_list, [0]*len(num_case_in_year_list), c=colors, s=50, marker="|")

    # ax.set_ylim([0, 800])
    ax.set_ylabel('number of samples needed')
    ax.set_xlabel('number of cases per year')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.set_title('Population size = %i; Confidence level = %.2f' % (pop_size, round(confidence_level, 2)), y=1.08)
    # plt.legend(title='Diagnostic method')



def diagnostic_probability_plotter(test_detection_probs, diagnostic_names, ax=None):
    """

    :param test_detection_probs: for each of the tests, gives the probabilities an individual who was infected X
        days ago will have a positive test result
    :param diagnostic_names: strings to use to describe each test
    :param ax: axes to use for the current plot
    """
    if ax is None:
        ax = plt.gca()
    # ax=plt.subplot(111)

    colormap0 = cm.tab10.colors
    color_order = [5, 9, 8, 7, 6, 2, 1, 0, 4, 8, 3]
    colormap = [colormap0[y] for y in color_order]

    for pp in range(len(test_detection_probs)):
        ax.plot(test_detection_probs[pp], 'k', color=colormap[pp], linewidth=4.0, alpha=0.75,
                label=diagnostic_names[pp])

    ax.set_ylabel('probability individual tests positive')
    ax.set_xlabel('days since most recent infection')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # plt.legend(title='Diagnostic method')



def time_since_infection_plotter(pop_size, num_case_in_year_list, sim_time_since_last_infection, ax=None):
    """
    Create plot of the distribution of times-since-infection under different transmission intensity scenarios

    :param pop_size: number of individuals in surveilled population
    :param num_case_in_year_list: list containing values used for simulations under different transmission
        intensities. Gives the average number of malaria cases expected to occur in a year among all individuals in the
        population (assumed to be evenly distributed through time and assume individuals are selected to be infected at
        random)
    :param sim_time_since_last_infection: the first element output from the simulate_over_transmission_intensities
        function. It describes the time since the last infection for individuals in the population.
        It is a list of lists of lists. The outer-most list is the transmission intensity scenario (each entry
        corresponds to a value from num_case_in_year_list);
        the middle list is the simulation; the inner list is the individuals in the population
    :param ax: axes to use for the current plot
    """
    if ax is None:
        ax = plt.gca()
    # fig = plt.figure()
    # ax = fig.add_subplot(111)  # plt.subplot(111)

    # get colors
    start = 0.0
    stop = 1.0
    number_of_lines = len(num_case_in_year_list)
    cm_subsection = np.linspace(start, stop, number_of_lines)

    colors = [cm.jet(x) for x in cm_subsection]

    num_data_points = len(sim_time_since_last_infection[0]) * len(sim_time_since_last_infection[0][0])

    for i1 in range(len(sim_time_since_last_infection)):
        # Across all simulations and individuals in the population, plot the distribution of time-since infection for
        #   a particular transmission intensity
        density_days_since_infection = [0] * 365
        density_more_than_364 = 0
        for i2 in range(len(sim_time_since_last_infection[i1])):
            for i3 in range(len(sim_time_since_last_infection[i1][i2])):
                if sim_time_since_last_infection[i1][i2][i3] >= 364:
                    density_more_than_364 += 1
                else:
                    density_days_since_infection[sim_time_since_last_infection[i1][i2][i3]] += 1
        ax.plot([y/num_data_points for y in density_days_since_infection], 'k', label=num_case_in_year_list[i1], color=colors[i1])

    # ax.set_title('Population size = %i' % pop_size, y=1.08)
    ax.set_ylabel('proportion of population')
    ax.set_xlabel('days since most recent infection')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # plt.legend(title='Number of cases per year')
    return ax.legend


def case_seasonality_plotter(pop_size, seasonal_scalar, num_case_in_year_list, surveillance_days, ax=None):
    """

    :param pop_size: number of individuals in surveilled population
    :param seasonal_scalar: relative intensity of transmission on this day of the year.
    :param num_case_in_year_list: list detailing the average number of malaria cases expected to occur in a year among
    all individuals in the population (list may contain several transmission intensity scenarios to plot)
    :param surveillance_days: list containing the day of the year surveillance was conducted
    :param ax: axes to use for the current plot
    """
    if ax is None:
        ax = plt.gca()

    # get colors
    start = 0.0
    stop = 1.0
    number_of_lines = len(num_case_in_year_list)
    cm_subsection = np.linspace(start, stop, number_of_lines)

    colors = [cm.jet(x) for x in cm_subsection]

    for i1 in range(len(num_case_in_year_list)):
        # Plot the expected number of cases on each day of the year for a particular transmission intensity
        ax.plot([num_case_in_year_list[i1]/365 * y for y in seasonal_scalar], 'k', label=num_case_in_year_list[i1],
                color=colors[i1])

    ax.axvline(x=surveillance_days[0], color='k')
    # ax.set_title('Population size = %i' % pop_size, y=1.08)
    ax.set_ylabel('expected number of cases')
    ax.set_xlabel('day of year')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # plt.legend(title='Number of cases per year')


def num_cases_legend_plotter(num_case_in_year_list, ax):
    """

    :param num_case_in_year_list: list detailing the average number of malaria cases expected to occur in a year among
    all individuals in the population (list may contain several transmission intensity scenarios to plot)
    :param ax: axes to use for the current plot
    """
    # get colors
    start = 0.0
    stop = 1.0
    number_of_lines = len(num_case_in_year_list)
    cm_subsection = np.linspace(start, stop, number_of_lines)

    colors = [cm.jet(x) for x in cm_subsection]

    ax.axis('off')
    ax.set_xlim([0, 4])
    ax.set_ylim([-0.5, len(num_case_in_year_list)])
    ax.text(0.5, len(num_case_in_year_list), 'Number of cases per year', horizontalalignment='left')
    # ax.text(0.5, len(num_case_in_year_list), 'Number of cases per year per capita', horizontalalignment='left')

    for i1 in range(len(num_case_in_year_list)):
        ax.plot([1, 2], [i1, i1], 'k', label=num_case_in_year_list[i1],
                color=colors[i1])
        ax.text(2.1, i1, num_case_in_year_list[i1], verticalalignment='center')


def diagnostic_legend_plotter(diagnostic_names, ax):
    """

    :param diagnostic_names: strings to use to describe each test
    :param ax: axes to use for the current plot
    """

    colormap0 = cm.tab10.colors
    color_order = [5, 9, 8, 7, 6, 2, 1, 0, 4, 8, 3]
    colormap = [colormap0[y] for y in color_order]

    ax.axis('off')
    ax.set_xlim([0, 4])
    ax.set_ylim([-0.5, len(diagnostic_names)])
    ax.text(0.5, len(diagnostic_names), 'Diagnostic test', horizontalalignment='left')

    for pp in range(len(diagnostic_names)):
        ax.plot([1, 2], [pp, pp], 'k', color=colormap[pp], linewidth=4.0, alpha=0.75)
        ax.text(2.1, pp, diagnostic_names[pp], verticalalignment='center')


# Older plotting functions:

def samples_needed_fraction_positive_plotter():
    """
    Plot of the number of samples that should be collected so that if transmission is occurring, there is a 75%, 85%,
    or 95% chance that an infected sample will be collected. Here, we assume large populations or sampling with
    replacement.
    """

    prob_detect_list = [0.95, 0.85, 0.75]  # limit for probability we will detect transmission
    line_color_list = ['k-', 'b-', 'g-']  # colors to use in plot for each detection probability
    fraction_pos_samples = np.linspace(0.01, 0.4, 100)  # fraction of samples that test positive

    # create lines on plot corresponding to each detection probability
    for jj in range(len(prob_detect_list)):
        plt.plot(fraction_pos_samples,
                 [math.log(1 - prob_detect_list[jj]) / math.log(1 - y) for y in fraction_pos_samples],
                 line_color_list[jj], label=prob_detect_list[jj]
                 )

    plt.ylabel('number of samples needed')
    plt.xlabel('fraction of individuals testing positive')
    plt.legend(title='Probability of detecting ongoing transmission')
    plt.show()


def samples_needed_number_positive_plotter(pop_size=100, times_sampled=1, prob_detect_list=None,
                                           line_color_list=None,
                                           fraction_pos_samples=None, ax=None, subplot=True):
    """
    Plot of the number of samples that should be collected so that if transmission is occurring, there is a 75%, 85%,
    or 95% chance (by default, different values can be provided as an argument) that an infected sample will be
    collected. Here, we assume that the population is of size pop_size and sampling is done without replacement. We
    assume that sampling efforts are carried out times_sampled time in a year.

    Arguments:
        pop_size - number of individuals in the sampled population
        times_sampled - number of times a year the population is surveyed (Assume surveys independent from one another)
        prob_detect_list - limit for probability we will detect transmission
        line_color_list - colors to use in plot for each detection probability
        fraction_pos_samples - fraction of samples that test positive if the entire population were surveyed
        ax - the axes in the subplot to create this figure
        subplot - flag that indicates whether a plot should be created upon function call or whether the plot should
                  be returned as a subplot in a multi-paneled figure
    """

    # insert default values if none are provided:
    if prob_detect_list is None:
        prob_detect_list = [0.95, 0.85, 0.75]  #
    if line_color_list is None:
        line_color_list = ['k-', 'b-', 'g-']  # colors to use in plot for each detection probability
    if fraction_pos_samples is None:
        fraction_pos_samples = np.linspace(0.01, 0.4, 100)  # fraction of samples that test positive
    if ax is None:
        ax = plt.gca()

    number_pos_samples = [math.floor(pop_size * y) for y in fraction_pos_samples]  # number that test positive

    # For each value in number_pos_samples, find the minimum number of samples so that probability of detection is
    #   above each of the prob_detect_list limits.
    # Store results in a list of lists (contains one list for each prob_detect_list value)
    # Assumes prob_detect_list is sorted from largest to smallest value
    min_samples_needed = [[None]*len(number_pos_samples) for n in range(len(prob_detect_list))]
    for ii in range(len(number_pos_samples)):
        upper_value = min((pop_size - number_pos_samples[ii]), math.ceil(math.log(1 - prob_detect_list[0])
                                                                         / math.log(1 - fraction_pos_samples[ii])))
        # probability infection would be detected if upper.value individuals were sampled
        if upper_value == 0:
            current_prob = 0
        else:
            current_prob = prob_detected_without_replacement(pop_size=pop_size,
                                                             number_pos=number_pos_samples[ii],
                                                             number_tested=upper_value,
                                                             times_sampled=times_sampled
                                                             )

        # currently coded to search each value; if population sizes are large, a more efficient search algorithm would
        #   be preferable
        for jj in range(len(prob_detect_list)):
            while current_prob > prob_detect_list[jj]:
                # decrease the number of samples taken and recalculate the new current_prob
                upper_value -= 1
                if upper_value == 0:
                    current_prob = 0
                else:
                    current_prob = prob_detected_without_replacement(pop_size=pop_size,
                                                                     number_pos=number_pos_samples[ii],
                                                                     number_tested=upper_value,
                                                                     times_sampled=times_sampled
                                                                     )

            min_samples_needed[jj][ii] = upper_value + 1  # upper.value is one below the threshold probability

    if subplot:
        # create lines on plot corresponding to each detection probability
        for jj in range(len(prob_detect_list)):
            cur_plot = ax.plot(fraction_pos_samples,
                               min_samples_needed[jj],
                               line_color_list[jj], label=prob_detect_list[jj]
                               )
        # vertical line showing population size
        # plt.plot(fraction_pos_samples, [pop_size]*len(fraction_pos_samples),
        #          '--r', label='population size'
        #          )
        ax.set_ylabel('number of samples needed')
        ax.set_xlabel('fraction of population positive')
        ax.set_ylim([0, 200])
        # ax.set_xlim([0, max(fraction_pos_samples)])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # plt.legend(title='Probability of detecting ongoing transmission')
        ax.set_title('N = %i ; Sampling times = %i' % (pop_size, times_sampled), y=1.08)
    else:
        # create lines on plot corresponding to each detection probability
        for jj in range(len(prob_detect_list)):
            plt.plot(fraction_pos_samples,
                     min_samples_needed[jj],
                     line_color_list[jj], label=prob_detect_list[jj]
                     )
        # vertical line showing population size
        # plt.plot(fraction_pos_samples, [pop_size]*len(fraction_pos_samples),
        #          '--r', label='population size'
        #          )
        plt.ylabel('number of samples needed')
        plt.xlabel('fraction of population positive')
        # plt.legend(title='Probability of detecting ongoing transmission')
        plt.title('Number of samples needed to detect ongoing transmission \n Pop size = %i ; Sampling %i times a year'
                  % (pop_size, times_sampled))
        plt.show()
        cur_plot = None
    return cur_plot


def prob_positive_sample_plotter(pop_size=100, times_sampled=1, fraction_pos_samples=None,
                                 line_color_list=None, numbers_sampled_list=None,
                                 ax=None, subplot=True):
    """
    Plot of the probability a positive sample will be observed given the number of individuals surveyed during each
    visit. Different lines correspond to different fractions of the population positive (representing transmission
    intensity). Here, we assume that the population is of size pop_size and sampling is done without replacement. We
    assume that sampling efforts are carried out times_sampled time in a year.

    Arguments:
        pop_size - number of individuals in the sampled population
        times_sampled - number of times a year the population is surveyed (Assume surveys independent from one another)
        fraction_pos_samples - fraction of samples that test positive if the entire population were surveyed
        line_color_list - colors to use in plot for each fraction of samples that test positive
        numbers_sampled_list - number of individuals sampled at each visit
        ax - the axes in the subplot to create this figure
        subplot - flag that indicates whether a plot should be created upon function call or whether the plot should
                  be returned as a subplot in a multi-paneled figure
    """

    # insert default values if none are provided:
    if line_color_list is None:
        # line_color_list = ['p-', 'r-', 'o-', 'g-']  # colors to use in plot for each detection probability
        line_color_list = ['darkblue', 'green', 'peru', 'maroon']  # colors to use in plot for each detection probability
    if fraction_pos_samples is None:
        fraction_pos_samples = [0.01, 0.05, 0.1, 0.2]  # fraction of samples that test positive
    if numbers_sampled_list is None:
        numbers_sampled_list = list(range(1, 101))
    if ax is None:
        ax = plt.gca()

    number_pos_samples = [math.floor(pop_size * y) for y in fraction_pos_samples]  # number that test positive

    # For each value in fraction_pos_samples, find the probability a positive sample is detected if
    # numbers_sampled_list samples are collected. Store results in a list of lists (contains one list for each
    # fraction_pos_samples value).
    prob_detected = [[None]*len(numbers_sampled_list) for n in range(len(number_pos_samples))]
    for ii in range(len(number_pos_samples)):
        for jj in range(len(numbers_sampled_list)):
            prob_detected[ii][jj] = prob_detected_without_replacement(pop_size=pop_size,
                                                                      number_pos=number_pos_samples[ii],
                                                                      number_tested=numbers_sampled_list[jj],
                                                                      times_sampled=times_sampled
                                                                      )

    if subplot:
        # create lines on plot corresponding to each detection probability
        for jj in range(len(number_pos_samples)):
            cur_plot = ax.plot(numbers_sampled_list,
                               prob_detected[jj],
                               line_color_list[jj], label=fraction_pos_samples[jj]
                               )

        ax.set_ylabel('probability of a positive test')
        ax.set_xlabel('number of individuals tested')
        ax.set_ylim([0, 1])
        # ax.set_xlim([0, max(fraction_pos_samples)])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # plt.legend(title='Probability of detecting ongoing transmission')
        ax.set_title('N = %i ; Sampling times = %i' % (pop_size, times_sampled), y=1.08)
    else:
        # create lines on plot corresponding to each detection probability
        for jj in range(len(number_pos_samples)):
            plt.plot(numbers_sampled_list,
                     prob_detected[jj],
                     line_color_list[jj], label=fraction_pos_samples[jj]
                     )
        # vertical line showing population size
        # plt.plot(fraction_pos_samples, [pop_size]*len(fraction_pos_samples),
        #          '--r', label='population size'
        #          )
        plt.ylabel('probability of a positive test')
        plt.xlabel('number of individuals tested')
        # plt.legend(title='Probability of detecting ongoing transmission')
        plt.title('Probability of detecting ongoing transmission as a function of number of samples collected \n Pop size = %i ; Sampling %i times a year'
                  % (pop_size, times_sampled))
        plt.show()
        cur_plot = None
    return cur_plot






# test calling the functions:
# samples_needed_fraction_positive_plotter()
# samples_needed_number_positive_plotter(pop_size=1000, times_sampled=3, fraction_pos_samples=np.linspace(0.01, 0.1, 20))

# # Create figure panels
#
# # Number of samples needed to detect ongoing transmission
# # create multiple subplots across a range of population sizes and times sampled in a year
# pop_size_list = [100, 500, 1000]
# times_sampled_list = [1, 2, 3]
#
# fig, ax_array = plt.subplots(nrows=len(pop_size_list), ncols=len(times_sampled_list), figsize=(10, 10))
# fig.subplots_adjust(hspace=0.7, wspace=0.7)
# for rr, ax_row in enumerate(ax_array):
#     for cc, axes in enumerate(ax_row):
#         samples_needed_number_positive_plotter(pop_size=pop_size_list[rr], times_sampled=times_sampled_list[cc],
#                                                fraction_pos_samples=np.linspace(0.001, 0.1, 20), ax=axes)
# # axes.ylabel('number of samples needed')
# # axes.xlabel('fraction of individuals testing positive')
# fig.suptitle('Number of samples needed to detect ongoing transmission')
#
# # fig.legend(title='Probability of detecting ongoing transmission')
# plt.show()
#
#
# # Probability of detecting ongoing transmission
# # create multiple subplots across a range of population sizes and times sampled in a year
# pop_size_list = [100, 500, 1000]
# times_sampled_list = [1, 2, 3]
#
# fig, ax_array = plt.subplots(nrows=len(pop_size_list), ncols=len(times_sampled_list), figsize=(10, 10))
# fig.subplots_adjust(hspace=0.7, wspace=0.7)
# for rr, ax_row in enumerate(ax_array):
#     for cc, axes in enumerate(ax_row):
#         prob_positive_sample_plotter(pop_size=pop_size_list[rr], times_sampled=times_sampled_list[cc], ax=axes)
# # axes.ylabel('probability detected')
# # axes.xlabel('number of individuals tested')
# fig.suptitle('Probability of detecting ongoing transmission give number of individuals tested')
#
# # fig.legend(title='Probability of detecting ongoing transmission')
# plt.show()