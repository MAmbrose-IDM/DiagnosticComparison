import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


from plotters_fromDTKSims import prob_dist_n_p_plotter, \
    prob_positive_sample_plotter, \
    likelihood_given_observation_plotter



#
# Probability distribution for number of individuals positive given circulation
#
def plot_panel_prob_num_pos_circulation(pop_size, test_names, filename_suffix, all_sampling_dates, all_LHs):
    """

    :param pop_size: population sized assumed for plots
    :param test_names: names of diagnostic tests
    :param filename_suffix: indicates which DTK simulation set to use
    :param all_sampling_dates: list of dates when sampling occurs
    :param all_LHs: list of larval habitat multipliers used in simulations
    """
    # Set up the axes with gridspec
    fig = plt.figure(figsize=(14, 10.5))
    # grid = plt.GridSpec(4, 4, hspace=2, wspace=1)

    num_sampling_dates = len(all_sampling_dates)  # will be the columns
    num_LHs = len(all_LHs)  # will be the rows
    grid = gridspec.GridSpec((num_LHs+1), (num_sampling_dates+1), hspace=1, wspace=0.6,
                             width_ratios=[0.1] + [2]*num_sampling_dates,
                             height_ratios=[0.1] + [2]*num_LHs
                             )

    for s1 in range(len(all_sampling_dates)):
        for s2 in range(len(all_LHs)):
            # load files describing the probability of each possible number of positive individuals for this scenario
            with open("simOutputs_DTK/prob_num_pos_circulation_samplingDate%i_xLH%i_%s.p" % (all_sampling_dates[s1],
                                                                                             round(all_LHs[s2] * 100),
                                                                                             filename_suffix), "rb") as f:
                prob_num_pos_circulation = pickle.load(f)

            ax_cur = fig.add_subplot(grid[s2+1, s1+1])
            prob_dist_n_p_plotter(prob_num_pos=prob_num_pos_circulation,
                                  diagnostic_names=test_names, pop_size=pop_size, ax=ax_cur)

    # row and column labels
    for s1 in range(len(all_sampling_dates)):
        ax = fig.add_subplot(grid[0, s1+1])
        ax.axis('off')
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        ax.text(0.5, 0.5, 'sampling day: %i' % all_sampling_dates[s1], horizontalalignment='center',
                verticalalignment='center', fontweight='bold')
    for s2 in range(len(all_LHs)):
        ax = fig.add_subplot(grid[s2+1, 0])
        ax.axis('off')
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        ax.text(0, 0.5, 'LH multiplier: %.1f' % all_LHs[s2], horizontalalignment='center',
                verticalalignment='center', rotation=90, fontweight='bold')

    fig.suptitle('Probability distribution for number of individuals positive given circulation is occurring')

    plt.show()

    fig.savefig('figures/ProbDistNumPositive_circulation_DTK_%s.pdf' % filename_suffix)


#
# Probability distribution for number of individuals positive given no circulation
#
def plot_panel_prob_num_pos_no_circulation(pop_size, test_names, filename_suffix, all_sampling_dates, all_LHs):
    """

    :param pop_size: population sized assumed for plots
    :param test_names: names of diagnostic tests
    :param filename_suffix: indicates which DTK simulation set to use
    :param all_sampling_dates: list of dates when sampling occurs
    :param all_LHs: list of larval habitat multipliers used in simulations
    """
    # Plot looks the same across all scenarios - just show one
    fig = plt.figure(figsize=(7, 7))
    ax_cur = plt.subplot(111)
    # load files describing the probability of each possible number of positive individuals for this scenario

    with open("simOutputs_DTK/prob_num_pos_no_circulation_samplingDate%i_xLH%i_%s.p" % (all_sampling_dates[0],
                                                                                        round(all_LHs[0] * 100),
                                                                                        filename_suffix), "rb") as f:
        prob_num_pos_no_circulation = pickle.load(f)

    prob_dist_n_p_plotter(prob_num_pos=prob_num_pos_no_circulation,
                          diagnostic_names=test_names, pop_size=pop_size, ax=ax_cur, xmax=10)

    fig.suptitle('Probability distribution for number of individuals \n positive given circulation is not occurring')

    plt.show()

    fig.savefig('figures/ProbDistNumPositive_no_circulation_DTK_%s.pdf' % filename_suffix)


#
# Probability of at least one positive sample given circulation
#
def plot_panel_prob_positive_sample_circulation(pop_size, test_names, filename_suffix, all_sampling_dates, all_LHs):
    """

    :param pop_size: population sized assumed for plots
    :param test_names: names of diagnostic tests
    :param filename_suffix: indicates which DTK simulation set to use
    :param all_sampling_dates: list of dates when sampling occurs
    :param all_LHs: list of larval habitat multipliers used in simulations
    """
    # Set up the axes with gridspec
    fig = plt.figure(figsize=(14, 10.5))
    # grid = plt.GridSpec(4, 4, hspace=2, wspace=1)

    num_sampling_dates = len(all_sampling_dates)  # will be the columns
    num_LHs = len(all_LHs)  # will be the rows
    grid = gridspec.GridSpec((num_LHs+1), (num_sampling_dates+1), hspace=1, wspace=0.6,
                             width_ratios=[0.1] + [2]*num_sampling_dates,
                             height_ratios=[0.1] + [2]*num_LHs
                             )

    for s1 in range(len(all_sampling_dates)):
        for s2 in range(len(all_LHs)):
            # load files describing the probability of each possible number of positive individuals for this scenario
            with open("simOutputs_DTK/prob_pos_sample_given_circulation_samplingDate%i_xLH%i_%s.p" % (all_sampling_dates[s1],
                                                                                                      round(all_LHs[s2] * 100),
                                                                                                      filename_suffix), "rb") as f:
                prob_pos_sample_given_circulation = pickle.load(f)

            ax_cur = fig.add_subplot(grid[s2+1, s1+1])
            prob_positive_sample_plotter(prob_pos_sample=prob_pos_sample_given_circulation, diagnostic_names=test_names, pop_size=pop_size,
                                         circulating_string='circulation is',
                                         detection_limit=0.95, ax=ax_cur, line_flag=True)

    # row and column labels
    for s1 in range(len(all_sampling_dates)):
        ax = fig.add_subplot(grid[0, s1+1])
        ax.axis('off')
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        ax.text(0.5, 0.5, 'sampling day: %i' % all_sampling_dates[s1], horizontalalignment='center',
                verticalalignment='center', fontweight='bold')
    for s2 in range(len(all_LHs)):
        ax = fig.add_subplot(grid[s2+1, 0])
        ax.axis('off')
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        ax.text(0, 0.5, 'LH multiplier: %.1f' % all_LHs[s2], horizontalalignment='center',
                verticalalignment='center', rotation=90, fontweight='bold')

    fig.suptitle('Probability at least one positive sample will be detected given circulation is occurring')

    plt.show()

    fig.savefig('figures/ProbAtLeastOnePositive_circulation_DTK_%s.pdf' % filename_suffix)


#
# Probability of at least one positive sample given no circulation
#
def plot_panel_prob_positive_sample_no_circulation(pop_size, test_names, filename_suffix, all_sampling_dates, all_LHs):
    """

    :param pop_size: population sized assumed for plots
    :param test_names: names of diagnostic tests
    :param filename_suffix: indicates which DTK simulation set to use
    :param all_sampling_dates: list of dates when sampling occurs
    :param all_LHs: list of larval habitat multipliers used in simulations
    """
    # Plot looks the same across all scenarios - just show one
    fig = plt.figure(figsize=(7, 7))
    ax_cur = plt.subplot(111)
    # load files describing the probability of each possible number of positive individuals for this scenario

    with open("simOutputs_DTK/prob_pos_sample_given_no_circulation_samplingDate%i_xLH%i_%s.p" % (all_sampling_dates[0],
                                                                                                 round(all_LHs[0] * 100),
                                                                                                 filename_suffix),
              "rb") as f:
        prob_pos_sample_given_no_circulation = pickle.load(f)

    prob_positive_sample_plotter(prob_pos_sample=prob_pos_sample_given_no_circulation, diagnostic_names=test_names,
                                 pop_size=pop_size,
                                 circulating_string='circulation is not',
                                 detection_limit=0.95, ax=ax_cur, line_flag=False, xmax=pop_size)

    fig.suptitle('Probability of at least one positive sample \n  given circulation is not occurring')

    plt.show()

    fig.savefig('figures/ProbAtLeastOnePositive_no_circulation_DTK_%s.pdf' % filename_suffix)


#
# Likelihood of circulation given observation
#
def plot_panel_likelihood_circulation(pop_size, test_names, filename_suffix, all_sampling_dates,
                                      all_LHs, ss_values, p_circulation):
    """

    :param pop_size: population sized assumed for plots
    :param test_names: names of diagnostic tests
    :param filename_suffix: indicates which DTK simulation set to use
    :param all_sampling_dates: list of dates when sampling occurs
    :param all_LHs: list of larval habitat multipliers used in simulations
    :param ss_values: list of number of individuals sampled during survey
    :param p_circulation: list of values for the probability of circulation
    """

    for ss_index, ss in enumerate(ss_values):
        for p_c in range(len(p_circulation)):

            # Set up the axes with gridspec
            fig = plt.figure(figsize=(14, 10.5))
            # grid = plt.GridSpec(4, 4, hspace=2, wspace=1)

            num_sampling_dates = len(all_sampling_dates)  # will be the columns
            num_LHs = len(all_LHs)  # will be the rows
            grid = gridspec.GridSpec((num_LHs + 1), (num_sampling_dates + 1), hspace=1, wspace=0.6,
                                     width_ratios=[0.1] + [2] * num_sampling_dates,
                                     height_ratios=[0.1] + [2] * num_LHs
                                     )

            for s1 in range(len(all_sampling_dates)):
                for s2 in range(len(all_LHs)):

                    # save the nested list containing the relevant probabilities for this scenario
                    with open("simOutputs_DTK/lik_circulation_numSamples%i_PCirculation%i_samplingDate%i_xLH%i_%s.p" % (ss, round(p_circulation[p_c] * 100), all_sampling_dates[s1],
                                                                                                      round(all_LHs[s2] * 100), filename_suffix), "rb") as f:
                        like_circulation_given_s_n = pickle.load(f)

                    ax_cur = fig.add_subplot(grid[s2+1, s1+1])
                    likelihood_given_observation_plotter(likelihood_given_s_n=like_circulation_given_s_n, diagnostic_names=test_names, ss=ss,
                                                         circulating_string='circulation is',
                                                         detection_limit=0.95, ax=ax_cur)

            # row and column labels
            for s1 in range(len(all_sampling_dates)):
                ax = fig.add_subplot(grid[0, s1+1])
                ax.axis('off')
                ax.set_xlim([0, 1])
                ax.set_ylim([0, 1])
                ax.text(0.5, 0.5, 'sampling day: %i' % all_sampling_dates[s1], horizontalalignment='center',
                        verticalalignment='center', fontweight='bold')
            for s2 in range(len(all_LHs)):
                ax = fig.add_subplot(grid[s2+1, 0])
                ax.axis('off')
                ax.set_xlim([0, 1])
                ax.set_ylim([0, 1])
                ax.text(0, 0.5, 'LH multiplier: %.1f' % all_LHs[s2], horizontalalignment='center',
                        verticalalignment='center', rotation=90, fontweight='bold')

            fig.suptitle('Likelihood of circulation, S = %i, P(circulation)=%.2f' % (ss, p_circulation[p_c]))

            plt.show()

            fig.savefig('figures/Likelihood_circulation_S%i_PCirculation%i_DTK_%s.pdf' % (ss, int(round(p_circulation[p_c]*100)), filename_suffix))

    #
    # Likelihood of circulation given observation - v2
    #
    # same as plot panel above, but now have all P(circulation) and ss values in the same plot (different plots for different LHs
    #    and sampling-day scenarios)

    for s1 in range(len(all_sampling_dates)):
        for s2 in range(len(all_LHs)):
            # Set up the axes with gridspec
            fig = plt.figure(figsize=(14, 10.5))
            # grid = plt.GridSpec(4, 4, hspace=2, wspace=1)

            num_sample_values = len(ss_values)  # will be the columns
            num_p_circulations = len(p_circulation)  # will be the rows
            grid = gridspec.GridSpec((num_p_circulations + 1), (num_sample_values + 1), hspace=1, wspace=0.6,
                                     width_ratios=[0.1] + [2] * num_sample_values,
                                     height_ratios=[0.1] + [2] * num_p_circulations
                                     )

            for ss_index, ss in enumerate(ss_values):
                for p_c in range(len(p_circulation)):

                    # save the nested list containing the relevant probabilities for this scenario
                    with open("simOutputs_DTK/lik_circulation_numSamples%i_PCirculation%i_samplingDate%i_xLH%i_%s.p" % (ss, round(p_circulation[p_c] * 100), all_sampling_dates[s1],
                                                                                                      round(all_LHs[s2] * 100), filename_suffix), "rb") as f:
                        like_circulation_given_s_n = pickle.load(f)

                    ax_cur = fig.add_subplot(grid[p_c+1, ss_index+1])
                    likelihood_given_observation_plotter(likelihood_given_s_n=like_circulation_given_s_n, diagnostic_names=test_names, ss=ss,
                                                         circulating_string='circulation is',
                                                         detection_limit=0.95, line_flag=False, ax=ax_cur)

            # row and column labels
            # row labels
            for ss_index in range(len(ss_values)):
                ax = fig.add_subplot(grid[0, ss_index+1])
                ax.axis('off')
                ax.set_xlim([0, 1])
                ax.set_ylim([0, 1])
                ax.text(0.5, 0.5, 'number samples \n collected: %i' % ss_values[ss_index], horizontalalignment='center',
                        verticalalignment='center', fontweight='bold')
            # column labels
            for p_c in range(len(p_circulation)):
                ax = fig.add_subplot(grid[p_c+1, 0])
                ax.axis('off')
                ax.set_xlim([0, 1])
                ax.set_ylim([0, 1])
                ax.text(0, 0.5, 'Probability of circulation: %.1f' % p_circulation[p_c], horizontalalignment='center',
                        verticalalignment='center', rotation=90, fontweight='bold')

            fig.suptitle('Likelihood of circulation \n Sampling day = %i, LH multiplier = %.2f' % (all_sampling_dates[s1], all_LHs[s2]))

            # plt.show()

            fig.savefig('figures/Likelihood_circulation_samplingDay%i_LH%i_DTK_%s.pdf' % (all_sampling_dates[s1], int(round(all_LHs[s2]*100)), filename_suffix))


#
# Likelihood of no circulation given observation - v1
#

def plot_panel_likelihood_no_circulation(pop_size, test_names, filename_suffix, all_sampling_dates,
                                      all_LHs, ss_values, p_circulation):
    """

    :param pop_size: population sized assumed for plots
    :param test_names: names of diagnostic tests
    :param filename_suffix: indicates which DTK simulation set to use
    :param all_sampling_dates: list of dates when sampling occurs
    :param all_LHs: list of larval habitat multipliers used in simulations
    :param ss_values: list of number of individuals sampled during survey
    :param p_circulation: list of values for the probability of circulation
    """

    for ss_index, ss in enumerate(ss_values):
        for p_c in range(len(p_circulation)):

            # Set up the axes with gridspec
            fig = plt.figure(figsize=(14, 10.5))
            # grid = plt.GridSpec(4, 4, hspace=2, wspace=1)

            num_sampling_dates = len(all_sampling_dates)  # will be the columns
            num_LHs = len(all_LHs)  # will be the rows
            grid = gridspec.GridSpec((num_LHs + 1), (num_sampling_dates + 1), hspace=1, wspace=0.6,
                                     width_ratios=[0.1] + [2] * num_sampling_dates,
                                     height_ratios=[0.1] + [2] * num_LHs
                                     )

            for s1 in range(len(all_sampling_dates)):
                for s2 in range(len(all_LHs)):

                    # save the nested list containing the relevant probabilities for this scenario
                    with open("simOutputs_DTK/lik_no_circulation_numSamples%i_PCirculation%i_samplingDate%i_xLH%i_%s.p" % (ss, round(p_circulation[p_c] * 100), all_sampling_dates[s1],
                                                                                                         round(all_LHs[s2] * 100), filename_suffix), "rb") as f:
                        like_no_circulation_given_s_n = pickle.load(f)

                    ax_cur = fig.add_subplot(grid[s2+1, s1+1])
                    likelihood_given_observation_plotter(likelihood_given_s_n=like_no_circulation_given_s_n, diagnostic_names=test_names, ss=ss,
                                                         circulating_string='circulation is not',
                                                         detection_limit=0.95, ax=ax_cur)

            # row and column labels
            for s1 in range(len(all_sampling_dates)):
                ax = fig.add_subplot(grid[0, s1+1])
                ax.axis('off')
                ax.set_xlim([0, 1])
                ax.set_ylim([0, 1])
                ax.text(0.5, 0.5, 'sampling day: %i' % all_sampling_dates[s1], horizontalalignment='center',
                        verticalalignment='center', fontweight='bold')
            for s2 in range(len(all_LHs)):
                ax = fig.add_subplot(grid[s2+1, 0])
                ax.axis('off')
                ax.set_xlim([0, 1])
                ax.set_ylim([0, 1])
                ax.text(0, 0.5, 'LH multiplier: %.1f' % all_LHs[s2], horizontalalignment='center',
                        verticalalignment='center', rotation=90, fontweight='bold')

            fig.suptitle('Likelihood of no circulation, S = %i, P(circulation)=%.2f' % (ss, p_circulation[p_c]))

            # plt.show()

            fig.savefig('figures/Likelihood_no_circulation_S%i_PCirculation%i_DTK_%s.pdf' % (ss, int(round(p_circulation[p_c]*100)), filename_suffix))


    #
    # Likelihood of no circulation given observation - v2
    #

    # same as plot panel, but now have all P(circulation) and SS values in the same plot (different plots for different LHs
    #    and sampling-day scenarios)

    for s1 in range(len(all_sampling_dates)):
        for s2 in range(len(all_LHs)):
            # Set up the axes with gridspec
            fig = plt.figure(figsize=(14, 10.5))
            # grid = plt.GridSpec(4, 4, hspace=2, wspace=1)

            num_sample_values = len(ss_values)  # will be the columns
            num_p_circulations = len(p_circulation)  # will be the rows
            grid = gridspec.GridSpec((num_p_circulations + 1), (num_sample_values + 1), hspace=1, wspace=0.6,
                                     width_ratios=[0.1] + [2] * num_sample_values,
                                     height_ratios=[0.1] + [2] * num_p_circulations
                                     )

            for ss_index, ss in enumerate(ss_values):
                for p_c in range(len(p_circulation)):

                    # save the nested list containing the relevant probabilities for this scenario
                    with open("simOutputs_DTK/lik_no_circulation_numSamples%i_PCirculation%i_samplingDate%i_xLH%i_%s.p" % (ss, round(p_circulation[p_c] * 100), all_sampling_dates[s1],
                                                                                                      round(all_LHs[s2] * 100), filename_suffix), "rb") as f:
                        like_no_circulation_given_s_n = pickle.load(f)

                    ax_cur = fig.add_subplot(grid[p_c+1, ss_index+1])
                    likelihood_given_observation_plotter(likelihood_given_s_n=like_no_circulation_given_s_n, diagnostic_names=test_names, ss=ss,
                                                         circulating_string='circulation is not',
                                                         detection_limit=0.95, ax=ax_cur)

            # row and column labels
            # row labels
            for ss_index in range(len(ss_values)):
                ax = fig.add_subplot(grid[0, ss_index+1])
                ax.axis('off')
                ax.set_xlim([0, 1])
                ax.set_ylim([0, 1])
                ax.text(0.5, 0.5, 'number samples \n collected: %i' % ss_values[ss_index], horizontalalignment='center',
                        verticalalignment='center', fontweight='bold')
            # column labels
            for p_c in range(len(p_circulation)):
                ax = fig.add_subplot(grid[p_c+1, 0])
                ax.axis('off')
                ax.set_xlim([0, 1])
                ax.set_ylim([0, 1])
                ax.text(0, 0.5, 'Probability of circulation: %.1f' % p_circulation[p_c], horizontalalignment='center',
                        verticalalignment='center', rotation=90, fontweight='bold')

            fig.suptitle('Likelihood of no circulation \n Sampling day = %i, LH multiplier = %.2f' % (all_sampling_dates[s1], all_LHs[s2]))

            # plt.show()

            fig.savefig('figures/Likelihood_no_circulation_samplingDay%i_LH%i_DTK_%s.pdf' % (all_sampling_dates[s1], int(round(all_LHs[s2]*100)), filename_suffix))


















