import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm


def prob_dist_n_p_plotter(prob_num_pos, diagnostic_names, pop_size, ax=None, xmax=None):
    """
    :param prob_num_pos: nested list containing the probability that a population where disease is
        circulating contains n_p individuals who would test positive on each diagnostic test (given circulation or
        given no circulation).
        prob_num_pos[test][n_p]
    :param diagnostic_names: list of character strings naming the diagnostic tests plotted
    :param pop_size: assumed number of individuals in the population
    :param ax: axes to use for the current plot
    :param xmax: xlim max to plot if specified
    """
    if ax is None:
        ax = plt.gca()

    colormap0 = cm.tab10.colors
    color_order = [5, 9, 8, 7, 6, 2, 1, 0, 4, 8, 3]
    colormap = [colormap0[y] for y in color_order]


    for test in range(len(prob_num_pos)):
        # create histogram instead of line plot
        data_0 = [[y] * int(np.round(prob_num_pos[test][y] * 400)) for y in
                  range(len(prob_num_pos[test]))]
        data = [val for sublist in data_0 for val in sublist]
        # if (np.max(data) - np.min(data)) > 0:
        #     ax.hist(data, bins=(np.max(data) - np.min(data)), density=True, alpha=0.6, color=colormap[test])
        # else:
        #     ax.hist(data, bins=pop_size, density=True, alpha=0.6, color=colormap[test])

        ax.plot(prob_num_pos[test], '.', color=colormap[test], alpha=0.9, markersize=3,
                label=diagnostic_names[test])  #  linewidth=1.0,

    # draw vertical line showing average population size across all simulations
    ax.axvline(x=pop_size, linestyle=':', color='k', alpha=0.5, linewidth=1)

    if xmax is None:
        ax.set_xlim(0, pop_size)
    else:
        ax.set_xlim(0, xmax)

    ax.set_ylim(0, 1)
    ax.set_ylabel('probability density')
    ax.set_xlabel('number of individuals who \n would test positive if sampled')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.set_title('Population size = %i; Confidence level = %.2f' % (pop_size, round(confidence_level, 2)), y=1.08)
    # plt.legend(title='Diagnostic method')


# prob_dist_n_p_plotter(prob_num_pos=prob_num_pos_circulation,
#                                    diagnostic_names=['RDT', 'hsRDT', 'Serology'], pop_size=271, ax=None)


def prob_positive_sample_plotter(prob_pos_sample, diagnostic_names, pop_size,
                                 circulating_string='circulation is not',
                                 detection_limit=0.95, ax=None, line_flag=False, xmax=None):
    """
    :param prob_pos_sample: nested list containing the probability that a at least one sampled
        individual will test positive for each diagnostic test given that circulation is occurring (or is not occurring)
        and ss individuals are sampled (to calculate this number, we integrated over different values of n_p)
        prob_pos_sample[test][ss]
    :param diagnostic_names: list of character strings naming the diagnostic tests plotted
    :param pop_size: assumed number of individuals in the population
    :param circulating_string: what should the y-axis label indicate: whether circulation is or is not occurring
    :param detection_limit: show where we want to draw a vertical line for minimum number of samples so that the
        probability we get at least one positive sample given circulation is at least detection_limit
    :param ax: axes to use for the current plot
    :param line_flag: should horizontal and vertical lines be drawn to show the number of samples at the detection limit
    """
    if ax is None:
        ax = plt.gca()


    colormap0 = cm.tab10.colors
    color_order = [5, 9, 8, 7, 6, 2, 1, 0, 4, 8, 3]
    colormap = [colormap0[y] for y in color_order]

    # how far out into samples should we look? not much further than detection_limit, probably
    first_index_max = 0

    for test in range(len(prob_pos_sample)):
        ax.plot(prob_pos_sample[test], 'k', color=colormap[test], alpha=0.7, linewidth=3,
                label=diagnostic_names[test])
        # draw vertical line where we exceed the detection_limit threshold
        if line_flag:
            first_index = next(x[0] for x in enumerate(prob_pos_sample[test]) if x[1] >= detection_limit)
            first_index_max = np.max([first_index_max, first_index])
            ax.axvline(x=first_index, linestyle=':', color=colormap[test], alpha=0.7, linewidth=1)

    if line_flag:
        # draw horizontal line showing detection limit
        ax.axhline(y=detection_limit, linestyle=':', color='k', alpha=0.7, linewidth=1)

    if xmax is None:
        ax.set_xlim(0, np.min([pop_size, (first_index_max+10)]))
    else:
        ax.set_xlim(0, xmax)

    # ax.set_ylabel('probability at least one \n sample will test positive \n if %s occurring' % circulating_string)
    ax.set_ylabel('probability at least one \n sample will test positive')
    ax.set_xlabel('number of individuals tested')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.set_title('Population size = %i; Confidence level = %.2f' % (pop_size, round(confidence_level, 2)), y=1.08)
    # plt.legend(title='Diagnostic method')


# # load the nested list containing the relevant probabilities for this scenario
# with open("simOutputs_DTK/prob_pos_sample_given_circulation_samplingDate%i_xLH%i.p" % (all_sampling_dates[0],
#                                                                                        round(all_LHs[0] * 100)),
#           "rb") as f:
#     prob_pos_sample_given_circulation = pickle.load(f)
# prob_positive_sample_given_circulation_plotter(prob_pos_sample_given_circulation=prob_pos_sample_given_circulation,
#                                        diagnostic_names=['RDT', 'hsRDT', 'Serology'], pop_size=271,
#                                        detection_limit=0.95, ax=None)


def likelihood_given_observation_plotter(likelihood_given_s_n, diagnostic_names, ss,
                                         circulating_string='circulation is not',
                                         detection_limit=0.95, ax=None):
    """
    :param likelihood_given_s_n: nested list containing the likelihood there is (or is not) circulation given the
    observed number of negative samples (s_n) for a given diagnostic test
        like_no_circulation_given_s_n[test][s_n]
    :param diagnostic_names: list of character strings naming the diagnostic tests plotted
    :param ss: number of samples collected
    :param circulating_string: what should the y-axis label indicate: whether circulation is or is not occurring
    :param detection_limit: show where we want to draw a vertical line for minimum number of samples that would have
        to be negative before we thought there was a likelihood of at least detection_limit that there was no circulation
    :param ax: axes to use for the current plot
    """
    if ax is None:
        ax = plt.gca()

    colormap0 = cm.tab10.colors
    color_order = [5, 9, 8, 7, 6, 2, 1, 0, 4, 8, 3]
    colormap = [colormap0[y] for y in color_order]

    # how far out into samples should we look? not much further than detection_limit, probably
    first_index_max = 0

    for test in range(len(likelihood_given_s_n)):
        ax.plot(likelihood_given_s_n[test], 'k', color=colormap[test], alpha=0.7, linewidth=3,
                label=diagnostic_names[test])

        # draw vertical line where we exceed the detection_limit threshold
        first_index = next(x[0] for x in enumerate(likelihood_given_s_n[test]) if x[1] >= detection_limit)
        first_index_max = np.max([first_index_max, first_index])
        ax.axvline(x=first_index, linestyle=':', color=colormap[test], alpha=0.7, linewidth=1)

    # draw horizontal line showing detection limit
    ax.axhline(y=detection_limit, linestyle=':', color='k', alpha=0.7, linewidth=1)

    ax.set_ylim(0, 1)
    # ax.set_xlim(0, np.min([ss, (first_index_max+10)]))
    ax.set_ylabel('likelihood %s \n occurring' % circulating_string)
    ax.set_xlabel('number of negative tests \n (out of %i)' % ss)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.set_title('Population size = %i; Confidence level = %.2f' % (pop_size, round(confidence_level, 2)), y=1.08)
    # plt.legend(title='Diagnostic method')


# # load the nested list containing the relevant probabilities for this scenario
# p_circulation = [0.1, 0.3, 0.5]
# pop_size = 271
# all_sampling_dates = [45, 190, 300, 360]
# all_LHs = [0.6, 0.7, 0.8, 1.0]
# ss_values = [int(round(pop_size * y)) for y in [0.1, 0.3, 0.5]]
# with open("simOutputs_DTK/lik_no_circulation_numSamples%i_PCirculation%i_samplingDate%i_xLH%i.p" % (ss_values[0], round(p_circulation[0] * 100), all_sampling_dates[0],
#                                                                                                     round(all_LHs[0] * 100)), "rb") as f:
#     like_no_circulation_given_s_n = pickle.load(f)
#
# with open("simOutputs_DTK/lik_circulation_numSamples%i_PCirculation%i_samplingDate%i_xLH%i.p" % (ss_values[0], round(p_circulation[0] * 100), all_sampling_dates[0],
#                                                                                                     round(all_LHs[0] * 100)), "rb") as f:
#     like_circulation_given_s_n = pickle.load(f)
#
# like_no_circulation_given_observation_plotter(like_no_circulation_given_s_n=like_no_circulation_given_s_n,
#                                               diagnostic_names=['RDT', 'hsRDT', 'Serology'], ss=ss_values[0],
#                                               detection_limit=0.95, ax=None)















