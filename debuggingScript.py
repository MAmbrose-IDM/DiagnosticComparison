# insert default values if none are provided:
    if line_color_list is None:
        line_color_list = ['p-', 'r-', 'o-']  # colors to use in plot for each detection probability
    if number_pos_samples is None:
        number_pos_samples = [1, 5, 10, 50]  # fraction of samples that test positive
    if numbers_sampled_list is None:
        numbers_sampled_list = list(range(1, 101))
    if ax is None:
        ax = plt.gca()

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
                               line_color_list[jj], label=number_pos_samples[jj]
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
                     line_color_list[jj], label=number_pos_samples[jj]
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