import numpy as np
import matplotlib.pyplot as plt
import pickle

from scenarioPlotter import samples_needed_each_test_plotter, diagnostic_probability_plotter, \
    time_since_infection_plotter, case_seasonality_plotter, num_cases_legend_plotter, diagnostic_legend_plotter
from likelihoodCirculation import calculate_likelihood_no_circulation_each_sim

# initialize parameters that describe the system
pop_size = 500
num_case_in_year_list = [round(pop_size * y) for y in [0.03, 0.06, 0.09, 0.22, 0.35, 0.48, 0.61, 0.74, 0.87, 1, 1.25]]
initialize_years = 2
test_detection_probs = [([0.0]*1 + [0.01*np.exp(np.log(70)/3)**y for y in range(4)] + [0.7]*10
                         + [0.7*np.exp(np.log(0.00007)/30)**y for y in range(31)] + [0]*250),
                        ([0.0]*1 + [0.01*np.exp(np.log(80)/2)**y for y in range(3)] + [0.8]*15
                         + [0.8*np.exp(np.log(0.00008)/150)**y for y in range(151)] + [0]*150),
                        ([0.0]*2 + [0.01*np.exp(np.log(90)/10)**y for y in range(11)] + [0.9]*20
                         + [0.9*np.exp(np.log(0.00009)/300)**y for y in range(301)])]


diagnostic_names = ['RDT (pretend)', 'uRDT (pretend)', 'serology (pretend)']
# TODO: replace with true values; redo full analysis
false_pos_prob = [0.001, 0.002, 0.001]
num_sims_each = 200
confidence_level = 0.95

# Lists of the three seasonality/sampling scenarios to explore
seasonal_scalar_list = [[1] * 365, [(1 + 0.5 * np.sin(2 * np.pi / 365 * y - 166)) for y in range(1, 366)],
                        [(1 + 0.5 * np.sin(2 * np.pi / 365 * y - 166)) for y in range(1, 366)]]
surveillance_days_list = [[250], [250], [90]]

# number of samples collected (list)
number_tested = 100
number_tested_list = [25, 50, 100, 200, 400]

# probability of circulation
prob_circulation = 0.2
prob_circulation_list = [0.01, 0.1, 0.2, 0.4, 0.6]


# Load output from previous run (if this parameter set has not already been saved, run it
#   from run_diagnosticTestComparisons_plotPanel.py

sim_output_list = pickle.load(open("sim_output_list_pop%i_CI%i.p" % (pop_size, round(confidence_level * 100)), "rb"))

# calculate the number of samples positive
# Calculate likelihood of circulation if all samples negative
calculate_likelihood_no_circulation_each_sim(sim_output_list=sim_output_list, pop_size=pop_size,
                                          number_tested=number_tested, false_pos_prob=false_pos_prob,
                                          prob_circulation=prob_circulation)






# TODO: haven't started converting the next section into likelihood of circulation figure panel - actually, I don't think I want to plot num_samples_needed
# At this point, I think I'm going to migrate over to working on DTK output primarily, so I'm not sure this section will be touched for a while (if ever...)







# Create multi-panel plot
# - Each columns of the new plot panel should correspond to a different assumption about seasonality of cases given
# # circulation is occurring.
#   The first row shows the expected number of cases through the year. The second row shows the proportion of the
#   population that has had their most recent infection within a certain number of days. The final row shows the
#   number of samples needed to achieve a particular certainty that disease would be detected if circulation were
#   occurring.
# - The last column (before the data) shows the legend for the different transmission intensity scenarios and the
#   plot of the diagnostic test detection probabilities

# Set up the axes with gridspec
fig = plt.figure(figsize=(14, 10.5))
grid = plt.GridSpec(6, 4, hspace=1, wspace=0.6)

# iterate through three seasonality/sampling scenarios
for scenario in range(len(surveillance_days_list)):
    # values for this seasonality / surveillance scenario
    seasonal_scalar = seasonal_scalar_list[scenario]
    surveillance_days = surveillance_days_list[scenario]
    sim_output = sim_output_list[scenario]
    sim_num_samples_needed = sim_num_samples_needed_list[scenario]

    # Plotting

    # first plotted row: seasonality in cases
    ax_seasonality = fig.add_subplot(grid[0:2, scenario])
    # plot expected number of cases on each day of year (seasonality of transmission)
    case_seasonality_plotter(pop_size=pop_size, seasonal_scalar=seasonal_scalar,
                             num_case_in_year_list=num_case_in_year_list, surveillance_days=surveillance_days,
                             ax=ax_seasonality)

    # second plotted row: time since last infection
    ax_time_since_infection = fig.add_subplot(grid[2:4, scenario])
    # plot simulation results - time since infection for all individuals in all simulations
    transmission_intensity_legend = time_since_infection_plotter(pop_size=pop_size,
                                                                 num_case_in_year_list=num_case_in_year_list,
                                                                 sim_time_since_last_infection=sim_output[0],
                                                                 ax=ax_time_since_infection)

    # third plotted row: number of samples needed
    ax_samples_needed = fig.add_subplot(grid[4:6, scenario])

    # Get the mean values as well as the 95CI for plotting lines and shaded CI regions
    # save results in lines_for_plot
    values_for_plot = [[[None]*len(num_case_in_year_list) for n in range(3)]
                       for m in range(len(test_detection_probs))]
    for i1 in range(len(sim_num_samples_needed)):
        for i2 in range(len(sim_num_samples_needed[i1])):
            # mean value
            values_for_plot[i1][0][i2] = np.percentile(sim_num_samples_needed[i1][i2], 50)
            # upper 95% CI
            values_for_plot[i1][1][i2] = min(pop_size, np.percentile(sim_num_samples_needed[i1][i2], 97.5))
            # lower 95% CI
            values_for_plot[i1][2][i2] = np.percentile(sim_num_samples_needed[i1][i2], 2.5)

    # Call the plotting function to create the mean lines and 95CI areas for each of the diagnostic tests
    samples_needed_each_test_plotter(x_values=num_case_in_year_list, values_for_plot=values_for_plot, pop_size=pop_size,
                                     diagnostic_names=diagnostic_names, confidence_level=confidence_level,
                                     num_case_in_year_list=num_case_in_year_list, ax=ax_samples_needed)


# Add legends for different transmission intensities and diagnostic tests.
ax_legend_intensity = fig.add_subplot(grid[0:3, -1])
num_cases_legend_plotter(num_case_in_year_list, ax_legend_intensity)

ax_legend_test = fig.add_subplot(grid[3, -1])
diagnostic_legend_plotter(diagnostic_names, ax_legend_test)


# plot diagnostic tests
ax_diagnostic_test = fig.add_subplot(grid[4:6, -1])
diagnostic_probability_plotter(test_detection_probs=test_detection_probs, diagnostic_names=diagnostic_names)

fig.suptitle('Population size = %i; Confidence level = %.2f' % (pop_size, round(confidence_level, 2)))

plt.show()

fig.savefig('LikelihoodCirculation_pop%i_CI%i.pdf' % (pop_size, round(confidence_level * 100)))
