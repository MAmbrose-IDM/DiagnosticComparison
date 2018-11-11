import numpy as np
import csv
import math
import matplotlib.pyplot as plt
from operator import mul
from functools import reduce

from simulateExposureAndTestResults import simulate_over_transmission_intensities, time_since_infection
from calculateNumSamplesNeeded import calculate_num_samples_each_sim
from scenarioPlotter import samples_needed_each_test_plotter, diagnostic_probability_plotter, \
    time_since_infection_plotter, case_seasonality_plotter

# initialize parameters that describe the system
pop_size = 500
num_case_in_year_list = [15, 30] + list(range(45, round(pop_size), round(pop_size/9))) + [round(pop_size * 1.25)]
# seasonal_scalar = [1] * 365
seasonal_scalar = [(1 + 0.5 * np.sin(2 * np.pi / 365 * y - 166)) for y in range(1, 366)]
surveillance_days = [250]  # [90] [250]
initialize_years = 2
test_detection_probs = [([0.0]*1 + [0.01*np.exp(np.log(70)/3)**y for y in range(4)] + [0.7]*10
                         + [0.7*np.exp(np.log(0.00007)/30)**y for y in range(31)] + [0]*250),
                        ([0.0]*1 + [0.01*np.exp(np.log(80)/2)**y for y in range(3)] + [0.8]*15
                         + [0.8*np.exp(np.log(0.00008)/150)**y for y in range(151)] + [0]*150),
                        ([0.0]*2 + [0.01*np.exp(np.log(90)/10)**y for y in range(11)] + [0.9]*20
                         + [0.9*np.exp(np.log(0.00009)/300)**y for y in range(301)])]

# test_detection_probs = [([0.0]*2 + list(np.linspace(0.001, 0.7, 5)) + list(np.linspace(0.7, 0.5, 20))
#                          + list(np.linspace(0.5, 0, 20))),
#                         ([0.0]*2 + list(np.linspace(0.001, 0.9, 5)) + list(np.linspace(0.9, 0.7, 50))
#                          + list(np.linspace(0.7, 0, 30))),
#                         ([0.0]*2 + list(np.linspace(0.001, 0.9, 14)) + list(np.linspace(0.9, 0.5, 100))
#                          + list(np.linspace(0.5, 0, 120)))]

diagnostic_names = ['RDT (pretend)', 'uRDT (pretend)', 'serology (pretend)']
false_pos_prob = [0.0]*len(test_detection_probs)
num_sims_each = 200
confidence_level = 0.95


# plot diagnostic tests
diagnostic_probability_plotter(test_detection_probs=test_detection_probs, diagnostic_names=diagnostic_names)

# plot expected number of cases on each day of year (seasonality of transmission)
case_seasonality_plotter(pop_size=pop_size, seasonal_scalar=seasonal_scalar,
                         num_case_in_year_list=num_case_in_year_list, surveillance_days=surveillance_days)


# Start by generating the datasets giving how long individuals had been infected when sampling occurred and whether
# they tested positive for each of the tests
sim_output = simulate_over_transmission_intensities(pop_size=pop_size,
                                                    num_case_in_year_list=num_case_in_year_list,
                                                    initialize_years=initialize_years,
                                                    seasonal_scalar=seasonal_scalar,
                                                    surveillance_days=surveillance_days,
                                                    test_detection_probs=test_detection_probs,
                                                    false_pos_prob=false_pos_prob,
                                                    num_sims_each=num_sims_each)

# if instead of simulating, we want to import simulation results directly from csv, do so here:
# read in csv file on simulation output
# with open('testData_MRA.csv') as data_file:
#     data_reader = csv.reader(data_file, delimiter=',')
#     line_count = 0
#     for row in data_reader:
#         if line_count == 0:
#             print(f'Column names are {", ".join(row)}')
#             line_count = line_count + 1
#         else:
#             print(f'\t{row[0]} had test result {row[1]} for month1, test1')
#             line_count = line_count + 1

# plot simulation results
time_since_infection_plotter(pop_size=pop_size, num_case_in_year_list=num_case_in_year_list,
                             sim_time_since_last_infection=sim_output[0])

# Calculate the number of samples needed for each of the fraction-of-tests-positive scenarios
sim_num_samples_needed = calculate_num_samples_each_sim(sim_number_positive=sim_output[1],
                                                        pop_size=pop_size,
                                                        confidence_level=confidence_level)
# In sim_num_samples_needed,
#    The outer-most level has one entry for each of the diagnostic tests
#    The middle level has one entry for each scenario on the number of cases in a year (from num_case_in_year_list)
#    The inner-most level has one entry for each simulation run.


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



# print(values_for_plot[0][0])
# print(values_for_plot[1][0])
# print(values_for_plot[2][0])
#
#
# print(sim_num_samples_needed[0][0])
# print(sim_num_samples_needed[0][1])
#
# print(np.sum(sim_output[0][0][0]))
# print(np.sum(sim_output[0][0][1]))
#
# print(np.sum(sim_output[1][0][0]))
# print(np.sum(sim_output[1][0][1]))


# Call the plotting function to create the mean lines and 95CI areas for each of the diagnostic tests
samples_needed_each_test_plotter(x_values=num_case_in_year_list, values_for_plot=values_for_plot, pop_size=pop_size,
                                 diagnostic_names=diagnostic_names, confidence_level=confidence_level)


