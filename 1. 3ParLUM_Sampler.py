"""
This is the base program used to sample Land-use model parameters for the
generation of a population of solutions for subsequent optimisation
using multi-objective optimisation.
"""

# Allows for multi-dimensional array handling.
import numpy as np
# Allows for maps to be stored as 2d arrays.
from read_map import read_map
# Determines the size (# of cells) of square neighbourhoods.
from considered_distances import considered_distances
# Used to generate a dictionary of rules for seeding.
from seed_rules2 import seed_rules2
# Set neighbourhood rules based on a four-point structure.
from set_NR import set_exp_rule
# Set the accessibility parameters.
from set_acc import set_acc
# Set the random number seed in the Metronamica file.
from set_rand import set_rand
# Run the Metronamica model to generate output.
from run_metro import run_metro
# Module for the calculation of Fuzzy Kappa
from fuzzy_kappa import fuzzy_kappa
# Module for calculation of Fuzzy Kappa Simulation.
from fuzzy_kappa import fks
# Module for calculation of Absolute Area Weighted Avg. Clumpiness Error (AAWCE)
from area_weighted_clu import area_weighted_clu_error
# Interaction with csv file format.
import csv
# Track the time.
import time

# Specify the base path to the directory containing this program.
base_path = "C:\\Users\\charl\\OneDrive\\Documents\\LUM_sampler\\"
# Set the case study
case_study = "Madrid"
# Set the path to the data.
data_path = base_path + "EU_data\\"
# Set the path to the directory containing the input parameters.
input_path = base_path + "EU_output\\" + case_study + "\\"
# Set the path to the directory where the output will be store.
output_path = base_path + "Sampler_output\\"

# Specify the original map path at time slice 1.
omap_path = data_path + case_study + "\\" + case_study.lower() + "_1990.asc"
# Specify the actual map path at time slice 2.
amap_path = data_path + case_study + "\\" + case_study.lower() + "_2000.asc"
# Specify the masking map.
mask_path = data_path + case_study + "\\" + case_study.lower() + "_mask.asc"

# Specify the fuzzy weights for the calculation of Fuzzy Kappa and FKS
fuzzy_coefficients = data_path + "coeff13.txt"
fuzzy_trans_coefficients = data_path + "coefficients13.txt"

# Set the working directory the contains the geoproject file.
working_directory = ("C:\\Geonamica\\Metronamica\\" + case_study)
# Set the project file path.
project_file = working_directory + "\\" + case_study + ".geoproj"
# Set the path to the command line version of Geonamica
geo_cmd = "C:\\Program Files (x86)\\Geonamica\\Metronamica\\GeonamicaCmd.exe"
# Set the path to the log file.
log_file = base_path + "LogSettings.xml"
# Set the path to the simulated output map
smap_path = (
    working_directory + "\\Log\\Land_use\\"
                        "Land use map_2000-Jan-01 00_00_00.rst"
)

# Set the land-use class names.
luc_names = ["Natural areas", "Arable land", "Permanent crops", "Pastures",
             "Agricultural areas", "Residential", "Industry & commerce",
             "Recreation areas", "Forest", "Road & rail", "Port area",
             "Airports", "Mine & dump sites", "Fresh water", "Marine water"]
# Set the land-use class parameters: number of land-use classes, passive,
# feature, and active.
luc = len(luc_names)
pas = 1
fea = 6
act = luc - (pas + fea)
# Specify the maximum neighbourhood size distance considered
max_distance = 5
# Specify the random number seed.
rseed = 1000
# Specify the minimum influence. Influence values below this are set to 0.
min_in = 10**-4

# Read in the map for time slice 1.
omap = read_map(omap_path)
# Read in the map for time slice 2.
amap = read_map(amap_path)
# Read in the masking map.
mask = read_map(mask_path)
# Analyse the input maps for evaluation purposes
map_dimensions = np.shape(omap)
rows = map_dimensions[0]
cols = map_dimensions[1]

# Count the presence of each land-use class in the actual map. This is
# used in the calculation of area-weighted average clumpiness across the
# active classes.
luc_count = [0] * luc
for i in range(0, rows):
    for j in range(0, cols):
        if mask[i, j] > 0:
            luc_count[amap[i, j]] = luc_count[amap[i, j]] + 1

# Determine the distances that will be analysed using the module considered
# distances.
temp = considered_distances(max_distance)
# Store the list of considered distances as a variable.
cd = temp[0]
# Store the total number of distances considered
cdl = temp[1]
# Determine the maximum neighbourhood size (unit) from considered distances
N_all = [1, 8, 12, 16, 32, 28, 40, 40, 20]
N = []
for c in range(0, max_distance):
    N.append(N_all[c])

# Specify the files with important inertia/conversion points and tails.
points_rule_file = input_path + "Rules\\conversion_pts.txt"
tails_rule_file = input_path + "Rules\\att_rules.txt"
# Store the information about points and tails into matrices.
att_rules = np.loadtxt(tails_rule_file)
convo_pts = np.loadtxt(points_rule_file)
# Analyse conversion points. Remove those for feature classes (cannot effect
# model).
for i in range(0, luc):
    for j in range(0, act):
        if i > (act + pas - 1):
            convo_pts[i, j] = 0

# Specify how the rules are going to be input.
meta_parameter_seeding = False
# Generate rules accordingly
if meta_parameter_seeding is True:
    theta_st = 0.040
    theta_cp = 0.040
    theta_it = 0.008
    seeding_rules = seed_rules2(omap, amap, mask, max_distance, luc_names,
                                luc, act, pas, att_rules, theta_st, theta_cp,
                                theta_it, project_file)
else:
    # Specify the input rules file
    seeding_rules_file = input_path + "Rules\\final_rules.csv"
    # Read in rules, store in a dictionary.
    seeding_rules = {}
    for i in range(0, luc):
        for j in range(0, act):
            key = "from " + luc_names[i] + " to " + luc_names[j + pas]
            seeding_rules[key] = [0, 0, 0, 5]
    # Read inputs from csv file
    with open(seeding_rules_file, 'r', newline='') as f:
        readCSV = csv.reader(f)
        next(f)  # This skips the header line
        for row in readCSV:
            i = row[0]
            j = row[1]
            key = "from " + i + " to " + j
            seeding_rules[key][0] = float(row[2])
            seeding_rules[key][1] = float(row[3])
            seeding_rules[key][2] = float(row[4])
            seeding_rules[key][3] = float(row[5])
# Convert the rules to the form ae^-bx. Each rule is parameterised by 3 parameters:
# C/I, the inertia or conversion point at distance 0;
# a, the initial strength of the function; &
# b, the rate of decay.

# Generate a rules dictionary to store values of this form.
rules = {}
# Populate the dictionary
for i in range(0, luc):
    for j in range(0, act):
        key = "from " + luc_names[i] + " to " + luc_names[j + pas]
        rules[key] = [0, 0, 0]
# Convert the rules generated via seeding to the required format.
for i in range(0, luc):
    for j in range(0, act):
        key = "from " + luc_names[i] + " to " + luc_names[j + pas]
        rules[key][0] = seeding_rules[key][0]
        rules[key][1] = seeding_rules[key][1] * 10
        rules[key][2] = 2.31
# Set the seeded accessibility parameters.
acc_dd = [1, 1, 1, 1, 20, 20, 20, 1]
acc_w = [0, 0, 0, 0, 0.5, 0.5, 0.5, 0]

# Initialise the iterative testing.
# Specify the number of iterations.
max_iterations = 3
# Specify the number of parameters that need to be generated (3 per rule)
# per iteration.
par_no = luc * act * 3 + 16
# Initialise a list specifying the parameter
par_index = [0] * par_no
# Specify the parameter indices for the inertia rules.
for i in range(0, act):
    for c in range(0, 3):
        index = i*3 + c
        if c == 0:
            par_index[index] = "IP: " + luc_names[i + pas]
        if c == 1:
            par_index[index] = "IT_a: " + luc_names[i + pas]
        if c == 2:
            par_index[index] = "IT_b: " + luc_names[i + pas]
# Specify the parameter indices for the conversion rules.
counter = act * 3
for i in range(0, act):
    for j in range(0, luc):
        # Skip inertia rules.
        if i + pas == j :
            pass
        else:
            for c in range(0, 3):
                index = index + 1
                if c == 0:
                    par_index[index] = ("CP, from " + luc_names[j] + " to " +
                                        luc_names[i + pas])
                if c == 1:
                    par_index[index] = ("CT_a, from " + luc_names[j] + " to " +
                                        luc_names[i + pas])
                if c == 2:
                    par_index[index] = ("CT_b, from " + luc_names[j] + " to " +
                                        luc_names[i + pas])
# Specify an array to store the parameter values generated.
pv_store = np.zeros(shape=(max_iterations, par_no))
# Specify an array to store the objective metrics calculated.
objective_store = np.zeros(shape=(max_iterations, 3))
# Specify the distribution and required parameters.
use_normal_distribution = False
use_triangle_distribution = True
# Specify required parameters for distribution.
if use_normal_distribution is True:
    # Specify a standard deviation for each type of parameter as a percentage
    # of the mean value.
    # Inertia point standard deviation percentage.
    inertia_point_sd_rate = 0.0
    # Inertia tail a and b standard deviation percentages.
    inertia_tail_a_sd_rate = 0.0
    inertia_tail_b_sd_rate = 0.0
    # Conversion point standard deviation percentage.
    conversion_point_sd_rate = 0.0
    # Inertia tail a and b standard deviation percentages.
    conversion_tail_a_sd_rate = 0.0
    conversion_tail_b_sd_rate = 0.0
if use_triangle_distribution is True:
    # Specify a width for each type of parameter as a percentage
    # of the mean value.
    inertia_point_width_rate = 0.25
    # Inertia tail a and b standard deviation percentages.
    inertia_tail_a_width_rate = 0.25
    inertia_tail_b_width_rate = 0.25
    # Conversion point standard deviation percentage.
    conversion_point_width_rate = 0.25
    # Inertia tail a and b standard deviation percentages.
    conversion_tail_a_width_rate = 0.25
    conversion_tail_b_width_rate = 0.25
    # Accessibility distance decay and weight sd %.
    acc_dd_width_rate = 0.25
    acc_w_width_rate = 0.25

# Start a timer to track the duration
start = time.time()

# Start the iterative generation of parameter values.
counter = 0
while counter < max_iterations:
    # Initialise a temp dictionary to store the iteration rules.
    iterations_rules = {}
    for i in range(0, luc):
        for j in range(0, act):
            key = "from " + luc_names[i] + " to " + luc_names[j + pas]
            iterations_rules[key] = [0, 0, 0, 5]
    # Generate the rules by sampling from a normal distribution.
    if use_normal_distribution is True:
        # Initialise an index to track the parameter value.
        index = 0
        # Generate sample values for the inertia point and tail values.
        for i in range(0, act):
            key = "from " + luc_names[i + pas] + " to " + luc_names[i + pas]
            # Generate the values
            for c in range(0, 3):
                # If c is zero, generate a point value.
                if c == 0:
                    # Check if seeded value is 0. If not, perform analysis.
                    if rules[key][0] != 0:
                        mu = rules[key][0]
                        # Calculate the standard deviation
                        sigma = mu * inertia_point_sd_rate
                        # Generate the parameter value.
                        sample_pv = np.random.normal(mu, sigma)
                        # Load the parameter value into the dictionary.
                        iterations_rules[key][0] = sample_pv
                        # Store the point in the tracking array.
                        pv_store[counter, index] = sample_pv
                        # Add 1 to index.
                        index = index + 1
                    # If not, add 1 to index.
                    else:
                        index = index + 1
                # If c is one, generate an 'a' value for the tail.
                if c == 1:
                    # Check if seeded value is 0. If not, perform analysis.
                    if rules[key][1] != 0:
                        mu = rules[key][1]
                        # Calculate the standard deviation
                        sigma = mu * inertia_tail_a_sd_rate
                        # Generate the parameter value.
                        sample_pv = np.random.normal(mu, sigma)
                        # Load the parameter value into the dictionary.
                        iterations_rules[key][1] = sample_pv
                        # Store the point in the tracking array.
                        pv_store[counter, index] = sample_pv
                        # Add 1 to index.
                        index = index + 1
                    # If not, add 1 to index.
                    else:
                        index = index + 1
                # If c is two, generate an 'b' value for the tail.
                if c == 2:
                    mu = rules[key][2]
                    # Calculate the standard deviation
                    sigma = mu * inertia_tail_b_sd_rate
                    # Generate the parameter value.
                    sample_pv = np.random.normal(mu, sigma)
                    # Load the parameter value into the dictionary.
                    iterations_rules[key][2] = sample_pv
                    # Store the point in the tracking array.
                    pv_store[counter, index] = sample_pv
                    # Add 1 to index.
                    index = index + 1
        # Generate sample values for the conversion point and tail values.
        for i in range(0, act):
            for j in range(0, luc):
                # Skip in an inertia interaction.
                if i + pas == j:
                    pass
                else:
                    # Specify the key
                    key = "from " + luc_names[j] + " to " + luc_names[i + pas]
                    for c in range(0, 3):
                        # If c is zero, generate a point value.
                        if c == 0:
                            # Check if seeded value is 0. If not, perform analysis.
                            if rules[key][0] != 0:
                                mu = rules[key][0]
                                # Calculate the standard deviation
                                sigma = mu * conversion_point_sd_rate
                                # Generate the parameter value.
                                sample_pv = np.random.normal(mu, sigma)
                                # Load the parameter value into the dictionary.
                                iterations_rules[key][0] = sample_pv
                                # Store the point in the tracking array.
                                pv_store[counter, index] = sample_pv
                                # Add 1 to index.
                                index = index + 1
                            # If not, add 1 to index.
                            else:
                                index = index + 1
                        # If c is one, generate an 'a' value for the tail.
                        if c == 1:
                            # Check if seeded value is 0. If not, perform analysis.
                            if rules[key][1] != 0:
                                mu = rules[key][1]
                                # Calculate the standard deviation
                                sigma = mu * conversion_tail_a_sd_rate
                                # Generate the parameter value.
                                sample_pv = np.random.normal(mu, sigma)
                                # Load the parameter value into the dictionary.
                                iterations_rules[key][1] = sample_pv
                                # Store the point in the tracking array.
                                pv_store[counter, index] = sample_pv
                                # Add 1 to index.
                                index = index + 1
                            # If not, add 1 to index.
                            else:
                                index = index + 1
                        # If c is two, generate an 'b' value for the tail.
                        if c == 2:
                            mu = rules[key][2]
                            # Calculate the standard deviation
                            sigma = mu * conversion_tail_b_sd_rate
                            # Generate the parameter value.
                            sample_pv = np.random.normal(mu, sigma)
                            # Load the parameter value into the dictionary.
                            iterations_rules[key][2] = sample_pv
                            # Store the point in the tracking array.
                            pv_store[counter, index] = sample_pv
                            # Add 1 to index.
                            index = index + 1
    # Generate the rules by sampling from a normal distribution.
    if use_triangle_distribution is True:
        # Initialise an index to track the parameter value.
        index = 0
        # Generate sample values for the inertia point and tail values.
        for i in range(0, act):
            key = "from " + luc_names[i + pas] + " to " + luc_names[i + pas]
            # Generate the values
            for c in range(0, 3):
                # If c is zero, generate a point value.
                if c == 0:
                    # Check if seeded value is 0. If not, perform analysis.
                    if rules[key][0] != 0:
                        # Specify the mode.
                        mode = rules[key][0]
                        # Calculate the width.
                        width = mode * inertia_point_width_rate
                        # Calculate the left & right.
                        left = mode - width
                        right = mode + width
                        # Generate the point.
                        sample_pv = np.random.triangular(left, mode, right)
                        # Load the point into the dictionary.
                        iterations_rules[key][0] = sample_pv
                        # Store the point in the tracking array.
                        pv_store[counter, index] = sample_pv
                        # Add 1 to index.
                        index = index + 1
                    # If not, add 1 to index.
                    else:
                        index = index + 1
                # If c is one, generate an 'a' value for the tail.
                if c == 1:
                    # Check if seeded value is 0. If not, perform analysis.
                    if rules[key][1] != 0:
                        # Specify the mode.
                        mode = rules[key][1]
                        # Calculate the width.
                        width = mode * inertia_tail_a_width_rate
                        # Calculate the left & right.
                        left = mode - width
                        right = mode + width
                        # Generate the point.
                        sample_pv = np.random.triangular(left, mode, right)
                        # Load the point into the dictionary.
                        iterations_rules[key][1] = sample_pv
                        # Store the point in the tracking array.
                        pv_store[counter, index] = sample_pv
                        # Add 1 to index.
                        index = index + 1
                    # If not, add 1 to index.
                    else:
                        index = index + 1
                # If c is two, generate an 'b' value for the tail.
                if c == 2:
                    # Specify the mode.
                    mode = rules[key][2]
                    # Calculate the width.
                    width = mode * inertia_tail_b_width_rate
                    # Calculate the left & right.
                    left = mode - width
                    right = mode + width
                    # Generate the point.
                    sample_pv = np.random.triangular(left, mode, right)
                    # Load the point into the dictionary.
                    iterations_rules[key][2] = sample_pv
                    # Store the point in the tracking array.
                    pv_store[counter, index] = sample_pv
                    # Add 1 to index.
                    index = index + 1
        # Generate sample values for the conversion point and tail values.
        for i in range(0, act):
            for j in range(0, luc):
                # Skip in an inertia interaction.
                if i + pas == j:
                    pass
                else:
                    # Specify the key
                    key = "from " + luc_names[j] + " to " + luc_names[i + pas]
                    for c in range(0, 3):
                        # If c is zero, generate a point value.
                        if c == 0:
                            # Check if seeded value is 0. If not, perform analysis.
                            if rules[key][0] != 0:
                                # Specify the mode.
                                mode = rules[key][0]
                                # Calculate the width.
                                width = mode * conversion_point_width_rate
                                # Calculate the left & right.
                                left = mode - width
                                right = mode + width
                                # Generate the point.
                                sample_pv = np.random.triangular(left, mode, right)
                                # Load the point into the dictionary.
                                iterations_rules[key][0] = sample_pv
                                # Store the point in the tracking array.
                                pv_store[counter, index] = sample_pv
                                # Add 1 to index.
                                index = index + 1
                            # If not, add 1 to index.
                            else:
                                index = index + 1
                        print(index)
                        # If c is one, generate an 'a' value for the tail.
                        if c == 1:
                            # Check if seeded value is 0. If not, perform analysis.
                            if rules[key][1] != 0:
                                # Specify the mode.
                                mode = rules[key][1]
                                # Calculate the width.
                                width = mode * conversion_tail_a_width_rate
                                # Calculate the left & right.
                                left = mode - width
                                right = mode + width
                                # Generate the point.
                                sample_pv = np.random.triangular(left, mode, right)
                                # Load the point into the dictionary.
                                iterations_rules[key][1] = sample_pv
                                # Store the point in the tracking array.
                                pv_store[counter, index] = sample_pv
                                # Add 1 to index.
                                index = index + 1
                            # If not, add 1 to index.
                            else:
                                index = index + 1
                        print(index)
                        # If c is one, generate an 'b' value for the tail.
                        if c == 2:
                            # Specify the mode.
                            mode = rules[key][2]
                            # Calculate the width.
                            width = mode * conversion_tail_b_width_rate
                            # Calculate the left & right.
                            left = mode - width
                            right = mode + width
                            # Generate the point.
                            sample_pv = np.random.triangular(left, mode, right)
                            # Load the point into the dictionary.
                            iterations_rules[key][2] = sample_pv
                            # Store the point in the tracking array.
                            pv_store[counter, index] = sample_pv
                            # Add 1 to index.
                            index = index + 1
            # Initialise a temp list to store the iteration distance
            # decays & weights.
            iteration_dd = acc_dd
            iteration_w = acc_w
            # Generate sample values for the accessibility distance decay
            # parameters.
            for i in range(0, act):
                # Generate the distance decay values.
                mode = acc_dd[i]
                # Calculate the width.
                width = mode * acc_dd_width_rate
                # Calculate the left & right.
                left = mode - width
                right = mode + width
                # Generate the point.
                sample_pv = np.random.triangular(left, mode, right)
                # Load the point into the list.
                iteration_dd[i] = sample_pv
                # Store the point in the tracking array.
                pv_store[counter, index] = sample_pv
                # Add 1 to index.
                index = index + 1
            # Generate sample values for the accessibility weight
            # parameters.
            for i in range(0, act):
                # Generate the distance decay values.
                mode = acc_w[i]
                # If zero skip..
                if mode == 0:
                    # Load the point into the list.
                    iteration_w[i] = sample_pv
                    # Store the point in the tracking array.
                    pv_store[counter, index] = sample_pv
                    # Add 1 to index.
                    index = index + 1
                else:
                    # Calculate the width.
                    width = mode * acc_w_width_rate
                    # Calculate the left & right.
                    left = mode - width
                    right = mode + width
                    # Generate the point.
                    sample_pv = np.random.triangular(left, mode, right)
                    # Load the point into the list.
                    iteration_w[i] = sample_pv
                    # Store the point in the tracking array.
                    pv_store[counter, index] = sample_pv
                    # Add 1 to index.
                    index = index + 1
    # Input the rules generated into the model.
    for i in range(0, luc):
        for j in range(0, act):
            key = "from " + luc_names[i] + " to " + luc_names[j + pas]
            fu_elem = j
            lu_elem = i
            y0 = iterations_rules[key][0]
            a = iterations_rules[key][1]
            b = iterations_rules[key][2]
            set_exp_rule(project_file, fu_elem, lu_elem, y0, a, b, min_in)
    # Input the accessibility parameters generated into the model.
    for i in range(0, act):
        fu_elem = i
        dd = iteration_dd[i]
        weight = iteration_w[i]
        set_acc(project_file, fu_elem, dd, weight)
    # Set the random seed
    set_rand(project_file, rseed)
    # Run the model to generate the simulated output.
    run_metro(project_file, log_file, working_directory, geo_cmd)
    # Read in the simulated map.
    smap = read_map(smap_path)
    # Store the calculated Fuzzy Kappa.
    objective_store[counter, 0] = fuzzy_kappa(amap_path, smap_path,
                                              fuzzy_coefficients)
    # Store the calculated Fuzzy Kappa Simulation.
    objective_store[counter, 1] = fks(omap_path, amap_path, smap_path,
                                      fuzzy_trans_coefficients)
    # Store the calculated area-weighted clumpiness error.
    objective_store[counter, 2] = area_weighted_clu_error(amap, smap, mask,
                                                          luc, pas, act,
                                                          luc_count)
    # Add 1 to prevent an infinite loop!
    counter = counter + 1
    # Provide user feedback.
    print("Iterations completed: " + str(counter))

# Track the end time of the calibration.
end = time.time()
# Determine the duration of the calibration method in hours.
duration = end - start
duration_h = duration/3600
# Record the duration of calibration and number of iterations performed.
output_duration_file = output_path + "\\duration.txt"
f = open(output_duration_file, "w")
f.write("duration: " + str(duration_h))
f.close()

# Specify the file to save the output
tested_parameters_file = output_path + "tested_parameters.csv"
# Now write to the file.
with open(tested_parameters_file, "w", newline='') as csv_file:
    writer = csv.writer(csv_file)
    # Write a header line as the parameter names.
    writer.writerow(par_index)
    # Write the parameter values.
    for i in range(0, max_iterations):
        store = pv_store[i, :]
        writer.writerow(store)

# Specify the file to save the metrics
output_metrics_file = output_path + "output_metrics.csv"
store = [0] * 4
# Now write to the file.
with open(output_metrics_file, "w", newline='') as csv_file:
    writer = csv.writer(csv_file)
    # Write a header line for the metrics.
    values = ["Set", "FK", "FKS", "AWCE"]
    writer.writerow(values)
    # Write the output metrics.
    for i in range(0, max_iterations):
        store[0] = i
        store[1] = objective_store[i, 0]
        store[2] = objective_store[i, 1]
        store[3] = objective_store[i, 2]
        writer.writerow(store)

# Beeping noise to mark completion
import winsound
Freq = 2500 # Set Frequency To 2500 Hertz
Dur = 1000 # Set Duration To 1000 ms == 1 second
winsound.Beep(Freq,Dur)

# Completed!
