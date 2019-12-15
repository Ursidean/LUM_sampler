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
from set_NR import set_lp_rule
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
meta_parameter_seeding = True
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

# Begin the iterative sampling of rules, create an array to store for each
# iteration.
max_iterations = 500
# The number of parameters is a function of the model land-use classes and
# parameterisation method used.
par_no = luc * act * 2
# Three storing arrays are used to track parameter indices, parameter values
# and objective metrics.
parameter_indices_store = [0] * par_no
parameter_values_store = np.zeros(shape=(max_iterations, par_no))
objective_store = np.zeros(shape=(max_iterations, 3))
# Populate the parameter indices store
for i in range(0, luc):
    for j in range(0, act):
        for c in range(0, 2):
            if c == 0:
                par = "point"
            else:
                par = "tail"
            parameter_indices_store[i * luc + j * 2 + c] = (
                "from " + luc_names[i] + " to " + luc_names[j + pas] +
                " " + par
            )
# Specify the distribution and required parameters.
use_normal_distribution = False
use_triangle_distribution = True
# Specify required parameters for distribution.
if use_normal_distribution is True:
    # For the normal distribution specify a standard deviation for each type of
    # parameter. These values are in absolute terms of influence.
    inertia_point_sd = 100
    conversion_point_sd = 10
    inertia_tail_sd = 2
    conversion_tail_sd = 2
if use_triangle_distribution is True:
    # For the triangular distribution specify the width for each type of
    # parameter. These values are in absolute terms of influence.
    inertia_point_width = 200
    conversion_point_width = 20
    inertia_tail_width = 4
    conversion_tail_width = 4

# Start the iterative generation of solutions.
counter = 0
while counter < max_iterations:
    # Initialise a temp dictionary to store the iteration rules.
    iterations_rules = {}
    for i in range(0, luc):
        for j in range(0, act):
            key = "from " + luc_names[i] + " to " + luc_names[j + pas]
            iterations_rules[key] = [0, 0, 0, 5]
    # Sample using the normal distribution.
    if use_normal_distribution is True:
        # Generate sample values.
        for i in range(0, luc):
            for j in range(0, act):
                # Specify the key
                key = "from " + luc_names[i] + " to " + luc_names[j + pas]
                # Switch between points and tails. Both are parameterized by
                # one value.
                for c in range(0, 2):
                    # If c is 0 a point is being analysed.
                    if c == 0:
                        # If a conversion point is included, sample as specified.
                        if convo_pts[i, j] == 1:
                            # Specify the mean.
                            mu = seeding_rules[key][0]
                            # Specify the standard deviation.
                            if i == j + pas:
                                sigma = inertia_point_sd
                            else:
                                sigma = conversion_point_sd
                            # Generate the point.
                            sample_pv = np.random.normal(mu, sigma)
                            # Load the point into the dictionary.
                            iterations_rules[key][0] = sample_pv
                            # Store the point in the tracking array.
                            parameter_values_store[
                                counter, (i * luc + j * 2 + c)
                            ] = sample_pv
                    # If c is 1 a tail is being analysed.
                    if c == 1:
                        # If an attraction rule is included, sample as specified.
                        if att_rules[i, j] == 1:
                            # Specify the mean.
                            mu = seeding_rules[key][1]
                            # Specify the standard deviation.
                            if i == j + pas:
                                sigma = inertia_tail_sd
                            else:
                                sigma = conversion_tail_sd
                            # Generate the tail point at distance 1.
                            sample_pv = np.random.normal(mu, sigma)
                            # Generate the tail point at distance 2.
                            sample_pv_2 = sample_pv * 0.1
                            # Load the point into the dictionary.
                            iterations_rules[key][1] = sample_pv
                            iterations_rules[key][2] = sample_pv_2
                            # Store the point in the tracking array.
                            parameter_values_store[
                                counter, (i * luc + j * 2 + c)
                            ] = sample_pv
    # Sample using the triangle distribution.
    if use_triangle_distribution is True:
        # Generate sample values.
        for i in range(0, luc):
            for j in range(0, act):
                # Specify the key
                key = "from " + luc_names[i] + " to " + luc_names[j + pas]
                # Switch between points and tails. Both are parameterized by
                # one value.
                for c in range(0, 2):
                    # If c is 0 a point is being analysed.
                    if c == 0:
                        # If a conversion point is included, sample as specified.
                        if convo_pts[i, j] == 1:
                            # Specify the mode.
                            mode = seeding_rules[key][0]
                            # Specify the bounds.
                            if i == j + pas:
                                left = (seeding_rules[key][0] -
                                        inertia_point_width)
                                right = (seeding_rules[key][0] +
                                         inertia_point_width)
                            else:
                                left = (seeding_rules[key][0] -
                                        conversion_point_width)
                                right = (seeding_rules[key][0] +
                                         conversion_point_width)
                            # Generate the point.
                            sample_pv = np.random.triangular(left, mode,
                                                             right)
                            # Load the point into the dictionary.
                            iterations_rules[key][0] = sample_pv
                            # Store the point in the tracking array.
                            parameter_values_store[
                                counter, (i * luc + j * 2 + c)
                            ] = sample_pv
                    # If c is 1 a tail is being analysed.
                    if c == 1:
                        # If an attraction rule is included, sample as specified.
                        if att_rules[i, j] == 1:
                            # Specify the mode.
                            mode = seeding_rules[key][1]
                            # Specify the bounds.
                            if i == j + pas:
                                left = (seeding_rules[key][1] -
                                        inertia_tail_width)
                                right = (seeding_rules[key][1] +
                                         inertia_tail_width)
                            else:
                                left = (seeding_rules[key][1] -
                                        conversion_tail_width)
                                right = (seeding_rules[key][1] +
                                         conversion_tail_width)
                            # Generate the tail point at distance 1.
                            sample_pv = np.random.triangular(left, mode,
                                                             right)
                            # Generate the tail point at distance 2.
                            sample_pv_2 = sample_pv * 0.1
                            # Load the point into the dictionary.
                            iterations_rules[key][1] = sample_pv
                            iterations_rules[key][2] = sample_pv_2
                            # Store the point in the tracking array.
                            parameter_values_store[
                                counter, (i * luc + j * 2 + c)
                            ] = sample_pv
    # Input the rules generated into the model.
    for i in range(0, luc):
        for j in range(0, act):
            key = "from " + luc_names[i] + " to " + luc_names[j + pas]
            fu_elem = j
            lu_elem = i
            y0 = iterations_rules[key][0]
            y1 = iterations_rules[key][1]
            y2 = iterations_rules[key][2]
            xe = max_distance
            set_lp_rule(project_file, fu_elem, lu_elem, y0, y1, y2, xe)
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

# Specify the file to save the tested parameter values.
tested_parameters_file = output_path + "tested_parameter.csv"
# Now write to the file.
with open(tested_parameters_file, "w", newline='') as csv_file:
    writer = csv.writer(csv_file)
    # Write a header line as the parameter names.
    writer.writerow(parameter_indices_store)
    # Write the parameter values.
    for i in range(0, max_iterations):
        store = parameter_values_store[i, :]
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

# Sounds to mark completion
import winsound
Freq = 2500 # Set Frequency To 2500 Hertz
Dur = 1000 # Set Duration To 1000 ms == 1 second
winsound.Beep(Freq,Dur)

# Completed!
