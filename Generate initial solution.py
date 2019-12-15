# Generate initial solution
# Use advanced mathematical functions.
import math
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
        if seeding_rules[key][1] < 10:
            rules[key][1] = seeding_rules[key][1] * 10
            rules[key][2] = 2.31
        else:
            rules[key][1] = 75
            rules[key][2] = (math.log(seeding_rules[key][1] / 100))/(-1)
# Set the seeded accessibility parameters.
acc_dd = [1, 1, 1, 1, 20, 20, 20, 1]
acc_w = [0, 0, 0, 0, 0.5, 0.5, 0.5, 0]

# Specify the number of parameters that need to be generated (3 per rule)
# per iteration.
par_no = luc * act * 3 + 16
# Initialise a list specifying the parameter
par_index = [0] * par_no

# Initialise a temp dictionary to store the iteration rules.
iterations_rules = {}
for i in range(0, luc):
    for j in range(0, act):
        key = "from " + luc_names[i] + " to " + luc_names[j + pas]
        iterations_rules[key] = [0, 0, 0]
# Initialise an index to track the parameter value.
index = 0
# Generate sample values for the inertia point and tail values.
for i in range(0, act):
    key = "from " + luc_names[i + pas] + " to " + luc_names[i + pas]
    # Generate the values
    for c in range(0, 3):
        # If c is zero, generate a point value.
        if c == 0:
            # Store the point in the tracking array.
            par_index[index] = rules[key][0]
            # Add 1 to index.
            index = index + 1
        if c == 1:
            # Store the point in the tracking array.
            par_index[index] = rules[key][1]
            # Add 1 to index.
            index = index + 1
        if c == 2:
            # Store the point in the tracking array.
            par_index[index] = rules[key][2]
            # Add 1 to index.
            index = index + 1
# Generate values for the conversion point and tail values.
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
                    # Store the point in the tracking array.
                    par_index[index] = rules[key][0]
                    # Add 1 to index.
                    index = index + 1
                # If c is one, generate an 'a' value for the tail.
                if c == 1:
                    # Store the point in the tracking array.
                    par_index[index] = rules[key][1]
                    # Add 1 to index.
                    index = index + 1
                # If c is one, generate an 'b' value for the tail.
                if c == 2:
                    # Store the point in the tracking array.
                    par_index[index] = rules[key][2]
                    # Add 1 to index.
                    index = index + 1
# Initialise a temp list to store the iteration distance
# decays & weights.
iteration_dd = [0] * 8
iteration_w = [0] * 8
# Generate sample values for the accessibility distance decay
# parameters.
for i in range(0, act):
    # Store the point in the tracking array.
    par_index[index] = acc_dd[i]
    # Add 1 to index.
    index = index + 1
# Generate sample values for the accessibility weight
# parameters.
for i in range(0, act):
    # Store the point in the tracking array.
    par_index[index] = acc_w[i]
    # Add 1 to index.
    index = index + 1
# Input the rules generated into the model.
for i in range(0, luc):
    for j in range(0, act):
        key = "from " + luc_names[i] + " to " + luc_names[j + pas]
        fu_elem = j
        lu_elem = i
        y0 = rules[key][0]
        a = rules[key][1]
        b = rules[key][2]
        set_exp_rule(project_file, fu_elem, lu_elem, y0, a, b, min_in)
# Input the accessibility parameters generated into the model.
for i in range(0, act):
    fu_elem = i
    dd = acc_dd[i]
    weight = acc_w[i]
    set_acc(project_file, fu_elem, dd, weight)
# Set the random seed
set_rand(project_file, rseed)
# Run the model to generate the simulated output.
run_metro(project_file, log_file, working_directory, geo_cmd)
# Read in the simulated map.
smap = read_map(smap_path)
# Calculate Fuzzy Kappa.
FK = fuzzy_kappa(amap_path, smap_path, fuzzy_coefficients)
# Store the calculated Fuzzy Kappa Simulation.
FKS = fks(omap_path, amap_path, smap_path,
          fuzzy_trans_coefficients)
# Store the calculated area-weighted clumpiness error.
AWCE = area_weighted_clu_error(amap, smap, mask,
                               luc, pas, act,
                               luc_count)

# Generate a list of parameter names.
par_names = ["0"] * par_no
for i in range(0, par_no):
    par_names[i] = "Par_" + str(i)

# Specify the file to save the output
tested_parameters_file = output_path + "original_parameters.csv"
# Now write to the file.
with open(tested_parameters_file, "w", newline='') as csv_file:
    writer = csv.writer(csv_file)
    # Write a header line as the parameter names.
    writer.writerow(par_names)
    # Write the parameter values.
    writer.writerow(par_index)

# Specify the file to save the metrics
output_metrics_file = output_path + "original_metrics.csv"
store = [0] * 3
# Now write to the file.
with open(output_metrics_file, "w", newline='') as csv_file:
    writer = csv.writer(csv_file)
    # Write a header line for the metrics.
    values = ["FKS", "FK", "AWCE"]
    writer.writerow(values)
    # Write the output metrics.
    store[0] = FKS
    store[1] = FK
    store[2] = AWCE
    writer.writerow(store)

# Beeping noise to mark completion
import winsound
Freq = 2500 # Set Frequency To 2500 Hertz
Dur = 1000 # Set Duration To 1000 ms == 1 second
winsound.Beep(Freq,Dur)

# Completed!
