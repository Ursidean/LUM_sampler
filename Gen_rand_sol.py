"""
This is the base program used to sample Land-use model parameters for the
generation of a population of solutions for subsequent optimisation
using multi-objective optimisation.
"""

# Perform advanced mathematical operations
import math
# Allows for multi-dimensional array handling.
import numpy as np
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
omap_path = ("C:\\Geonamica\\Metronamica\\Madrid\\Data\\"
             + case_study.lower() + "_1990.asc")
# Specify the actual map path at time slice 2.
amap_path = ("C:\\Geonamica\\Metronamica\\Madrid\\Data\\"
             + case_study.lower() + "_2000.asc")
# Specify the masking map.
mask_path = ("C:\\Geonamica\\Metronamica\\Madrid\\Data\\"
             + case_study.lower() + "_mask.asc")

# Specify the fuzzy weights for the calculation of Fuzzy Kappa and FKS
fuzzy_coefficients = base_path + "coeff12.txt"
fuzzy_trans_coefficients = base_path + "coefficients12.txt"

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
             "Recreation areas", "Forest", "Road & rail",
             "Airports", "Mine & dump sites", "Fresh water"]
# Set the land-use class parameters: number of land-use classes, passive,
# feature, and active.
luc = len(luc_names)
pas = 1
fea = 4
act = luc - (pas + fea)
# Specify the maximum neighbourhood size distance considered
max_distance = 8
# Specify the lower and upper bounds for generating the solutions.
lb_file_path = ("C:\\Users\\charl\\Dropbox\PhD\\Journal Papers"
                "\\3. Integrated Automated Calibration"
                "\\1. Optimisation inputs\\"
                "Bounds_lower.txt")
ub_file_path = ("C:\\Users\\charl\\Dropbox\PhD\\Journal Papers"
                "\\3. Integrated Automated Calibration"
                "\\1. Optimisation inputs\\"
                "Bounds_upper.txt")

# Read in the lower bound values.
lower_bounds = np.loadtxt(lb_file_path, delimiter="\t")
# Read in the upper bound values.
upper_bounds = np.loadtxt(ub_file_path, delimiter="\t")

# Generate random solutions
no_solns = 298
no_pars = len(upper_bounds)
store = np.zeros(shape=(no_solns, no_pars))
# Begin iterative generation.
for i in range(0, no_solns):
    for j in range(0, no_pars):
        lb = lower_bounds[j]
        ub = upper_bounds[j]
        if lb != ub:
            store[i, j] = np.random.uniform(lb, ub)
        else:
            store[i, j] = lb

header_list = ["0"] * no_pars
dummy = [0] * no_pars

# Write output to a .csv file.
output_file = ("C:\\Users\\charl\\OneDrive\\Documents\\"
               "LUM_sampler\\Sampler_output\\Parameters_random1.csv")
with open(output_file, "w", newline='') as csv_file:
    writer = csv.writer(csv_file)
    # Write a header line for the parameters.
    for i in range(0, no_pars):
        header_list[i] = "Par " + str(i + 1)
    writer.writerow(header_list)
    for i in range(0, no_solns):
        for j in range(0, no_pars):
            dummy[j] = store[i, j]
        writer.writerow(dummy)

# Input parameter and check output model.
