"""
# Set the accessiblity distance decay and weighting
parameter values.
"""

# Multi-dimensional array handling.
import numpy as np
# Editing of xml tree structure.
import xml.etree.ElementTree as ET

# Function.


def set_acc(project_path, fu_elem, dd, weight):
    # Open tree for editing.
    source = open(project_path)
    tree = ET.parse(source)  # Parse georpoject file as xml type data structure for modification
    root = tree.getroot()
    # Access distance-decay parameter.
    spline_dd = root[2][1][3][0][0][8][0][1][0][0][0][0][fu_elem]
    spline_dd.text = str(dd)
    spline_w = root[2][1][3][0][0][8][0][2][0][0][0][0][fu_elem]
    spline_w.text = str(weight)
    tree.write(project_path)

# Test.
"""
project_path = ("C:\\Geonamica\\Metronamica\\Madrid\\"
                "Madrid.geoproj")
fu_elem = 0
dd = 17.5
w = 0.33
# Run function.
set_acc(project_path, fu_elem, dd, w)
fu_elem = 5
dd = 4
w = 0.80
# Run function 2.0.
set_acc(project_path, fu_elem, dd, w)
"""