'''Adjust neighbourhood rules'''
import math
'''
Used to set neighbourhood rules in Metronamica project file. Neighbourhood rules
are assumed to be structured as follows:
|
o
|\
| \
|  \
|   \
|    o
|     \
|      \
|       \ 
|        o
|         \
-----------o---------------------------

Where the circles are:
The influence at 0 (inertia or conversion)
The influence at 1 
The influence at 2
The location where influence is 0
'''

#Modules
import numpy as np
import xml.etree.ElementTree as ET

#Function
def set_lp_rule(project_path,fu_elem,lu_elem,y0,y1,y2,xe):
    #y0 is the influence at a distance of 0
    #y1 is the influence at a distance of 1
    #y2 is the influence at a distance of 2
    #xe is the distance for an influence of 0
    source=open(project_path)
    tree=ET.parse(source)                                                 #Parse georpoject file as xml type data structure for modification
    root=tree.getroot()
    #Access specified neighbourhood rule
    spline=root[2][1][3][0][0][1][0][1][0][0][fu_elem][0][lu_elem][0]
    #Store array of values
    arrange=np.array([[0,y0],[1,y1],[2,y2],[xe,0]])
    #Find length of arrange
    size4=len(arrange)
    #Remove current neighbourhood rule values
    for point in spline.findall('point'):
        spline.remove(point)
    #Dictionary to store inputs for neighbourhood rules
    inputs={}                                                             
    for i in range(0,size4):
        key='new_value'+str(i)
        inputs[key]={'y':str(arrange[i,1]),'x':str(arrange[i,0])}
    #Insert points, commit to case study    
    for i in range(0,size4):                                              
        key='new_value'+str(i)
        child=ET.SubElement(spline, 'point', attrib=inputs[key]) 
    
    tree.write(project_path)
    source.close()


def set_exp_rule(project_path, fu_elem, lu_elem, y0, a, b, min_in):
    # y0 is the influence at a distance of 0
    # a is the maximum influence.
    # b is the rate of decay.
    source = open(project_path)
    tree = ET.parse(source)  # Parse georpoject file as xml type data structure for modification
    root = tree.getroot()
    # Calculate the exponential decay values.
    y1 = a * math.exp(-b * 1);
    y2 = a * math.exp(-b * 2);
    y3 = a * math.exp(-b * 3);
    y4 = a * math.exp(-b * 4);
    y5 = a * math.exp(-b * 5);
    y6 = a * math.exp(-b * 6);
    y7 = a * math.exp(-b * 7);
    # Evaluate the values. If below the minimum influence set to 0.
    if abs(y1) < min_in:
        y1 = 0
    if abs(y2) < min_in:
        y2 = 0
    if abs(y3) < min_in:
        y3 = 0
    if abs(y4) < min_in:
        y4 = 0
    if abs(y5) < min_in:
        y5 = 0
    if abs(y6) < min_in:
        y6 = 0
    if abs(y7) < min_in:
        y7 = 0
    # Access specified neighbourhood rule
    spline = root[2][1][3][0][0][1][0][1][0][0][fu_elem][0][lu_elem][0]
    # Store array of values
    arrange = np.array([[0, y0], [1, y1], [2, y2], [3, y3], [4, y4], [5, y5],
                        [6, y6], [7, y7], [8, 0]])
    # Find length of arrange
    size4 = len(arrange)
    # Remove current neighbourhood rule values
    for point in spline.findall('point'):
        spline.remove(point)
    # Dictionary to store inputs for neighbourhood rules
    inputs = {}
    for i in range(0, size4):
        key = 'new_value' + str(i)
        inputs[key] = {'y': str(arrange[i, 1]), 'x': str(arrange[i, 0])}
    # Insert points, commit to case study
    for i in range(0, size4):
        key = 'new_value' + str(i)
        child = ET.SubElement(spline, 'point', attrib=inputs[key])

    tree.write(project_path)
    source.close()