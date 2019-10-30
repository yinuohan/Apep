'''
Gives a quick overview of the darks and flats needed. 
'''

import os
import csv
import numpy as np

values_dark = []
values_flat = []

with open('Path//Info.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    
    for row in csv_reader:
        row = np.array(row)
        value = "_".join(row[[7,8,15,16,17,18,14]])
        
        if value not in values_dark and 'DARK' in row[1]:
            values_dark.append(value)
            
        if value not in values_flat and 'FLAT' in row[1]:
            values_flat.append(value)

for value in sorted(values_dark):
    print(value)

print()

for value in sorted(values_flat):
    print(value)