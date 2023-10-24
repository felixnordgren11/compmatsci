import csv
import numpy as np
import matplotlib.pyplot as plt

def load_csv(filename: str):
    with open(filename) as csvFile:
        reader = csv.reader(csvFile)
        rows = []
        for row in reader:
            rows.append(row)
    return rows

#data1 = load_csv('')
#edata1 = {v[1].lower(): [float(k) for k in v[3:]] for v in data}

filenames = ['fcc100a108.txt']
data108 = np.loadtxt(filenames[0], delimiter=',', skiprows=1, dtype=str)
newdata = []
for i in range(len(data108)):
    newdata.append(data108[i].split(' '))
    
def pos(i):
    x = newdata[i][0]
    y = newdata[i][1]
    z = newdata[i][2]

x1 = newdata[0][0] #First index is the atom number, second index is the column (0 = x, 1 = y, 2 = z)

print(x1)
