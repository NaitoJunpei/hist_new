import numpy as np
import pandas as pd

def make_csv(filename) :
    file1 = open(filename)
    file2 = open(filename + ".csv", "w")
    file2.write("time\n")
    for line in file1 :
        arr = line.split()
        newstr = "\n".join(arr)
        file2.write(newstr + "\n")

    file1.close()
    file2.close()
