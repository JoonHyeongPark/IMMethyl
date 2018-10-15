import pandas as pd
import matplotlib.pyplot as plt

gap = [0.05, 0.025, 0.01, 0.001, 0.0001]
                
def Draw(cutoff_gap) :

    input_file = open(str(cutoff_gap) + "cutoff.count_array.txt", 'r')
    input_table = input_file.readlines()

    histogram_vector = []

    for line in input_table : histogram_vector.append(int(line.split()[2]))

    plt.hist(histogram_vector)
    plt.show()

    input_file.close()

    return

Draw(0.05)
