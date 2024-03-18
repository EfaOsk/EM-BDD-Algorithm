import matplotlib.pyplot as plt
from itertools import groupby
import statistics
from scipy.stats import pearsonr
import numpy as np


def parse_line_to_dict(line):
    # Splitting the line by commas
    parts = line.split(',')

    # Extracting values from each part
    model = parts[0].split(':')[1].strip()
    obs = int(parts[1].split(':')[1].strip())
    run = int(parts[2].split(':')[1].strip())
    l_org = float(parts[3].split(':')[1].strip())
    num_iter = int(parts[4].strip())
    l_final = float(parts[5].strip())
    runtime = float(parts[6].strip())

    # Constructing the result dictionary
    result = {
        'name': model,
        'Obs': obs,
        'Lorg': l_org,
        'num iter': num_iter,
        'Lfinal': l_final,
        'runtime': runtime
    }

    return result

def read_file_to_dict(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            if '|' not in line:
                data_dict = parse_line_to_dict(line)
                data.append(data_dict)
    return data




# Function to create a plot like the one in the image
def plot_bdd_size_vs_observation(data_x, data_y):
    """
    Plot BDD Size vs Observation Sequences length.

    Parameters:
    - data_x: Array-like, the x-axis data (Observation Sequences length)
    - data_y: Array-like, the y-axis data (Number of nodes in A)
    """
    # Scatter plot of the data
    plt.scatter(data_x, data_y, color='blue', label='Data')

    # Calculate trendline
    z = np.polyfit(data_x, data_y, 1)
    p = np.poly1d(z)
    plt.plot(data_x, p(data_x), linestyle='dashed', color='red', label='Trendline')

    # Labeling the plot
    plt.title('BDD Size vs Number of States')
    plt.xlabel('Number of states (N)')
    plt.ylabel(r'Number of nodes in $\Delta$')
    plt.legend()

    # Show plot
    plt.show()

data = []
with open("new_g_n.txt", 'r') as file:
# with open("test.txt", 'r') as file:
    for line in file:
        data.append(int(line.strip()))

# Example data (this would be your actual data)
data_x = np.arange(5, 33)
# data_x = np.arange(5, 50)  # This represents the observation sequence length
data_y = data  # This represents the BDD size

# Call the function with the example data
plot_bdd_size_vs_observation(data_x, data_y)


correlation_coefficient, _ = pearsonr(data_x, data_y)
print(correlation_coefficient)