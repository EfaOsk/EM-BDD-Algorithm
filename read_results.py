import matplotlib.pyplot as plt

def read_results_file(filename):
    """
    Read the results from a file and return a list of dictionaries containing the data.
    """
    results = []
    with open(filename, 'r') as file:
        for line in file:
            data = line.strip().split(' | ')
            result = {}
            for item in data:
                key, value = item.split(': ')
                result[key.strip()] = value.strip()
            results.append(result)
    return results

def print_results(results):
    """
    Print the results in a formatted table.
    """
    print("{:<5} {:<5} {:<10} {:<20} {:<20} {:<20} {:<20}".format("N", "M", "T", "Encode", "DdManager vars", "DdManager nodes", "DdManager memory"))
    for result in results:
        print("{:<5} {:<5} {:<5} {:<10} {:<20} {:<20} {:<20} {:<20}".format(result["N"], result["M"], result["T"], result["Encode"], result["DdManager vars"], result["DdManager nodes"], result["DdManager reorderings"], result["DdManager memory"]))


def plot_nodes_vs_m(results, results1):
    """
    Plot the number of DdManager nodes vs M for N=2 and T=3.
    """
    encode_true_data = []
    encode_false_data = []

    for result in results:
        if result["M"] == "2" and result["T"] == "2" and int(result["N"]) <= 10:
            if result["Encode"] == "TRUE":
                encode_true_data.append((int(result["N"]), int(result["DdManager nodes"])))
            elif result["Encode"] == "FALSE":
                encode_false_data.append((int(result["N"]), int(result["DdManager nodes"])))

    encode_true_data.sort(key=lambda x: x[0])
    encode_false_data.sort(key=lambda x: x[0])

    n_true, nodes_true = zip(*encode_true_data)
    n_false, nodes_false = zip(*encode_false_data)

    # plt.plot(n_true, nodes_true, color='blue', label='Order encoding')
    # plt.plot(n_false, nodes_false, color='red', label='Direct encoding')

    print(nodes_true)
    print(nodes_false)

    encode_true_data = []
    encode_false_data = []

    for result in results1:
        if result["M"] == "2" and result["T"] == "2" and int(result["N"]) <= 10:
            if result["Encode"] == "TRUE":
                encode_true_data.append((int(result["N"]), int(result["DdManager nodes"])))
            elif result["Encode"] == "FALSE":
                encode_false_data.append((int(result["N"]), int(result["DdManager nodes"])))

    encode_true_data.sort(key=lambda x: x[0])
    encode_false_data.sort(key=lambda x: x[0])

    plt.plot(n_true, nodes_true, color='green', label='Order encoding + Reordering')
    plt.plot(n_false, nodes_false, color='orange', label='Direct encoding + Reordering')


    n_true, nodes_true = zip(*encode_true_data)
    n_false, nodes_false = zip(*encode_false_data)

    plt.xlabel('N')
    plt.ylabel('DdManager Nodes')
    plt.title('Number of Nodes vs N for M=2 and T=2')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_nodes_vs_n(results):#, results1):
    """
    Plot the number of DdManager nodes vs N for M=2 and T=2.
    """
    encode_true_data = []
    encode_false_data = []

    for result in results:
        if result["M"] == "5" and result["T"] == "2":
            # if result["Encode"] == "TRUE":
            #     encode_true_data.append((int(result["N"]), int(result["DdManager nodes"])))
            # elif result["Encode"] == "FALSE":
                encode_false_data.append((int(result["N"]), int(result["DdManager nodes"])))

    # encode_true_data.sort(key=lambda x: x[0])
    encode_false_data.sort(key=lambda x: x[0])

    # n_true, nodes_true = zip(*encode_true_data)
    n_false, nodes_false = zip(*encode_false_data)

    # plt.plot(n_true, nodes_true, color='blue', label='Order encoding')
    plt.plot(n_false, nodes_false) #, color='red', label='Direct encoding')

    # encode_true_data = []
    # encode_false_data = []

    # for result in results1:
    #     if result["M"] == "2" and result["T"] == "2" and int(result["N"]) < 6:
    #         if result["Encode"] == "TRUE":
    #             encode_true_data.append((int(result["N"]), int(result["DdManager nodes"])))
    #         elif result["Encode"] == "FALSE":
    #             encode_false_data.append((int(result["N"]), int(result["DdManager nodes"])))

    # encode_true_data.sort(key=lambda x: x[0])
    # encode_false_data.sort(key=lambda x: x[0])

    # n_true, nodes_true = zip(*encode_true_data)
    # n_false, nodes_false = zip(*encode_false_data)

    # plt.plot(n_true, nodes_true, color='green', label='Order encoding + Reordering')
    # plt.plot(n_false, nodes_false, color='orange', label='Direct encoding + Reordering')

    plt.xlabel('Number of states')
    plt.ylabel('Number of Nodes')
    plt.title('Number of Nodes vs N for M=5 and T=2')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    filename_without_reordering = "results_new_not_RO.txt"
    filename_with_reordering = "results_new_RO.txt"
    
    results_without_reordering = read_results_file(filename_without_reordering)
    results_with_reordering = read_results_file(filename_with_reordering)

    # print_results(results_without_reordering)

    # Plot the data for M=5 and Encode=True/False
    # plot_nodes_vs_n(results_with_reordering)
    # plot_nodes_vs_n(results_without_reordering, results_with_reordering)
    plot_nodes_vs_m(results_without_reordering, results_with_reordering)
