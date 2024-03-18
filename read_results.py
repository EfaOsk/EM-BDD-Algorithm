import matplotlib.pyplot as plt
from itertools import groupby
import statistics
import numpy as np
from scipy.stats import ttest_rel


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
        'Run': run,
        'Lorg': l_org,
        'num iter': num_iter,
        'Lfinal': l_final,
        'runtime': runtime
    }

    return result

def read_file_to_dict(filename, tmp):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            if '|' not in line:
                data_dict = parse_line_to_dict(line)
                data_dict['lr']= list()
                if tmp:
                    logfilename = f"experiments\experiments_10_BDD\experiments_10_BDD\{data_dict['name']}_Obs0_Hypo{data_dict['Run']}\logs\log.txt"
                else:
                    logfilename = f"experiments\experiments_10_BW\experiments_10_BW\{data_dict['name']}_Obs0_Hypo{data_dict['Run']}\logs\log.txt"
                with open(logfilename, 'r') as logs:
                    for log in logs:
                        data_dict['lr'].append(parse_log_file_line(log))
                data.append(data_dict)

    return data

def parse_log_file_line(line):
    parts = line.split(',')
    return float(parts[3].split(':')[1].strip())

def average_learning_rate(data):
    max_length = max(len(d['lr']) for d in data)
    avg_lr = []
    for i in range(max_length):
        rates = [data[k]['lr'][i] for k in range(len(data)) if i < len(data[k]['lr'])]
        avg_lr.append(np.mean(rates))
    return avg_lr



if __name__ == '__main__':
    experiment_folder = "experiments\experiments_10_BW\experiments_10_BW/results.txt"

    dataBW = read_file_to_dict(experiment_folder, False)
    experiment_folder = "experiments\experiments_10_BDD\experiments_10_BDD/results.txt"

    dataBDD = read_file_to_dict(experiment_folder, True)



    filename = "tmp_BW.txt"
    i = 0
    with open(filename, 'r') as file:
        for line in file:
            parts = line.split(',')
            dataBW[i]["test Lorg"] = float(parts[0].split(':')[1].strip())
            dataBW[i]["test Lfin"] = float(parts[1].split(':')[1].strip())
            i+=1
        

    filename = "tmp_BDD.txt"
    i = 0
    with open(filename, 'r') as file:
        for line in file:
            parts = line.split(',')
            dataBDD[i]["test Lorg"] = float(parts[0].split(':')[1].strip())
            dataBDD[i]["test Lfin"] = float(parts[1].split(':')[1].strip())
            i+=1

    models = dict()

    for i in range(50):
        model_name = dataBDD[i]['name']
        if model_name not in models:
            models[model_name] = {"test Lorg": [], "test LfinBW": [], "test LfinBDD": []}

        models[model_name]["test Lorg"].append(dataBDD[i]["test Lorg"])
        models[model_name]["test LfinBDD"].append(dataBDD[i]["test Lfin"])
        models[model_name]["test LfinBW"].append(dataBW[i]["test Lfin"])

    # for i in range(50):
    #     plt.figure()
    #     plt.plot(range(1, 1+len(dataBW[i]['lr'])), dataBW[i]['lr'])
    #     plt.plot(range(1,1+len(dataBDD[i]['lr'])), dataBDD[i]['lr'])
    #     plt.hlines(dataBDD[i]['Lorg'], 0, len(dataBDD[i]['lr']),  color='r', linestyle='--')
    #     plt.show()
        
    
    # models_timeBW = [list(),]*10
    # models_timeBDD = [list(),]*10
    # for i in range(10):
    #     for j in range(5):
    #         models_timeBW[i]+= dataBW[i*5+j]['lr']
    #         models_timeBDD[i]+= dataBDD[i*5+j]['lr']
    #     print(f"Model {i} & {statistics.mean(models_timeBW[i]):.4f} ( {statistics.stdev(models_timeBW[i]):.4f}) & {statistics.mean(models_timeBDD[i]):.4f} ( {statistics.stdev(models_timeBDD[i]):.4f})")

    models_timeBW = list() #[list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list()]
    models_timeBDD = list() #[list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list()]
    for i in range(10):
        for j in range(5):
            tmp = min(len(dataBW[i*5+j]['lr']), len(dataBDD[i*5+j]['lr']), 30)
            models_timeBW+=dataBW[i*5+j]['lr'][:tmp-1]
            models_timeBDD+= dataBDD[i*5+j]['lr'][:tmp-1]
        
        # Perform an independent t-test
    t_statistic, p_value = ttest_rel(models_timeBW, models_timeBDD)


        # Results
    print(f"model {i}, t-statistic {t_statistic}, p-value {p_value}")
    
    # print(models_timeBW)
    # print(models_timeBDD)



    avg_lr_bdd = average_learning_rate(dataBDD)[:60]
    avg_lr_bw = average_learning_rate(dataBW)[:60]

    # for i in range(len(avg_lr_bw)):
    #     print(f"\t{i+1} & {avg_lr_bw[i]:.3f} & {avg_lr_bdd[i]:.3f}\\\\hline")
    #     if (i%2):
    #         print("\t\\rowcolor{lightgray}")
    # Plotting the average learning rates
    # plt.figure(figsize=(10, 6))
    # plt.plot(avg_lr_bdd, label='EM-BDD', linestyle='--')
    # plt.plot(avg_lr_bw, label='BW', linestyle='--')
    # plt.xlabel('Iteration')
    # plt.ylabel('Average Learning Rate')
    # plt.title('Average Learning Rate per Iteration for EM-BDD and BW Algorithms')
    # plt.legend()
    # # plt.grid(True)
    # plt.show()


#     print(models)    
#     for m in models:
#         # print(f'{m} & {statistics.mean(models[m]["test Lorg"]):.3f} & {statistics.mean(models[m]["test LfinBW"]):.3f} & {statistics.mean(models[m]["test LfinBDD"]) -statistics.mean(models[m]["test LfinBW"]):.3f}\\\\ ')
#         for i in range(5):
#             print(f'{m}, {abs(models[m]["test Lorg"][i]-models[m]["test LfinBW"][i])}, {abs(models[m]["test Lorg"][i]-models[m]["test LfinBDD"][i]) }')
        

    
