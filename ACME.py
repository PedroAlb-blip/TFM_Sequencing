#### The idea behind this file is to produce a list of results and values regarding the presence or absence of peaks.
import matplotlib.pyplot as plt
import numpy as np
import ast
from Functions import filterer, peak_discovery, rename

file = open("Sp105-VP7.txt", "r")
data = file.read()
all = list(filter(('').__ne__, data.split('\n')))
channels_fw = ast.literal_eval(all[10])
channels_rv = ast.literal_eval(all[12])
guide_fw = all[14]
guide_rv = all[16]
channels_fw = rename(channels_rv, guide_rv)
# ploc fw ast.literal_eval(all[6]), 

ploc_fw_l = peak_discovery(channels_fw)[1]
ploc_fw_d = peak_discovery(channels_fw)[0]
length = len(ploc_fw_l)
color=["red", "green", "yellow", "blue"]
maxes = []
nt = []
results = []

for i in range(0, length, 50):
    ins_1 = " "
    j=0
    try:
        for c in list(channels_fw.keys()):
            plt.plot(channels_fw[c][ploc_fw_l[i]:ploc_fw_l[i+50]], color=color[j], label=c)
            j=j+1
        plt.legend(loc="upper left")
        list_of_string = []
        for k in ploc_fw_l[i:i+50]:
            try:
                min_lst = [ploc_fw_d[c][0] for c in list(ploc_fw_d.keys())]
            except IndexError:
                rmv = [c for c in list(ploc_fw_d.keys()) if len(ploc_fw_d[c]) == 0]
                ploc_fw_d.pop(str(rmv[0]))
                min_lst = [ploc_fw_d[c][0] for c in list(ploc_fw_d.keys())]
            min_key = list(ploc_fw_d.keys())[min_lst.index(min(min_lst))]
            string = str(min_key) + str(k)
            list_of_string.append(string)
            plt.text(k - ploc_fw_l[i], channels_fw[min_key][k], string)
            ploc_fw_d[min_key].pop(0)
        plt.show(block=False)
        peak = 0
        while ins_1 != "" and ins_1 != "BREAK" and peak < len(list_of_string):
            print(list_of_string[peak])
            ins_1 = input()
            if ins_1 != "" and ins_1 != "ALL" and ins_1 != "NONE" and ins_1 != "BREAK":
                maxes.append(int(list_of_string[peak][1:]))
                nt.append(list_of_string[peak][0])
                results.append(int(ins_1))
            peak = peak + 1
        plt.close()
        if ins_1 == "BREAK":
            break
    except IndexError:
        break


#### Instead of removing the training parameters created by jiggling try expanding the results 
# new_params = filterer(maxes, channels_fw)
# redundance = [new_params[i][-1] for i in range(0,len(new_params))]
# jig_max = jiggler(maxes, channels_fw)[1]
# print(redundance, '\n', jig_max)
# for i in maxes:
#     rem = []
#     for j in redundance:
#         if j in range(i-3, i+3):
#             rem.append(j)
#             results_1.append(results[maxes.index(i)])
#     for k in rem:
#         redundance.pop(redundance.index(k))
# param = maxes[0:end] + new_params
# print(len(results_1), len(new_params))

rest = open("results.txt", "a")
rest.write(str(maxes))
rest.write('\n')
rest.write(str(nt))
rest.write('\n')
rest.write(str(results))

param_list = filterer(channels_fw, peak_discovery(channels_fw)[0])
print(param_list)
print(len(results))


train = open("training_params.txt", "a")
train.write('\n')
train.write(str(param_list))        
locs = open("locations.txt", "a")
locs.write('\n')
locs.write(str(maxes))
