#### The idea behind this file is to produce a list of results and values regarding the presence or absence of peaks.
import matplotlib.pyplot as plt
import numpy as np
import ast
from functions import jiggler, filterer, peak_discovery, rename

file = open("Sp105-VP7.txt", "r")
data = file.read()
all = list(filter(('').__ne__, data.split('\n')))
channels_fw = ast.literal_eval(all[10])
channels_rv = ast.literal_eval(all[12])
guide_fw = all[14]
guide_rv = all[16]
print(guide_rv, guide_fw)
channels_rv = rename(channels_rv, guide_rv)
# ploc fw ast.literal_eval(all[6]), 

ploc_rv_l = peak_discovery(channels_rv)[1]
ploc_rv_d = peak_discovery(channels_rv)[0]
length = len(ploc_rv_l)
locs = ast.literal_eval(all[24])
color=["red", "green", "yellow", "blue"]
maxes = []
nt = []
results = []

for i in range(0, length, 50):
    ins_1 = " "
    j=0
    try:
        for c in list(channels_rv.keys()):
            plt.plot(channels_rv[c][ploc_rv_l[i]:ploc_rv_l[i+50]], color=color[j], label=c)
            j=j+1
        plt.legend(loc="upper left")
        list_of_string = []
        for k in ploc_rv_l[i:i+50]:
            min_lst = [ploc_rv_d[c][0] for c in list(channels_rv.keys())]
            min_key = list(channels_rv.keys())[min_lst.index(min(min_lst))]
            string = str(min_key) + str(k)
            list_of_string.append(string)
            plt.text(k - ploc_rv_l[i], channels_rv[min_key][k], string)
            ploc_rv_d[min_key].pop(0)
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

param_list = filterer(channels_rv, peak_discovery(channels_rv)[0])
print(len(param_list[1]))
print(len(results))


train = open("training_params.txt", "a")
train.write('\n')
train.write(str(param_list))        
locs = open("locations.txt", "a")
locs.write('\n')
locs.write(str(maxes))