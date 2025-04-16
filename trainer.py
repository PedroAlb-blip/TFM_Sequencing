#### The idea behind this file is to produce a list of results and values regarding the presence or absence of peaks.
import matplotlib.pyplot as plt
import numpy as np
import ast
from functions import jiggler, filterer

file = open("Sp101b-VP7.txt", "r")
data = file.read()
all = list(filter(('').__ne__, data.split('\n')))
channels_fw = ast.literal_eval(all[10])
channels_rv = ast.literal_eval(all[12])
ploc_fw = jiggler(ast.literal_eval(all[6]), channels_fw)[1]
ploc_rv = jiggler(ast.literal_eval(all[8]), channels_rv)[1]
length = len(ploc_rv)
locs = ast.literal_eval(all[24])
color=["red", "green", "yellow", "blue"]
maxes_file = open("training_params.txt", "r")
maxes_symb = maxes_file.read()
maxes = list(ast.literal_eval(maxes_symb))
results_file = open("results.txt", "r")
results_symb = results_file.read()
results = list(ast.literal_eval(results_symb))
end = len(maxes)
for i in range(0, length, 50):
    ins_1 = " "
    ins_2 = " "
    j=0
    try:
        for c in list(channels_rv.keys()):
            plt.plot(channels_rv[c][ploc_rv[i]:ploc_rv[i+50]], color=color[j])
            j=j+1
        for k in ploc_rv[i:i+50]:
            string = str(k)
            plt.text(k - ploc_rv[i], max([channels_rv[c][k] for c in channels_rv.keys()]), string)
        plt.show(block=False)
        while ins_1 != "":
            ins_1 = input()
            ins_2 = input()
            if ins_1 != "" and ins_1 != "ALL":
                maxes.append(int(ins_1))
                results.append(int(ins_2))
            if ins_1 == "ALL":
                maxes = maxes + ploc_rv[i:i+50]
                results = results + 50*[1]
        plt.close()
    except IndexError:
        break
param = maxes[0:end] + filterer(channels_rv, maxes[end + 1:-1])

rest = open("results.txt", "w")
rest.write(str(results))
train = open("training_params.txt", "w")
train.write(str(param))        