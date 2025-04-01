#### The idea behind this file is to produce a list of results and values regarding the presence or absence of peaks.
import matplotlib.pyplot as plt
import numpy as np
import ast
from functions import jiggler, cross_check, peak_discovery
file = open("Sp101b-VP7.txt", "r")
data = file.read()
all = list(filter(('').__ne__, data.split('\n')))
channels_fw = ast.literal_eval(all[10])
channels_rv = ast.literal_eval(all[12])
ploc_fw = jiggler(ast.literal_eval(all[6]), channels_fw)
ploc_rv = jiggler(ast.literal_eval(all[8]), channels_rv)
length = len(ploc_fw)
color=["red", "green", "yellow", "blue"]
for i in range(0, length, 50):
    j=0
    maxes = {}
    try:
        for c in list(channels_fw.keys()):
            plt.plot(channels_fw[c][ploc_fw[i]:ploc_fw[i+50]], color=color[j])
            j=j+1
        for k in ploc_fw[i:i+50]:
            string = str(k)
            plt.text(k - ploc_fw[i], max([channels_fw[c][k] for c in channels_fw.keys()]), string)
        plt.show()
    except IndexError:
        break

        