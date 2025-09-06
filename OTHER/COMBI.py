import matplotlib.pyplot as plt
import numpy as np

file = open("OTHER\\Combinations.txt", "r")
content = file.read()
split_1 = [i for i in content.split("\n") if len(i) != 0]
genotype = []
farm = []
perc = []
for i in split_1:
    farm.append(i.split(">")[1].split("_")[0] + i.split(">")[1].split("_")[2][0] + i.split(">")[2].split("_")[2][0])
    genotype.append([i.split(">")[1].split("_")[1], i.split(">")[2].split("_")[1]])
    perc.append([float(i.split(">")[1].split(" ")[1].replace("\t","").replace("%",""))/100, float(i.split(">")[2].split(" ")[1].replace("\t","").replace("%",""))/100])
g_key = ['G4', 'G5', 'G9', 'G1', 'G3', 'GX']
p_key = ['P7', 'P13', 'P23', 'P6', 'P19', 'P32', 'PX']
g_dict = {key : [] for key in g_key}
p_dict = {key : [] for key in p_key}
intersec = {}
for i in range(0, len(farm)):
    g_dict[genotype[i][0]].append(farm[i])
    p_dict[genotype[i][1]].append(farm[i])
for key_1 in g_key:
    for key_2 in p_key:
        for i in g_dict[key_1]:
            if i in p_dict[key_2]:
                try:
                    intersec[key_1+key_2] = intersec[key_1+key_2] + [i]
                except KeyError:
                    intersec.update({key_1+key_2 : [i]})
fig, ax = plt.subplots()
fig_str, ax_str = plt.subplots()
arrowprops=dict(arrowstyle='<-', color='black', linewidth=1, mutation_scale=2)
for i in range(0,len(farm)):
    if farm[i].startswith("Sp"):
        col = "violet"
    if farm[i].startswith("It"):
        col = "green"
    if farm[i].startswith("Pl"):
        col = "red"
    offset_x = 0.8 * np.cos(6.28*(intersec[genotype[i][0]+genotype[i][1]].index(farm[i])+1)/len(intersec[genotype[i][0]+genotype[i][1]]))
    offset_y = 0.9 * np.sin(6.28*(intersec[genotype[i][0]+genotype[i][1]].index(farm[i])+1)/len(intersec[genotype[i][0]+genotype[i][1]]))
    if perc[i][0] > 0 and perc[i][1] > 0:
        ax.scatter(g_key.index(genotype[i][0])*2 + (1 - perc[i][0])*4, p_key.index(genotype[i][1])*2 + (1 - perc[i][1])*4, s=60, c=col, alpha=0.5)
        ax.annotate(farm[i], xy = (g_key.index(genotype[i][0])*2 + (1 - perc[i][0])*4, p_key.index(genotype[i][1])*2 + (1 - perc[i][1])*4), xytext= (g_key.index(genotype[i][0])*2 + offset_x, p_key.index(genotype[i][1])*2 + offset_y), fontsize = 5, arrowprops=arrowprops)
    else:
        ax.scatter(g_key.index(genotype[i][0])*2, p_key.index(genotype[i][1])*2, s=60, c=col, alpha=0.5)
        ax.annotate(farm[i], xy = (g_key.index(genotype[i][0])*2, p_key.index(genotype[i][1])*2), xytext= (g_key.index(genotype[i][0])*2 + offset_x, p_key.index(genotype[i][1])*2 + offset_y), fontsize = 5, arrowprops=arrowprops)
    ax_str.scatter(g_key[g_key.index(genotype[i][0])], p_key[p_key.index(genotype[i][1])], s=4, c=col)
plt.show()



