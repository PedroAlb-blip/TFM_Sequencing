import ast
import re
import subprocess
import matplotlib.pyplot as plt
from Functions import jiggler, confidence, filterer, rename

file = open("results.txt", "r")
data = file.read()
file.close()
seqtxt = open("Sp101b-VP7.txt", "r")
seqrd = seqtxt.read()
vals = list(filter(('').__ne__, seqrd.split('\n')))
channels_fw = rename(ast.literal_eval(vals[10]), vals[14])
channels_rv = rename(ast.literal_eval(vals[12]), vals[16])
all = list(filter(('').__ne__, data.split('\n')))
peak_fw = ast.literal_eval(all[2])
peak_rv = ast.literal_eval(all[5])
seq_fw = ast.literal_eval(all[1])
seq_rv = ast.literal_eval(all[4])
ploc_fw = ast.literal_eval(all[0])
ploc_rv = ast.literal_eval(all[3])
seq_1 = ""
ploc_1 =[]
index = 0
for i in peak_fw:
    if i == 0:
        seq_1 = seq_1 + seq_fw[index]
        ploc_1.append(ploc_fw[index])
    index = index + 1
seq_2 = ""
ploc_2 = []
index = 0
for i in peak_rv:
    if i == 0:
        seq_2 = seq_2 + seq_rv[index]
        ploc_2.append(ploc_rv[index])
    index = index + 1

align = f">\n{seq_1}\n>\n{seq_2}"

cmd = [r"C:\Mafft\mafft.bat", "--clustalout", "--localpair", "--maxiterate", "1000", "--op", "2.0", "--ep", "0.1", "-"]
process = subprocess.run(
    cmd,
    input=align,
    capture_output=True,
    text=True,
    shell=True  # May be needed on Windows
)
align = process.stdout.split('\n')
align = [i for i in align if len(i) > 1]
aln = align.copy()
only = ""
seq_fw_al = ""
seq_rv_al = ""
for i in aln:
    if align.index(i) % 3 == 0 and align.index(i) != 0:
        only = only + i[16:]
        align[align.index(i)] = 0
    elif align.index(i) % 3 == 1:
        seq_fw_al = seq_fw_al + i[16:]
        align[align.index(i)] = 0
    elif align.index(i) % 3 == 2:
        seq_rv_al = seq_rv_al + i[16:]
        align[align.index(i)] = 0
# print(seq_2, seq_rv_al.replace("-", "").upper())
only = only.replace(" ", "_").replace(".","_")
ylno = list(only)
ylno.reverse()
ylno = ''.join(ylno)
indices = (only.find("*"), len(ylno) - ylno.find("*"))
corrections = [ (i.start(), i.end()) for i in list(re.finditer('\*_+\*', only))]
color=["red", "green", "yellow", "blue"]
for i in corrections:
    fw = seq_fw_al[i[0] - 3 :i[1] + 3].replace("-","").upper() ## fw is the mismatched sequence for the forward sequence in a specific region of incongruence
    rv = seq_rv_al[i[0] - 3 :i[1] + 3].replace("-","").upper() 
    attribute_fw = [list(channels_fw[c][ploc_1[seq_1.find(fw)] : ploc_1[seq_1.find(fw)+len(fw)]]) for c in list(channels_fw.keys())] ## These are the electropherogram values within the limits of mismatch
    attribute_rv = [list(channels_rv[c][ploc_2[seq_2.find(rv)] : ploc_2[seq_2.find(rv)+len(rv)]]) for c in list(channels_rv.keys())] ## The order of the guides is different
    attribute_rv.reverse()
    MPD_fw = [ploc_1[seq_1.find(fw) + 1 + k] - ploc_1[seq_1.find(fw) + k] for k in range(0,len(fw))]
    MPD_fw.sort()
    MPD_rv = [ploc_2[seq_2.find(rv) + 1 + k] - ploc_2[seq_2.find(rv) + k] for k in range(0,len(rv))]
    MPD_rv.sort()
    if MPD_rv[int(round(len(MPD_rv)/2,0))] - MPD_fw[int(round(len(MPD_fw)/2,0))] > 0:
        print(len(attribute_fw[0]))
        suml = []
        itr = 0
        for i in range(0, len(attribute_fw[0])):
            suma = 0
            for c in range(0, len(attribute_fw)):
                suma = suma + attribute_fw[c][i]
            suml.append(suma)
            if len(suml) > len(attribute_fw[0])/len(fw):
                diff = 0
                c = 0
                while c in range(0, len(attribute_fw)) and diff < MPD_rv[int(round(len(MPD_rv)/2,0))] - MPD_fw[int(round(len(MPD_fw)/2,0))]:
                    attribute_fw[c].insert(int(round(suml.index(min(suml)) + itr*len(attribute_fw[0])/len(fw),0)), (attribute_fw[c][suml.index(min(suml)) + int(round(itr*len(attribute_fw[0])/len(fw),0)) -1] + attribute_fw[c][1 + suml.index(min(suml)) + int(round(itr*len(attribute_fw[0])/len(fw),0))])/2)
                    c = c + 1
                    diff = diff + 1
                itr=itr+1
                suml = []
        print(len(attribute_fw[0]))
    elif  MPD_rv[int(round(len(MPD_rv)/2,0))] - MPD_fw[int(round(len(MPD_fw)/2,0))] < 0:
        print(len(attribute_rv[0]))
        suml = []
        itr = 0
        for i in range(0, len(attribute_rv[0])):
            suma = 0
            for c in range(0, len(attribute_rv)):
                suma = suma + attribute_rv[c][i]
            suml.append(suma)
            if len(suml) > len(attribute_rv[0])/len(rv):
                diff = 0
                c = 0
                while c in range(0, 4) and diff < MPD_fw[int(round(len(MPD_fw)/2,0))] -  MPD_rv[int(round(len(MPD_rv)/2,0))]:
                    attribute_rv[c].insert(int(round(suml.index(min(suml)) + itr*len(attribute_rv[0])/len(rv),0)), (attribute_rv[c][suml.index(min(suml)) + int(round(itr*len(attribute_rv[0])/len(rv),0)) -1] + attribute_rv[c][1 + suml.index(min(suml)) + int(round(itr*len(attribute_rv[0])/len(rv),0))])/2)
                    c = c + 1
                    diff = diff + 1
                itr=itr+1
                suml = []
        print(len(attribute_rv[0]))
    for j in range(0, 4):
        attribute_comb = [attribute_fw[j][n] + attribute_rv[j][n] for n in range(0, min(len(attribute_rv[j]), len(attribute_fw[j])))]
        plt.plot(attribute_comb, color=color[j])
        print(attribute_fw[j][min(len(attribute_rv[j]), len(attribute_fw[j])):])
        print(attribute_rv[j][min(len(attribute_rv[j]), len(attribute_fw[j])):])
    plt.show()
