import sys
from Bio import SeqIO
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

# handle_F=open("c:\\Users\\Pedro\\Downloads\\secuenciasvp7_sp101bsp105sp106sp109sp111sp113sp116\\sec2025-016_60_Sp101b-VP7_RV-VP7-F_2025-02-24.ab1", "rb")
# handle_R=open("c:\\Users\\Pedro\\Downloads\\secuenciasvp7_sp101bsp105sp106sp109sp111sp113sp116\\sec2025-016_91_Sp101b-VP7_RV-VP7-R_2025-02-24.ab1", "rb")
# readable_F=open("c:\\Users\\Pedro\\Downloads\\Read_F","a")
# readable_R=open("c:\\Users\\Pedro\\Downloads\\Read_R","a")
# for record in SeqIO.parse(handle_F, "abi"):
#     readable_F.write(str(record))
# for record in SeqIO.parse(handle_R, "abi"):
#     readable_R.write(str(record))

record = SeqIO.read("c:\\Users\\Pedro\\Downloads\\secuenciasvp7_sp101bsp105sp106sp109sp111sp113sp116\\sec2025-016_60_Sp101b-VP7_RV-VP7-F_2025-02-24.ab1", "abi")
print(record.annotations.keys())
# dict_keys(["dye", "abif_raw", "sample_well", "run_finish", "machine_model", "run_start", "polymer"])
list(record.annotations["abif_raw"].keys())
# dict_keys(["DATA5", "DATA8", "RUNT1", "phAR1", ..., "DATA6"])
channels = ["DATA9", "DATA10", "DATA11", "DATA12"]
trace = defaultdict(list)
for c in channels:
    trace[c] = record.annotations["abif_raw"][c]
peaks = record.annotations["abif_raw"]["PLOC2"]

# Values between 0 and 1200
# plt.plot(trace["DATA1"], color="blue")
# plt.plot(trace["DATA2"], color="red")
# plt.plot(trace["DATA3"], color="green")
# plt.plot(trace["DATA4"], color="yellow")
# plt.show()

# print(record.annotations["abif_raw"].keys())
# print(record.annotations["abif_raw"]["PLOC1"])

try: 
    assert len(trace["DATA1"]) == len(trace["DATA2"]) == len(trace["DATA3"]) == len(trace["DATA4"])
except AssertionError: 
    print("Different lengths of input")
    sys.exit(1)

dat1=[]
dat2=[]
dat3=[]
dat4=[]
for i in range(0,len(trace["DATA1"])):
    mean = (trace["DATA1"][i] + trace["DATA2"][i] + trace["DATA3"][i] + trace["DATA4"][i])/4
    sd = np.sqrt(((trace["DATA1"][i]-mean)**2 +(trace["DATA2"][i]-mean)**2 + (trace["DATA3"][i]-mean)**2 + (trace["DATA4"][i]-mean)**2)/3)
    dat1.append(np.exp((trace["DATA1"][i]-mean)/sd))
    dat2.append(np.exp((trace["DATA2"][i]-mean)/sd))
    dat3.append(np.exp((trace["DATA3"][i]-mean)/sd))
    dat4.append(np.exp((trace["DATA4"][i]-mean)/sd))

# plt.plot(dat1[2300:3032], color="blue")
# plt.plot(dat2[2300:3032], color="red")
# plt.plot(dat3[2300:3032], color="green")
# plt.plot(dat4[2300:3032], color="yellow")

plt.plot(trace["DATA9"][2300:3032], color="blue")
plt.plot(trace["DATA10"][2300:3032], color="red")
plt.plot(trace["DATA11"][2300:3032], color="green")
plt.plot(trace["DATA12"][2300:3032], color="yellow")

print(len(trace["DATA10"]))

peaks_purge = [num -2300 for num in peaks if num >= 2300 and num <= 3032]
for i in peaks_purge:
    maxes=max(trace["DATA9"][i+2300], trace["DATA10"][i+2300], trace["DATA11"][i+2300], trace["DATA12"][i+2300])
    plt.plot(i, maxes, "co")

### Proof of concept that the peak locations refer to the x axis of the channels 

plt.show()