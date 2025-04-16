import ast
from functions import jiggler, confidence, filterer

file = open("Sp101b-VP7.txt", "r")
data = file.read()
all = list(filter(('').__ne__, data.split('\n')))
channels_fw = ast.literal_eval(all[10])
channels_rv = ast.literal_eval(all[12])
ploc_fw = jiggler(ast.literal_eval(all[6]), channels_fw)
ploc_rv = jiggler(ast.literal_eval(all[8]), channels_rv)


print(filterer(channels_fw, list(ploc_fw[1])))

length = len(ploc_rv)
locs = ast.literal_eval(all[24])
align = all[18]





# record = SeqIO.read("c:\\Users\\Pedro\\Downloads\\secuenciasvp7_sp101bsp105sp106sp109sp111sp113sp116\\sec2025-016_91_Sp101b-VP7_RV-VP7-R_2025-02-24.ab1","abi")
# c = ["DATA9", "DATA10", "DATA11", "DATA12"]
# channels = dict()
# for i in c:
#     channels[i] = record.annotations["abif_raw"][i][::-1]
#     print(len(channels[i]))
# ploc = np.subtract(list(itertools.repeat(len(channels["DATA9"]), len(record.annotations["abif_raw"]["PLOC2"]))), list(record.annotations["abif_raw"]["PLOC2"][::-1]))

# plt.plot(channels["DATA9"][4000:4500], color = "blue")
# plt.plot(channels["DATA10"][4000:4500], color = "green")
# plt.plot(channels["DATA11"][4000:4500], color = "red")
# plt.plot(channels["DATA12"][4000:4500], color = "yellow")
# peaks_purge = [num - 4000 for num in ploc if num >= 4000 and num <= 4500]
# for i in peaks_purge:
#     maxes=max(channels["DATA9"][i+4000], channels["DATA10"][i+4000], channels["DATA11"][i+4000], channels["DATA12"][i+4000])
#     plt.plot(i, maxes, "co")

# plt.show()