from Functions import rename, peak_discovery, filterer
import ast
import numpy as np
# file_names = ["Sp101b-VP7.txt", "Sp105-VP7.txt", "Sp106-VP7.txt", "Sp109-VP7a.txt", "Sp111-VP7a.txt"]
# file_app = open("training_params.txt", "a")
# for i in file_names:
#     data = open(i, "r")
#     read = data.read()
#     all = list(filter(('').__ne__, read.split('\n')))
#     channels_fw = ast.literal_eval(all[10])
#     channels_rv = ast.literal_eval(all[12])
#     guide_fw = all[14]
#     guide_rv = all[16]
#     channels_fw = rename(channels_fw, guide_fw)
#     channels_rv = rename(channels_rv, guide_rv)
#     file_app.write(str(filterer(channels_fw, peak_discovery(channels_fw)[0])[0]))
#     file_app.write("\n")
#     file_app.write(str(filterer(channels_rv, peak_discovery(channels_rv)[0])[0]))
#     file_app.write("\n")
#     data.close()
# file_app.close()

