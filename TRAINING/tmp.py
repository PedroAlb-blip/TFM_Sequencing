from Functions import rename, peak_discovery, filterer
import ast
import numpy as np
file_names = ["TRAINING\\Sp101b-VP7.txt", "TRAINING\\Sp105-VP7.txt", "TRAINING\\Sp106-VP7.txt", "TRAINING\\Sp109-VP7a.txt", "TRAINING\\Sp111-VP7a.txt", "TRAINING\\Sp113-VP7a.txt", "TRAINING\\Sp116-VP7a.txt"]
file_app = open("TRAINING\\training_params.txt", "a")
for i in file_names:
    data = open(i, "r")
    read = data.read()
    all = list(filter(('').__ne__, read.split('\n')))
    channels_fw = ast.literal_eval(all[10])
    channels_rv = ast.literal_eval(all[12])
    guide_fw = all[14]
    guide_rv = all[16]
    channels_fw = rename(channels_fw, guide_fw)
    channels_rv = rename(channels_rv, guide_rv)
    file_app.write(str(filterer(channels_fw, peak_discovery(channels_fw)[0])))
    file_app.write("\n")
    file_app.write(str(filterer(channels_rv, peak_discovery(channels_rv)[0])))
    file_app.write("\n")
    data.close()
file_app.close()


file_res = open("TRAINING\\results.txt", "r")
file_par = open("TRAINING\\training_params.txt", "r")
lines_res = file_res.read()
lines_par = file_par.read()
results = list(filter(('').__ne__, lines_res.split('\n')))
params = list(filter(('').__ne__, lines_par.split('\n')))
ones = []
train = []
for i in range(0, len(results)):
    print("\t", i)
    if results.index(results[i]) % 3 == 1:
        res = list(ast.literal_eval(results[i]))
        par_lst = list(ast.literal_eval(params[int(round(i/3, 0))]))
        for j in res:
            ok = ""
            while j != par_lst[res.index(j)][-1] and ok != "ok":
                print(j, par_lst[res.index(j) - 1][-1], par_lst[res.index(j) + 1][-1])
                print(res.index(j))
                ok = input()
    elif results.index(results[i]) % 3 == 0:
        ones = ones + list(ast.literal_eval(results[i]))
        train = train + list(ast.literal_eval(params[int(round(i/3,0)) - 1]))
        if len(ones) < len(train):
            ones = ones + [1]*(len(train)-len(ones))