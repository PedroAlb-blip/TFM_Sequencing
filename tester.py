import subprocess

lis = ['nnnnnnnn','actgatctacagatca']
# lis.reverse()
print(lis[::-1])

print(str(lis[0:-1][4:-1]))

for i in lis:
    print(i.find("nnn"))