




import subprocess
import os
import numpy as np
import itertools
import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import defaultdict

def tmp_rm():
   for i in range(0, len(os.listdir())):
      check=str("aln_tmp_" + str(i))
      if os.path.exists(check):
         os.remove(check)
   return 

tmp_rm()

### This can be further automated by searching for -F_ or -R_ in the path
def reader (path, fmt="abi"):
   ### This function takes the path of the forward or reverse, extracts the sequences and produces an object containing the fasta sequence
   record = SeqIO.read(path, fmt)
   c = ["DATA9", "DATA10", "DATA11", "DATA12"]
   channels = dict()
   guide = str(record.annotations["abif_raw"]["FWO_1"]).replace("b","").replace("'", "")
   if "-R_" in path:
      sequence = str(record.annotations["abif_raw"]["PBAS2"]).replace("b","").replace("'","")[::-1]
      sequence = sequence.replace("A","Z").replace("G","X").replace("T","A").replace("C","G").replace("Z","T").replace("X","C")
      sequence = str(">Reverse" + '\n' + sequence + '\n')
      for i in c:
         channels[i] = record.annotations["abif_raw"][i][::-1]
      ploc = np.subtract(list(itertools.repeat(len(channels["DATA9"]), len(record.annotations["abif_raw"]["PLOC2"]))), list(record.annotations["abif_raw"]["PLOC2"][::-1]))
      
   elif "-F_" in path:
      sequence = str(record.annotations["abif_raw"]["PBAS2"]).replace("b","").replace("'","")
      sequence = str(">Forward" + '\n' + sequence +'\n')
      for i in c:
         channels[i] = record.annotations["abif_raw"][i] 
      ploc = record.annotations["abif_raw"]["PLOC2"]
   return sequence, ploc, channels, guide

def tmp_aln(fw, rv):
   ### This function creates a file with both sequences
   files=os.listdir()
   for i in range(0, len(os.listdir())):
      check=str("aln_tmp_" + str(i))
      if check in files:
         continue
      else:
         new_tmp=open(check, 'x')
         new_tmp.write(fw)
         new_tmp.write(rv)
         new_tmp.close()
         break
   return check

# tmp_aln(reader(path_fw), reader(path_rv, True))

def aligner(file, v=False):
   cmd = [r"C:\Mafft\mafft.bat", "--clustalout", "--localpair", "--maxiterate", "1000", "--op", "2.0", "--ep", "0.1", file]
   process = subprocess.run(cmd, capture_output=True, text=True)
   matches = process.stdout.split('\n')
   if v==True:
      print(process.stdout)
   fw_seq=""
   rv_seq=""
   match_seq=""
   for i in matches:
      if "Forward" in i:
         fw_seq=fw_seq+(i[16:])
      elif "Reverse" in i:
         rv_seq=rv_seq+(i[16:])
      elif len(i)>50:
         match_seq=match_seq+(i[16:])
   match_seq=match_seq.replace(" ", "_")
   return match_seq, fw_seq, rv_seq


def locator(path_fw, path_rv):
   all=aligner(tmp_aln(reader(path_fw)[0], reader(path_rv)[0]))
   align=all[0]
   fw=all[1]
   rv=all[2]
   Forward=reader(path_fw)[0].removeprefix(">Forward\n")
   Reverse=reader(path_rv)[0].removeprefix(">Reverse\n")
   i=0
   while i < (len(align)-20):
      kmer=align[i:i+20]
      if kmer.count('*') > 6:
         pos_1=Forward.index(fw[i + kmer.index('*'):i + 25].upper())
         # print(fw[i + kmer.index('*'):i + 25].upper())
         pos_2=Reverse.index(rv[i + kmer.index('*'):i + 25].upper())
         # print(rv[i + kmer.index('*'):i + 25].upper())
         break
      else:
         i=i+1
   return pos_1, pos_2

##### VERY IMPORTANT: POSITIONS OF PLOC CHANGE RELATIVE TO THE SEQUENCE WHEN DOING A REVERSE COMPLEMENT 
##### Perchance doing the absolute value of the subtraction between the PLOCs and the total length 14000

path_fw="c:\\Users\\Pedro\\Downloads\\secuenciasvp7_sp101bsp105sp106sp109sp111sp113sp116\\sec2025-016_60_Sp101b-VP7_RV-VP7-F_2025-02-24.ab1"
path_rv="c:\\Users\\Pedro\\Downloads\\secuenciasvp7_sp101bsp105sp106sp109sp111sp113sp116\\sec2025-016_91_Sp101b-VP7_RV-VP7-R_2025-02-24.ab1"
print(locator(path_fw, path_rv))

