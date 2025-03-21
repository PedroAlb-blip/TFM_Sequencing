import subprocess
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

def tmp_rm():
   for i in range(0, len(os.listdir())):
      check=str("aln_tmp_" + str(i))
      if os.path.exists(check):
         os.remove(check)
   return check

tmp_rm()

### This can be further automated by searching for -F_ or -R_ in the path
def reader (path, fmt="abi"):
   ### This function takes the path of the forward or reverse, extracts the sequences and produces an object containing the fasta sequence
   record = SeqIO.read(path, fmt)
   if "-R_" in path:
      sequence = str(record.annotations["abif_raw"]["PBAS1"]).replace("b","").replace("'","")[::-1]
      sequence=sequence.replace("A","Z").replace("G","X").replace("T","A").replace("C","G").replace("Z","T").replace("X","C")
      sequence = str(">Reverse" + '\n' + sequence + '\n')
   elif "-F_" in path:
      sequence = str(record.annotations["abif_raw"]["PBAS1"]).replace("b","").replace("'","")
      sequence = str(">Forward" + '\n' + sequence +'\n')
   return sequence

path_fw="c:\\Users\\Pedro\\Downloads\\secuenciasvp7_sp101bsp105sp106sp109sp111sp113sp116\\sec2025-016_60_Sp101b-VP7_RV-VP7-F_2025-02-24.ab1"
path_rv="c:\\Users\\Pedro\\Downloads\\secuenciasvp7_sp101bsp105sp106sp109sp111sp113sp116\\sec2025-016_91_Sp101b-VP7_RV-VP7-R_2025-02-24.ab1"

def tmp_aln(fw, rv):
   ### This function creates 
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

def aligner(file):
   cmd = [r"C:\Mafft\mafft.bat", "--clustalout", "--localpair", "--maxiterate", "1000", "--op", "2.0", "--ep", "0.1", file]
   process = subprocess.run(cmd, capture_output=True, text=True)
   matches = process.stdout.split('\n')
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


aligner(tmp_aln(reader(path_fw), reader(path_rv)))

# print(record_F.annotations.keys())
# dict_keys(["dye", "abif_raw", "sample_well", "run_finish", "machine_model", "run_start", "polymer"])
# list(record_F.annotations["abif_raw"].keys())
# dict_keys(["DATA5", "DATA8", "RUNT1", "phAR1", ..., "DATA6"])
channels = ["DATA9", "DATA10", "DATA11", "DATA12"]
trace_F = defaultdict(list)
trace_R = defaultdict(list)
record_F = SeqIO.read(path_fw, "abi")
record_R = SeqIO.read(path_rv, "abi")
for c in channels:
    trace_F[c] = record_F.annotations["abif_raw"][c]
    trace_R[c] = record_R.annotations["abif_raw"][c]


#### Extract the normalised data for the chromatogram from the .ab1 file and create a list with it
