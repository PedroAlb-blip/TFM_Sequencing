import subprocess

lis = ['nnnnnnnn','actgatctacagatca']
# lis.reverse()
print(lis[::-1])

real_sq=open("aln_tmp_0", "r")
douce=real_sq.readlines()
Forward=douce[1]
Reverse=douce[3]
align = '___________________________________________________________________________________________________________________________________________________________________________________*._****_*.***__.*___*_**********_*_**_*_*************_***_****_************_***************************************_*__**__******************************_*************_*******************_***************************_*************_**********************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************_***********************************************************_***********************_*****_***********_********************_*********_********_******_*********_**_*__*__._*__*_*****_*_*_**_**_*_**_**_****__*_*__**._.__*_***_*_*__*.._*_**__________________________________________________________________________________________________________________________________________________________________________________________________'
fw='---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------cnggaagctccttttatgt---atggtattgaatataccacaattctaacttttttgatatcgcttgtatttataaa-ctatatattgaaatcagtaactagaacaatggactttataa--tttacagatttctactggttatagttgtacttgcaccatttataaaaa-cccaaaattatggaattaa-tttgccgataactggatccatggatacaccgtatgcaaatt-ctacagcaagtgagacatttttaacttcaacattatgtttatattatccgaatgaggcagctactgaaattgcagacagtaaatggacagagacattgtcacagttgtttttaacgaaaggatggccgacaggttcagtttactttaaaggttatgcagatattgcatcattttctgtagaaccgcagttatactgtgactataacattgtattaatgaaatatgatgtaagcttgcaattagatatgtctgaattggctgatctaatattaaatgaatggttatgcaatccaatggatataacgctatattattatcaacaaactgatgaggcgaacaaatggatatctatgggttcttcatgtacaattaaagtatgtcccctaaatacacaaacccttggaataggatgttcaaccacagacactaactcatttgaaatggtggctaatgcagagaagttagttataacagatgttgtcgatggagtcaatcacaaactggacgtaacgacaaacacgtgtacnatacgaaattgtaaaaaacttggaccaagggaaaacgttgctgtaattcaggtaggaggnccagatgtacttgatataactgcngatccnactactgcgccccagactgaaagaatgatgcgnataaactggnaaagatgggggcaagncttttatacnatngnngnacnccntnatcaantggngcnagntntgnccnagcgnncccnntccncnnantccccngnnttctttcccnnnnttnnanntntctnagnaaaatgnngnngnncccccannnnnnnnnnnannnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn'
rv='nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnttnnnnnnnnnnnnnttggcctttaaaaaanaagaatttccgnccggncttnngganagctnccttttaatggnagggtattgaat-ttcccccattctaactttttngatntcgcntgtatttataaanctatatattgaaatcagtaactagaacaatggactttatnaatttaccagatttctactggttatagttgtacttgcnccatttataaaaaccccaaaattatggaattaattttgccgataactggatccatggataccccgtatgcaaattnctacagcaagtgagacatttttaacttcaacattatgtttatattatccgaatgaggcagctactgaaattgcagacagtaaatggacagagacattgtcacagttgtttttaacgaaaggatggccgacaggttcagtttactttaaaggttatgcagatattgcatcattttctgtagaaccgcagttatactgtgactataacattgtattaatgaaatatgatgtaagcttgcaattagatatgtctgaattggctgatctaatattaaatgaatggttatgcaatccaatggatataacgctatattattatcaacaaactgatgaggcgaacaaatggatatctatgggttcttcatgtacaattaaagtatgtcccctaaatacacaaacccttggaataggatgttcaaccacagacactaactcatttgaaatggtggctaatgcagagaagttagttataacagatgttgtcgatggagtcaatcacaaactggacgtaacgacaaacacgtgtacaatacgaaattgtaaaaaacttggaccaagggaaaacgttgctgtaattcaggtaggaggtccagatgtacttgatataactgcagatccaactactgcgccacagactgaaagaatgatgcgtataaactggaaaagatggtggcaagtcttttatacaatagtagactacgttaatcaaattgtgcaagttatgtccaagcgatcacgttctctagattccgctgcttctatacc----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------gagtagantt'
i=0
while i < (len(align)-20):
    kmer=align[i:i+20]
    if kmer.count('*') > 6:
        pos_1=Forward.index(fw[i + kmer.index('*'):i + 25].upper())
        print(fw[i + kmer.index('*'):i + 25].upper())
        pos_2=Reverse.index(rv[i + kmer.index('*'):i + 25].upper())
        print(rv[i + kmer.index('*'):i + 25].upper())
        break
    else:
        i=i+1
print(pos_1, pos_2)