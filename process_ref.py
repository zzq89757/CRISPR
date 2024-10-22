from pysam import FastaFile


def reverse_complement(seq: str) -> str:
    trantab = str.maketrans("ACGTNacgtnRYMKrymkVBHDvbhd", "TGCANtgcanYRKMyrkmBVDHbvdh")
    return seq.translate(trantab)[::-1]

ref = FastaFile("/home/wayne/DataBase/Ref/Homo/hg38.fa")
output = open("./sgRNA.tsv",'w')
print("no\tchr\tstart\tend\tstrand\tsgRNA",file=output)
# print(dir(ref))
# print(ref.get_reference_length('chr1'))
# print(ref.fetch('chr1',1,10))
# ref name
# print(ref.references) 

# chr1_str = ref.fetch('chr1')

# print(len(chr1_str))
# exit()
no = 1

for chr in ref.references:
    chr_seq = ref.fetch(chr)
    chr_len = len(chr_seq)
    for i in range(chr_len):
        if chr_seq[i] == 'N':continue
        # NGG or NAG forward
        if chr_seq[i].upper() == 'G' and (chr_seq[i - 1].upper() == 'A' or chr_seq[i - 1].upper() == 'G'):
            sgRNA = chr_seq[i - 22:i + 1].upper()
            print(no,end="\t",file=output)
            print(chr,end="\t",file=output)
            print(i - 21,end="\t",file=output)
            print(i - 1,end="\t",file=output)
            print("+",end="\t",file=output)
            print(sgRNA,end="\n",file=output)
            # print(chr1_str[i - 2 :i + 1],end="\t")
            # print(chr1_str[i - 22 :i + 1],end="\n")
            no += 1

        # NGG or NAG reverse
        if chr_seq[i].upper() == 'C' and (chr_seq[i - 1].upper() == 'T' or chr_seq[i - 1].upper() == 'C'):
            sgRNA = reverse_complement(chr_seq[i:i + 23].upper())
            print(no,end="\t",file=output)
            print(chr,end="\t",file=output)
            print(i + 4,end="\t",file=output)
            print(i + 23,end="\t",file=output)
            print("-",end="\t",file=output)
            print(sgRNA,end="\n",file=output)
            # print(chr_seq[i + 3:i + 23])
            # print(chr1_str[i - 2 :i + 1],end="\t")
            # print(chr1_str[i - 22 :i + 1],end="\n")
            no += 1
            # exit()