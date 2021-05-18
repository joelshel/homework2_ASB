#!/usr/bin/env python3

"""
CONVERT

"
>IADE3
ACCCCACTATGCTAAGCCATAAATATTGATAGATAA-ATTACAATACTTTCCGCCAGAGA
ACTACAAGTGAAAAACTTGAAACTCAAAGGACTTGGCGGTGTCCCACATTCAGCCTAGAG
GAGCCTGTCCTATAATCGATACCCCACGTTTTACCTCACCATCACTAGCACTAA-CTCAG
CCTATATACCGCCGTCGA-CAGCTTACCCCATGAGGGAAAAATAGTAAGCAAAATAGCCC
TC---CCCGCTAATACGTCAGGTCAAGGTGTAGCTCATGTGACGGAAGAGATTGGCTACA
TTTTTTATATTAAAAAACACGGAATGCTACATGA--AAAATAACATGAAGGCGAATTTAG
TAGTAAGACAGACAAGAGAACCTGTCTTAATAATGCTCTGGGACGCGCACACACCGCCCG
TCACCC

>IADE3
ACCCCACTATGCTAAGCCATAAATATTGATAGATAA-ATTACAATACTTTCCGCCAGAGA
ACTACAAGTGAAAAACTTGAAACTCAAAGGACTTGGCGGTGTCCCACATTCAGCCTAGAG
GAGCCTGTCCTATAATCGATACCCCACGTTTTACCTCACCATCACTAGCACTAA-CTCAG
CCTATATACCGCCGTCGA-CAGCTTACCCCATGAGGGAAAAATAGTAAGCAAAATAGCCC
TC---CCCGCTAATACGTCAGGTCAAGGTGTAGCTCATGTGACGGAAGAGATTGGCTACA
TTTTTTATATTAAAAAACACGGAATGCTACATGA--AAAATAACATGAAGGCGAATTTAG
TAGTAAGACAGACAAGAGAACCTGTCTTAATAATGCTCTGGGACGCGCACACACCGCCCG
TCACCC
"

TO

"
#NEXUS

BEGIN DATA;
DIMENSIONS NTAX=38 NCHAR=426;
FORMAT DATATYPE=DNA MISSING=N GAP=-;
MATRIX
    IADE3  ACCCCACTATGCTAAGCCATAAATATTGATAGATAA-ATTACAATACTTTCCGCCAGAGAACTACAAGTGAAAAACTTGAAACTCAAAGGACTTGGCGGTGTCCCACATTCAGCCTAGAGGAGCCTGTCCTATAATCGATACCCCACGTTTTACCTCACCATCACTAGCACTAA-CTCAGCCTATATACCGCCGTCGA-CAGCTTACCCCATGAGGGAAAAATAGTAAGCAAAATAGCCCTC---CCCGCTAATACGTCAGGTCAAGGTGTAGCTCATGTGACGGAAGAGATTGGCTACATTTTTTATATTAAAAAACACGGAATGCTACATGA--AAAATAACATGAAGGCGAATTTAGTAGTAAGACAGACAAGAGAACCTGTCTTAATAATGCTCTGGGACGCGCACACACCGCCCGTCACCC
  ;
END;

begin mrbayes;
  set autoclose=yes;
  outgroup Podarcis;
  mcmcp ngen=200000 printfreq=1000 samplefreq=100 diagnfreq=1000 nchains=4 savebrlens=yes filename=MyRun01;
  mcmc;
  sumt filename=MyRun01;
end;
"

"""




import sys

def open_fasta(fasta_file=''):
    try:
        with open(fasta_file, 'r') as fasta:
            return fasta.readlines()
    except FileNotFoundError:
        pass


# quantas condições????
# > até 99 caracteres
# > com mais de 99 caracteres
# sem >


def max_length(lines=''):
    max_len = 0
    for line in lines:
        if max_len < len(line):
            max_len = len(line)
    return max_len


def nchar(lines=''):
    count = 0
    length = 0
    for line in lines:
        if ">" in line:
            count += 1
        elif count == 2:
            break
        else:
            length += len(line.strip())
    return length



def write_nexus(nex_file='', lines='', lines_formated=False, name_lines_formated='', ntax = 0, nchar=0):    
    with open(nex_file, 'w') as nexus:
        nexus.write(f"""#NEXUS

BEGIN DATA;
DIMENSIONS NTAX={ntax} NCHAR={nchar};
FORMAT DATATYPE=DNA MISSING=N GAP=-;
MATRIX""")

        if lines_formated:
            counter = 0
        for line in lines:
            if '>' in line:
                if lines_formated:
                    nexus.write('\n' + name_lines_formated[counter] + '  ')
                    counter += 1
                else:
                    nexus.write('\n' + line[:99].strip().replace(' ', '_').replace('>', ' ') + '  ')
            elif line == "\n":
                pass
            else:
                nexus.write(line.strip())

        nexus.write(f"""
  ;
END;

begin mrbayes;
  set autoclose=yes;
  outgroup {sys.argv[2]};
  mcmcp ngen={sys.argv[3]} printfreq=1000 samplefreq=100 diagnfreq=1000 nchains=4 savebrlens=yes filename=MyRun01;
  mcmc;
  sumt filename=MyRun01;
end;
""")


if __name__ == '__main__':
    lines = open_fasta(sys.argv[1])
    name_lines = tuple(filter(lambda lines='': '>' in lines.strip()[0:1], lines))
    max_len = max_length(name_lines)
    max_len = max_len if max_len < 100 else 99
    name_lines_formated = tuple(map(lambda line='': ' '*(max_len-len(line)) + line.strip()[:100].replace(' ', '_').replace('>', ' '), name_lines))
    #  first_seq....
    # second seq....
    #  third_seq....
    make_nexus = lambda fasta_file='': fasta_file.split('.')[0] + '.nex'
    nexus = make_nexus(sys.argv[1])
    nchar = nchar(lines)
    write_nexus(nexus, lines, True, name_lines_formated, ntax=len(name_lines), nchar=nchar)

