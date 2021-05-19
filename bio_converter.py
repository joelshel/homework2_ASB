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


def handle_errors(filename, text=''):
    """A function to handle some errors and warnings to a file and to the stderr.

    Args:
        filename (str): filename you want to output the text.
        text (str, optional): text you want to print in the stderr and save in the file. Defaults to ''.
    """
    print(text, file = sys.stderr)
    sys.stderr = open(filename, 'w')
    print(text, file = sys.stderr)
    sys.stderr = sys.__stderr__ # setting the stderr configurations to the original

def open_fasta(fasta_file):
    """A function that returns the read lines of a fasta file.

    Args:
        fasta_file (str): the name of your fasta file.

    Returns:
        list of str: If successfully executed returns a list with all the lines of the fasta file.
        Otherwise prints a error mensage and creates a error file.
    """
    try:
        with open(fasta_file, 'r') as fasta:
            return fasta.readlines()
    except FileNotFoundError:
        error_file = fasta_file.rsplit('.', 1)[0] + "_errors.txt" # creating the error filename
        handle_errors(error_file, "We cannot find the file. Verify if you have typed the file name correctly.")
        exit()

# how many conditions????
# > until 99 characters
# > more than 99 characters
# without >


def max_length(lines):
    """Function that returns the greatest length of a list of lines/strings.

    Args:
        lines (list of str): A list of lines/strings.

    Returns:
        int: The greatest length of the element of the list 'lines'.
    """
    max_len = 0
    for line in lines:
        if max_len < len(line):
            max_len = len(line)
    return max_len


def nchar(lines):
    """Function to calculate the length of each sequence of a fasta file to put that value in the variable
    NCHAR of the correspondent nexus file.

    Args:
        lines (list of str): A list with all the lines of the fasta file.

    Returns:
        int: The length of each sequence.
    """
    count = 0
    length = 0
    for line in lines:
        if ">" in line:
            count += 1
        elif count == 2: # between the sequence names (or between 2 > symbols)
            break
        else:
            length += len(line.strip()) # counting the sequence length
    return length



def write_nexus(nex_file, lines, ntax, nchar, lines_formatted=False, name_lines_formatted=''):
    """A function to write a nexus file via a fasta file.

    Args:
        nex_file (str): name of the nexus file.
        lines (list of str): a list with all the lines of the fasta file (to convert to the nexus format).
        ntax (int): the value of ntax.
        nchar (int): the value of nchar.
        lines_formatted (bool, optional): If true the name lines will be formatted according to name_lines_formatted. Defaults to False.
        name_lines_formatted (list of str, optional): your name sequences with your formatation. Defaults to ''.
    """    
    with open(nex_file, 'w') as nexus:
        nexus.write(f"""#NEXUS

BEGIN DATA;
DIMENSIONS NTAX={ntax} NCHAR={nchar};
FORMAT DATATYPE=DNA MISSING=N GAP=-;
MATRIX""")

        if lines_formatted:
            counter = 0
        for line in lines:
            if '>' in line:
                if lines_formatted:
                    nexus.write('\n' + name_lines_formatted[counter] + '  ') # finishing formating the sequence names
                    counter += 1
                else:
                    nexus.write('\n' + line[:99].strip().replace(' ', '_').replace('>', ' ') + '  ') # printing normally the sequence names
            elif line == "\n":
                continue
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
    lines = open_fasta(sys.argv[1]) # getting the lines of the fasta file
    name_lines = tuple(filter(lambda lines='': '>' in lines.strip()[0:1], lines)) # getting all lines with the character ">"
    max_len_lines = max_length(name_lines)
    if max_len_lines > 99:
        error_file = sys.argv[1].rsplit('.', 1)[0] + "_errors.txt" # telling the user that the sequences names have more than 99 characters
        handle_errors(error_file, "WARNING: some of your data names had passed the character limit and they'll be limited to 99 characters.")
    max_len_lines = max_len_lines if max_len_lines < 100 else 99 # variable that limits the length of the sequence names
    # formating the name of the sequences
    name_lines_formatted = tuple(map(lambda line='': ' '*(max_len_lines-len(line)) + line.strip()[:100].replace(' ', '_').replace('>', ' '), name_lines))
    #  first_seq....
    # second seq....
    #  third_seq....
    nexus = sys.argv[1].rsplit('.')[0] + '.nex' # name of the nexus file
    write_nexus(nexus, lines, ntax=len(name_lines), nchar=nchar(lines), lines_formatted=True, name_lines_formatted=name_lines_formatted)
