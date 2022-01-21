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


def find_outgroup_index(outgroup='', seq_names=''):
    """Function to return the index of the outgroup (because of the cut of the name sequences).

    Args:
        outgroup (str, optional): The name of the outgroup. Defaults to ''.
        seq_names (list of str, optional): A list with all sequence names in order. Defaults to ''.

    Returns:
        int: The index of the outgroup.
    """
    for c in range(len(seq_names)):
        if seq_names[c].replace('>', '').strip().lower() == outgroup.lower():
            return c


def verify_equal_names(name_list):
    """Function to edit the equal sequence names.

    Args:
        name_list (list of str): A list with all the sequence names.

    Returns:
        list of str: The sequence names edited if are equal names.
    """
    for index_name in range(len(name_list)):
        counter = 2
        for next_index_name in range(index_name+1, len(name_list)):
            if name_list[index_name] == name_list[next_index_name]:
                name_list[next_index_name] = name_list[next_index_name].strip()[:99-len(str(counter))] + str(counter)
                counter += 1
    return name_list


def write_nexus(nex_file, lines, ntax, nchar, outgroup_index, name_lines=''):
    """A function to write a nexus file via a fasta file.

    Args:
        nex_file (str): name of the nexus file.
        lines (list of str): a list with all the lines of the fasta file (to convert to the nexus format).
        ntax (int): the value of ntax.
        nchar (int): the value of nchar.
        name_lines (list of str, optional): your sequence names. Defaults to ''.
    """    
    mb_file = nex_file.rsplit('.')[0]
    with open(nex_file, 'w') as nexus:
        nexus.write(f"""#NEXUS

BEGIN DATA;
DIMENSIONS NTAX={ntax} NCHAR={nchar};
FORMAT DATATYPE=DNA MISSING=N GAP=-;
MATRIX""")

        counter = 0
        for line in lines:  # writing the lines
            if '>' in line:
                nexus.write(name_lines[counter])
                counter += 1
            elif line == "\n":
                continue
            else:
                nexus.write(line.strip())

        nexus.write(f"""
  ;
END;

begin mrbayes;
  set quitonerror=no;
  set autoclose=yes;""")
        if outgroup_index:
            nexus.write("  outgroup {name_lines[outgroup_index].strip()};")
        nexus.write(f"""
  mcmcp ngen={sys.argv[2]} printfreq=1000 samplefreq=100 diagnfreq=1000 nchains=4 savebrlens=yes filename={mb_file};
  mcmc;
  sumt filename={mb_file};
end;
""")


def print_stdout(filename):
    """Function to write to stdout the lines of a file.

    Args:
        filename (str): name of the file.
    """
    with open(filename, "r") as nexus:
        print(nexus.read())


def main():
    """
    The main function.
    """
    lines = open_fasta(sys.argv[1]) # getting the lines of the fasta file
    name_lines = list(filter(lambda lines='': '>' in lines.strip()[0:1], lines)) # getting all lines with the character ">"
    if len(sys.argv) > 3:
        outgroup_index = find_outgroup_index(sys.argv[3], name_lines)
    else:
        outgroup_index = None
    name_lines = list(map(lambda line='': line.replace('>', '').strip()[:99], name_lines)) # cutting all sequence names to 99 characters
    max_len_lines = max_length(name_lines)
    if max_len_lines >= 99:
        error_file = sys.argv[1].rsplit('.', 1)[0] + "_messages.txt" # telling the user that the sequence names have more than 99 characters
        handle_errors(error_file, "WARNING: some of your data names had passed the character limit and they'll be limited to 99 characters.")
        max_len_lines = 99
    # formating the sequence names
    name_lines = verify_equal_names(name_lines)
    name_lines_formatted = list(map(lambda line='': '\n' + ' '*(max_len_lines-len(line) + 1) + line.strip().replace(' ', '_') + '  ', name_lines))
    #  first_seq....
    # second seq....
    #  third_seq....
    nexus = sys.argv[1].rsplit('.')[0] + '.nex' # nexus file name
    write_nexus(nexus, lines, ntax=len(name_lines), nchar=nchar(lines), outgroup_index=outgroup_index, name_lines=name_lines_formatted)
    print_stdout(nexus)


if __name__ == '__main__':
    main()
