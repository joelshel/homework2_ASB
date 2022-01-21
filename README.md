# homework2_ASB
The second homework of my ASB group.  
A program to convert fasta files in nexus files. It also prints the output to **stdout**.

## How to use
This terminal program is very simple to use. All you need is call it with all its arguments. It has 3 arguments. The first argument is the fasta filename,
and the second and the third arguments are respectively the **ngen** and **outgroup** variables of the **MrBayes Block** of the nexus file. The **outgroup** variable is optional.

An example below:
```
# filename is example.fasta, ngen is 200000 and outgroup is Pocardis.
./bio_converter.py example.fasta 200000 Podarcis
```
