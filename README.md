# LarGo
Large Scale Genome assembler

Installation

after install these dependencies:
 * [GCC version 4.1 or larger] (http://gcc.gnu.org)
 * [MPICH version 1.4 or larger] (http://www.mpich.org)
make

Parameters

 * 'k': size of k-mer (bp), odd number less than 127 ['23']
 * 'c': cutoff threshold for edges and k-melecules. ['0']
 * 'i': the dataset file in fasta format.
 * 'o': the directory for all the output files.
 * 'h': help information for the usage of SWAP-Assembler.
 * 'v': version information of SWAP-Assembler.
 * 's': output the kmer Graph in file kmerGraph.txt.
 * 'j': output the Jung Graph in file JungGraph_arc.txt and JungGraph_mul.txt.
 * 'd': output the contig graph in file contigGraph.txt.

Example

mpirun -np 24 ./Largo -k 63 -c 5 -i ./data/S.aureus.fasta -o Saur_k31_c5
