#ifndef _KMERGRAPH_H_
#define _KMERGRAPH_H_

#include <map>
#include <string>
#include <ext/hash_map>
#include <unistd.h>
#include "sequence.h"
#include "mympi.h"

using namespace std;

struct LKmer
{
public:
	unsigned long long seq0,seq1,seq2,seq3;    //0为高位  3为低位
	LKmer():seq3(0), seq2(0), seq1(0), seq0(0){}
    LKmer(unsigned long long a):seq3(a), seq2(0), seq1(0), seq0(0){}
    LKmer(unsigned long long a,unsigned long long b,unsigned long long c,unsigned long long d):seq3(a), seq2(b), seq1(c), seq0(d){}
    // operator unsigned long long() const{return seq3;}
    friend void operator<<=(LKmer& node, int op_num);
    friend void operator>>=(LKmer& node, int op_num);
    friend LKmer operator>>(const LKmer& node, int k);
    friend LKmer operator<<(const LKmer& node, int k);
    friend bool operator==(const LKmer& node0, const LKmer& node1);
    friend bool operator!=(LKmer& node0, const LKmer& node1);
    friend unsigned int operator&(const LKmer& node0, unsigned int op_num);
    friend int operator&(const LKmer& node0, int op_num);
    friend LKmer operator&(const LKmer& node0, const LKmer& node1);
    friend int operator|(LKmer& node0, int op_num);
    friend void operator|=(LKmer& node0, unsigned long long op_num);
    friend LKmer operator~(const LKmer& node0);
    friend bool operator>(LKmer& node0, LKmer& node1);
    friend bool operator<(LKmer& node0, LKmer& node1);

};

class kmerGraph
{
public:
	LKmer *kmers;
	unsigned char *arcs;
	unsigned long long size;
	unsigned long long read_pos;
	double commtime, commtimetot;
	double cuttime, storagetime;
	MPI_Datatype commType;

	unsigned long long readRound;
	char *readBuf;
	unsigned long long readStart;
	unsigned long long readEnd;
	unsigned long long readLen;
	
	static LKmer reverseComplement(LKmer kmerDescriptor, parameter *parameters);
    static LKmer stringToLongLong(const char *buf, int start, int end, parameter *parameters);
	static string longLongToString(LKmer a, parameter *parameters);
	static unsigned long long getProcsID(LKmer kmerID, int hashLength, MPIEnviroment *MPIcontrol);
	int arcPos(LKmer &A, LKmer &B, char directA, char directB, int hashLength);

	kmerGraph()
	{
		kmers    = NULL;
		arcs     = NULL;
		size     = 0;
		read_pos = 0;
		commtime = 0;
		commtimetot = 0;
		cuttime = 0;
		storagetime = 0;	
		
		readRound = 0;
		readBuf  = NULL;
		readStart = 0;
		readEnd   = 0;

		MPI_Type_contiguous(KMER_COMM_TYPE_LEN, MPI_UNSIGNED_LONG_LONG, &commType);
		MPI_Type_commit(&commType);
	}
	
	~kmerGraph()
	{
		if(kmers!=NULL)		delete(kmers);
		if(arcs!=NULL)		delete(arcs);
		if(readBuf!=NULL)	delete(readBuf);
		kmers = NULL;
		arcs = NULL;
		readBuf = NULL;
	}

public:
	int  constructKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol);
	void printKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol);
	void distributeKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol);
};

#endif
