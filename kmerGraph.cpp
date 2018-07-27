#include "sequence.h"
#include "mympi.h"
#include "kmerGraph.h"

LKmer operator>>(const LKmer& node, int k){
    LKmer tmp;
	tmp.seq3 = node.seq3 >> k%64 | node.seq2 << (64 - k%64);
	tmp.seq2 = node.seq2 >> k%64 | node.seq1 << (64 - k%64);
	tmp.seq1 = node.seq1 >> k%64 | node.seq0 << (64 - k%64);
	tmp.seq0 = node.seq0 >> k%64;
	while(k >= 64)
	{
		tmp.seq3=tmp.seq2;
		tmp.seq2=tmp.seq1;
		tmp.seq1=tmp.seq0;
		tmp.seq0=0;
		k -= 64;
	}
    return tmp;
}
LKmer operator<<(const LKmer& node, int k){ //仿佛有问题
    LKmer tmp;
	if(k%64 != 0)
	{
		tmp.seq0 = node.seq0 << k%64 | node.seq1 >> (64 - k%64);
		tmp.seq1 = node.seq1 << k%64 | node.seq2 >> (64 - k%64);
		tmp.seq2 = node.seq2 << k%64 | node.seq3 >> (64 - k%64);
		tmp.seq3 = node.seq3 << k%64 ;
	}
	else
		tmp = node;
	while(k >= 64)
	{
		tmp.seq0=tmp.seq1;
		tmp.seq1=tmp.seq2;
		tmp.seq2=tmp.seq3;
		tmp.seq3=0;
		k -= 64;
	}
    return tmp;
}

bool operator==(const LKmer& node0, const LKmer& node1){
    if (node0.seq3 == node1.seq3 && node0.seq2 == node1.seq2 && node0.seq1 == node1.seq1 && node0.seq0 == node1.seq0)
    	return true;
    else
    	return false;
}
bool operator!=(LKmer& node0, const LKmer& node1){
    if (node0.seq3 == node1.seq3 && node0.seq2 == node1.seq2 && node0.seq1 == node1.seq1 && node0.seq0 == node1.seq0)
    	return false;
    else
    	return true;
}
unsigned int operator&(const LKmer& node0, unsigned int op_num)
{
	return node0.seq3&op_num;
}
int operator&(const LKmer& node0, int op_num)
{
	return node0.seq3&op_num;
}
LKmer operator&(const LKmer& node0, const LKmer& node1)
{
	LKmer tmp;
	tmp.seq3 = node0.seq3 & node1.seq3;
	tmp.seq2 = node0.seq2 & node1.seq2;
	tmp.seq1 = node0.seq1 & node1.seq1;
	tmp.seq0 = node0.seq0 & node1.seq0;
	return tmp;
}

int operator|(LKmer& node0, int op_num)
{
	return node0.seq3|op_num;
}

LKmer operator~(const LKmer& node0)
{
	LKmer tmp;
	tmp.seq3 = ~node0.seq3;
	tmp.seq2 = ~node0.seq2;
	tmp.seq1 = ~node0.seq1;
	tmp.seq0 = ~node0.seq0;
	return tmp;
}

bool operator>(LKmer& node0, LKmer& node1)
{
	if (node0.seq0>node1.seq0)
		return true;
	else if(node0.seq0<node1.seq0)
		return false;
	else if(node0.seq1>node1.seq1)
		return true;
	else if(node0.seq1<node1.seq1)
		return false;
	else if(node0.seq2>node1.seq2)
		return true;
	else if(node0.seq2<node1.seq2)
		return false;
	else if(node0.seq3>node1.seq3)
		return true;
	else 
		return false;
}

bool operator<(LKmer& node0, LKmer& node1)
{
	if (node0.seq0<node1.seq0)
		return true;
	else if(node0.seq0>node1.seq0)
		return false;
	else if(node0.seq1<node1.seq1)
		return true;
	else if(node0.seq1>node1.seq1)
		return false;
	else if(node0.seq2<node1.seq2)
		return true;
	else if(node0.seq2>node1.seq2)
		return false;
	else if(node0.seq3<node1.seq3)
		return true;
	else 
		return false;
}

void operator<<=(LKmer& node, int k)
{
	node = node << k;
}

void operator>>=(LKmer& node, int k)
{
	node = node >> k;
}
void operator|=(LKmer& node0, unsigned long long op_num)
{
	node0.seq3 = node0.seq3 | op_num;
}

void operator&=(LKmer& node0, LKmer& node1)
{
	node0.seq3 = node0.seq3 & node1.seq3;
	node0.seq2 = node0.seq2 & node1.seq2;
	node0.seq1 = node0.seq1 & node1.seq1;
	node0.seq0 = node0.seq0 & node1.seq0;

}

LKmer kmerGraph::reverseComplement(LKmer a, parameter *parameters)
{
        LKmer rev = 0;
        for(int i=0; i<parameters->hashLength; i++){
                rev <<= 2;
                rev |= (3-a.seq3&3);
                a >>= 2;
        }
        return rev;
}

LKmer kmerGraph::stringToLongLong(const char *buf, int start, int end, parameter *parameters)
{
        LKmer ret = 0;
        for(int i=start; i<end; i++){
                ret = ret << 2;
                ret |= (unsigned long long)parameters->nucleotideValue[buf[i]]&3;
        }

        return ret;
}

string kmerGraph::longLongToString(LKmer a, parameter *parameters)
{
	string descriptor;
        descriptor.clear();
        for(int i=0;i<parameters->hashLength;i++)
        {
                descriptor += parameters->nucleotideArray[a&3];
                a = a>>2;
        }
        reverse(descriptor.begin(), descriptor.end());
	return descriptor;
}

int kmerGraph::arcPos(LKmer &A, LKmer &B, char directA, char directB, int hashLength)
{
	int lastB;
	if(directB=='+')
		lastB = B & 3;
	else
		lastB = (~(B>>((hashLength-1)*2) & 3)) & 3;   //******** by zhu.
	
	if(directA=='-')
		lastB += 4;
	return lastB;
}

unsigned long long kmerGraph::getProcsID(LKmer kmerID, int hashLength, MPIEnviroment *MPIcontrol)
{

        LKmer revKmer = 0, tmpID = kmerID;
        for(int i=0; i<hashLength; i++){
                revKmer <<= 2;
                revKmer = revKmer | (~(tmpID&3))&3;
                tmpID >>= 2;
        }

        if(revKmer>kmerID) revKmer = kmerID;
        // printf("kmer: %llu\n", kmerID.seq3);
        // printf("rev_kmer: %llu\n", revKmer.seq3); 

         unsigned int factor = 19;
         unsigned int numBytes = (hashLength + 3) / 4;

         unsigned int sum = 0;
         for(int i = 0; i < numBytes; i++){
                 sum = sum * factor + (revKmer & 0xFF);
                 revKmer >>= 8;
         }
         return sum % MPIcontrol->nprocs;
}


int kmerGraph::constructKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol)
{
	clock_t t0, t1, t2;
	t0 = clock();
	LKmer MASK = ~((LKmer)3<<(2*parameters->hashLength));

	unsigned long long readBases = MPIcontrol->nprocs*1024*16;
	if(readBases>BUF_SIZE)	readBases = BUF_SIZE;		
	readRound ++;

	if(readBuf == NULL)	
	{	
		readBuf = new char [BUF_SIZE*2+2];
		MPIcontrol->File_open(parameters->fastaPath);
    		MPIcontrol->File_locate();
	
		readLen = MPIcontrol->File_read(readBuf,BUF_SIZE);
		
	        unsigned long long i=0;
        	while(readBuf[i]!='>' && i<readLen) {i++; }
	        assert(readBuf[i]=='>');
	
		readStart = i;
		readEnd = readLen;	
		unsigned long long readDelta;
				
		MPI_Status MPIrecvstatus;
		MPI_Send(&readStart, 1, MPI_UNSIGNED_LONG_LONG, (MPIcontrol->rank-1+MPIcontrol->nprocs)%MPIcontrol->nprocs, 0, MPI_COMM_WORLD);
		MPI_Recv(&readDelta, 1, MPI_UNSIGNED_LONG_LONG, (MPIcontrol->rank+1)%MPIcontrol->nprocs, 0, MPI_COMM_WORLD, &MPIrecvstatus);
//		MPI_Sendrecv(&readStart, 1, MPI_UNSIGNED_LONG_LONG, (MPIcontrol->rank-1+MPIcontrol->nprocs)%MPIcontrol->nprocs, 0, &readDelta, 1, MPI_UNSIGNED_LONG_LONG, (MPIcontrol->rank+1)%MPIcontrol->nprocs, 0, MPI_COMM_WORLD, &MPIrecvstatus);	
	
 		MPIcontrol->datasize += readDelta; 		

		printf("Proc %d: %llu File_readLen=%llu, read_offset=%llu(%llu)\n", MPIcontrol->rank, readRound, readLen, MPIcontrol->read_offset-MPIcontrol->start_pos, MPIcontrol->datasize);
	}

	printf("Proc %d: %llu readStart=%llu, readEnd=%llu\n", MPIcontrol->rank, readRound, readStart, readEnd);	

	int endTag = 0, totalTag=0;
	if(readBases*readRound>readEnd)	endTag = 1;
	MPI_Allreduce(&endTag, &totalTag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	if( totalTag == MPIcontrol->nprocs )	
	{
		for(int i=0; i<readEnd-readStart; i++)	readBuf[i] = readBuf[readStart+i];
		readBuf[readEnd-readStart] = 0;
		readEnd = readEnd-readStart;			
		readStart = 0;
		
		readLen = MPIcontrol->File_read(readBuf+readEnd, BUF_SIZE);	
		readEnd += readLen;
		readBuf[readEnd] = 0;
		printf("Proc %d: %llu File_readLen=%llu, read_offset=%llu(%llu)\n", MPIcontrol->rank, readRound, readLen, MPIcontrol->read_offset-MPIcontrol->start_pos, MPIcontrol->datasize);

		//Check the collective I/O finished in all process.
		unsigned long long totLen = 0;	
		MPI_Allreduce(&readLen, &totLen, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);	

//		printf("Proc %d: %llu readEnd=%llu totLen=%llu, readStart=%llu, readEnd=%llu\n",MPIcontrol->rank, readRound, readEnd, totLen, readStart, readEnd);	

		readRound=1;
		if(readEnd == 0 && totLen == 0)	return 0;	//There is no reads left in the buf.	
	}

//	printf("Proc %d: %llu Processing data Block (%llu, %llu)\n", MPIcontrol->rank, readRound, MPIcontrol->read_offset-readEnd, MPIcontrol->read_offset);	
	t1 = clock();
	char *readPos = strstr(readBuf+readStart, ">");
	if(readPos == NULL)
	{
		printf("Proc %d: Can not find >\n",MPIcontrol->rank); 
		return 1;
	}
	assert(readPos == readBuf+readStart);

//	printf("Proc %d: readBuf=%s", MPIcontrol->rank, readBuf);

	unsigned long long tot_kmer_num = 0, totSize=0;
	while(readPos < readBuf + readBases*readRound)
	{
		char *readStr = strstr(readPos, "\n");
		if(readStr==NULL)	break;
		readStr++;
	
		char *readDelimiter = strstr(readStr, "\n");
		if( !(readDelimiter < readBuf + readBases*readRound) )	break;		
		if(readDelimiter==NULL)	break;

		int len = readDelimiter - readStr;
		if(len<parameters->hashLength)	
		{
			readDelimiter++;
			readPos = readDelimiter;	
			continue;
		}

		len = len - parameters->hashLength+1;
		tot_kmer_num += len;	
		readPos = readDelimiter+1;	
	}		
	
	printf("Proc %d: %llu tot_kmer_num=%d\n", MPIcontrol->rank, readRound, tot_kmer_num);


	size  = tot_kmer_num;
	kmers = new LKmer 	[size];
	arcs  = new unsigned char 	[size];

	readPos = readBuf+readStart;
	unsigned long long kmer_index=0;
	while(readPos < readBuf + readRound*readBases)
	{
		
		char *readStr = strstr(readPos, "\n");
		if(readStr==NULL)	break;
		readStr++;
	
		char *readDelimiter = strstr(readStr, "\n");
		if( !(readDelimiter < readBuf + readRound*readBases) )	break;		
		if(readDelimiter==NULL)	break;
	
		int len = readDelimiter - readStr;
		if(len<parameters->hashLength)	
		{
			readDelimiter++;
			readPos = readDelimiter;	
			continue;
		}

		for(char *p=readStr; p<readDelimiter;p++)
		{
			if(*p=='A' || *p=='C' || *p=='G' || *p=='T')	;
			else if(*p=='a' || *p=='c' || *p=='g' || *p=='t')	*p=*p+'A'-'a';
			else	*p = 'A';
		}	

	//	printf("Proc %d: readStart=%llu, readPos=%llu, readStr=%llu, readDelimiter=%llu\n", MPIcontrol->rank, readStart, readPos-readBuf, readStr-readBuf, readDelimiter-readBuf);
                
		LKmer curNodeID=0,     primeNodeID=0,   	twinNodeID=0;
		LKmer preNodeID=0,     prePrimeNodeID=0; 

                curNodeID = this->stringToLongLong(readStr, 0, parameters->hashLength-1, parameters);
                // printf("CurNode: %s\n", kmerGraph::longLongToString(curNodeID, parameters).c_str());
                for(int j=0; j<=len-parameters->hashLength; j++) {
                	curNodeID <<= 2;
                	// printf("1:%s\n", kmerGraph::longLongToString(curNodeID, parameters).c_str());
                	curNodeID |= (parameters->nucleotideValue[readStr[j+parameters->hashLength-1]]&(3ull));
                	// printf("2:%s\n", kmerGraph::longLongToString(curNodeID, parameters).c_str());
                	curNodeID &= MASK;
                	// printf("3:%s\n", kmerGraph::longLongToString(curNodeID, parameters).c_str());
                	twinNodeID  = this->reverseComplement(curNodeID, parameters);
                	primeNodeID = (curNodeID>twinNodeID) ? curNodeID : twinNodeID;

			if(j==0)	{
				preNodeID = curNodeID;
				prePrimeNodeID = primeNodeID; 
				kmers[kmer_index] = primeNodeID;	arcs[kmer_index] = 0;
				kmer_index ++;
				continue;
			}

//			kmers[kmer_index]     = prePrimeNodeID;	arcs[kmer_index] = 0;
			kmers[kmer_index]   = primeNodeID;	arcs[kmer_index]   = 0;
			
			int pos_ret;
			//preNode -> curNode  add(A, B, +, +, last(B+))
			if( preNodeID == prePrimeNodeID && curNodeID == primeNodeID)
			{
				pos_ret = arcPos(prePrimeNodeID, primeNodeID, '+', '+', parameters->hashLength);	
				arcs[kmer_index-1]    |= ((unsigned char) 1<<pos_ret);
				
				pos_ret = arcPos(primeNodeID, prePrimeNodeID, '-', '-', parameters->hashLength); 
				arcs[kmer_index] |=  ((unsigned char) 1<<pos_ret);
			}
			//preNode -> ~curNode
			else if(preNodeID == prePrimeNodeID && curNodeID != primeNodeID)
			{
				pos_ret = arcPos(prePrimeNodeID, primeNodeID, '+', '-', parameters->hashLength);	
				arcs[kmer_index-1]    |= ((unsigned char) 1<<pos_ret);
				
				pos_ret = arcPos(primeNodeID, prePrimeNodeID, '+', '-', parameters->hashLength); 
				arcs[kmer_index] |=  ((unsigned char) 1<<pos_ret);
			}
			//~preNode -> curNode
			else if(preNodeID != prePrimeNodeID && curNodeID == primeNodeID)
			{
				pos_ret = arcPos(prePrimeNodeID, primeNodeID, '-', '+', parameters->hashLength);	
				arcs[kmer_index-1]    |= ((unsigned char) 1<<pos_ret);
				
				pos_ret = arcPos(primeNodeID, prePrimeNodeID, '-', '+', parameters->hashLength); 
				arcs[kmer_index] |=  ((unsigned char) 1<<pos_ret);
			}			
			//~preNode -> ~curNode
			else if(preNodeID != prePrimeNodeID && curNodeID != primeNodeID)
			{
				pos_ret = arcPos(prePrimeNodeID, primeNodeID, '-', '-', parameters->hashLength);	
				arcs[kmer_index-1]    |= ((unsigned char) 1<<pos_ret);
				
				pos_ret = arcPos(primeNodeID, prePrimeNodeID, '+', '+', parameters->hashLength); 
				arcs[kmer_index] |=  ((unsigned char) 1<<pos_ret);
			}
			preNodeID = curNodeID;
			prePrimeNodeID = primeNodeID;
			kmer_index += 1;
		}
		readPos = readDelimiter+1;	
	}
	assert(kmer_index == size);

	readStart = readPos-readBuf;
//	printf("Proc %d: ConstructKmer Finished\n", MPIcontrol->rank);
	t2 = clock();
	cuttime += t2-t1;
	return 1;
}

void kmerGraph::printKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol)
{
}

void kmerGraph::distributeKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol)
{
	clock_t startTime, endTime;
	startTime = clock();
	LKmer *kmers_send = new LKmer [size];
	assert(kmers_send != NULL);
	unsigned char      *arcs_send  = new unsigned char      [size];
	assert(arcs_send != NULL);
	for(unsigned long long i=0; i < size; i++)	kmers_send[i] = arcs_send[i] = 0;

	int nproc = MPIcontrol->nprocs;
	
//	printf("Proc %d: AlltoAll start\n", MPIcontrol->rank);
//****	
	int  *send_size = new int [nproc];
	assert(send_size!=NULL);
	int  *recv_size = new int [nproc];
	assert(recv_size!=NULL);
	int  *send_pos  = new int [nproc];
	assert(send_pos!=NULL);
	int  *recv_pos  = new int [nproc];
	assert(recv_pos!=NULL);

	for(int i=0; i<nproc; i++)	send_size[i] = recv_size[i] = 0;

	unsigned long long send_sum=0,  recv_sum=0;
	unsigned long long ProcsID;	

	for(unsigned long long i=0; i<size; i++)
	{
		ProcsID = getProcsID(kmers[i], parameters->hashLength, MPIcontrol);
		send_size[ProcsID]++;
	}

//****	
	MPI_Alltoall(send_size,1, MPI_INT, recv_size, 1, MPI_INT, MPI_COMM_WORLD);

	for(int i=0; i<nproc; i++)	send_sum += send_size[i];
	for(int i=0; i<nproc; i++)	recv_sum += recv_size[i];
	assert(send_sum == size);
	
//****	
//	printf("Proc %d: send_sum=%llu recv_sum=%llu\n", MPIcontrol->rank, send_sum, recv_sum);

	send_pos[0] = 0;
	for(int i=1; i<nproc; i++)	send_pos[i] = send_pos[i-1] + send_size[i-1];
	recv_pos[0] = 0;
	for(int i=1; i<nproc; i++)	recv_pos[i] = recv_pos[i-1] + recv_size[i-1];

/*
        printf("proc%d: ", MPIcontrol->rank);
        for(int i=0;i<size;i++) printf("send_size[%d]=%d ", i, send_size[i]);
        printf("\t");
	printf("proc%d: ", MPIcontrol->rank);
        for(int i=0;i<size;i++) printf("recv_size[%d]=%d ", i, recv_size[i]);
        printf("\n");
*/
	for(int i=0; i<size; i++)
	{
		ProcsID = getProcsID(kmers[i], parameters->hashLength, MPIcontrol);
		kmers_send[send_pos[ProcsID]]  = kmers[i]; 
		arcs_send[send_pos[ProcsID]++] = arcs[i];	
	}

	send_pos[0] = 0;
	for(int i=1; i<nproc; i++)	send_pos[i] = send_pos[i-1] + send_size[i-1];

	delete kmers;
	delete arcs;
	kmers = new LKmer [recv_sum];
	assert(kmers!=NULL);
	arcs  = new unsigned char      [recv_sum];
	assert(arcs!=NULL);

	clock_t t1 = clock();
	MPI_Alltoallv(kmers_send,send_size,send_pos, MPIcontrol->MPI_LKmer, kmers,recv_size,recv_pos, MPIcontrol->MPI_LKmer, MPI_COMM_WORLD);	
	clock_t t2 = clock();
	MPI_Alltoallv(arcs_send, send_size,send_pos, MPI_CHAR,          arcs, recv_size,recv_pos, MPI_CHAR,          MPI_COMM_WORLD);
//	printf("Proc %d: AlltoAll Finished\n", MPIcontrol->rank);
	// printf("%llu %llu %llu %llu\n", kmers[888].seq3,kmers[0].seq2,kmers[0].seq1,kmers[0].seq0);
	// printf("%llu %llu %llu %llu\n", recv_size[0],recv_size[1],recv_size[2],recv_size[3]);
	LKmer tmmm(4611685998516522936,1020,0,0);
	for (int xx =0;xx< send_size[MPIcontrol->rank];xx++)
	{
		if (kmers_send[xx].seq3 == 18446744073709488385ull)
			printf("prod:%d, search...................",MPIcontrol->rank);
	}
	delete kmers_send;
	delete arcs_send;
	kmers_send = NULL;
	arcs_send = NULL;
	size  = recv_sum;
	
	delete send_size;
	delete recv_size;
	delete send_pos;
	delete recv_pos;
	send_size = recv_size = send_pos = recv_pos = NULL;
	
	endTime = clock();
	commtime += t2 - t1;
	commtimetot += t2-startTime;
}
/*
int main (int argc, char *argv[])
{
   	int i, size, namelen, lgsize;
	MPIEnviroment MPIcontrol;
	MPIcontrol.init(argc, argv);

	parameter parameters;
	sequence sequences;
	parameters.getParameters(argc, argv);	
	sequences.getSequences(&parameters, &MPIcontrol);

	kmerGraph mygraph;

	while(mygraph.constructKmerGraph(&sequences, &parameters, &MPIcontrol)!=0)
	{
		//mygraph.printKmerGraph(&parameters, &MPIcontrol);

	//	char message[100];
		//sprintf(message,"kmer number is %d\n", mygraph.kmers.size());	
		//MPIcontrol.print(message);

		mygraph.distributeKmerGraph(&parameters, &MPIcontrol);	

		delete mygraph.kmers;
		delete mygraph.arcs;
	
	//	sprintf(message, "distribute kmers finished\n");
	//	MPIcontrol.print(message);
	}

	MPIcontrol.finalize();
    	return (0);
}*/
