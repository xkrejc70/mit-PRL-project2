#include <mpi.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>

#define TAG 69
#define NIL -69
#define ROOT 0
#define PROC0 0	// processor 0

using namespace std;

// struct that represents edge between two nodes
struct Edge {
	int id;			// edge id
	int forward;	// forward edge (1-yes, 0-no)
	bool hasNext;	// has next edge
	int next;		// next edge
	int nd1;		// node1 (from)
	int nd2;		// node2 (to)
	int rever;		// reverse edge
};

// Edge from i-th node to his aprent (i+len-2)
Edge parentEdge(int i, int len, int graphLength) {
	bool hasNext = false;
	int next = NIL;

	if (2*i + 1 < len) {
		hasNext = true;
		next = graphLength + 1;
	}
	
	return {
		i+len-2,		// id
		0,				// forward
		hasNext,		// hasNext
		next,			// next edge
		i,				// from (nd1)
		(i-1)/2,		// to (nd2)
		i-1				// reverse edge
	};
}

// Edge from i-th node to his left son (2*i)
Edge leftSonEdge(int i, int len, int graphLength) {
	bool hasNext = false;
	int next = NIL;

	if (2*i + 2 < len) {
		hasNext = true;
		next = graphLength + 1;
	}
	
	return {
		2*i,			// id
		1,				// forward
		hasNext,		// hasNext
		next,			// next edge
		i,				// from (nd1)
		2*i + 1,		// to (nd2)
		2*i - 1 + len	// reverse edge
	};
}

// Edge from i-th node to his right son (2*i+1)
Edge rightSonEdge(int i, int len, int graphLength) {
	bool hasNext = false;
	int next = NIL;

	return {
		2*i + 1,		// id
		1,				// forward
		hasNext,		// hasNext
		next,			// next edge
		i,				// from (nd1)
		2*i + 2,		// to (nd2)
		2*i + len		// reverse edge
	};
}

bool isNotRoot(int node) {
	return (node != ROOT);
}

// Find parent, left subtree and right subtree of every given node
void calcAdjList(vector<Edge> &edges, int len, int *adjList) {
	bool hasNext = false;
	int next = NIL;

	for (int i = 0; i < len; i++) {
		adjList[i] = edges.size();
		// every node has a parent (except the root)
		if (isNotRoot(i)) {
			edges.push_back(parentEdge(i, len, edges.size()));
		}
		// left subtree
		if (2*i + 1 < len) {
			edges.push_back(leftSonEdge(i, len, edges.size()));
		}
		// right subtree
		if (2*i + 2 < len) {
			edges.push_back(rightSonEdge(i, len, edges.size()));
		}
	}

	// update reverse edges ids
	int *tmp = new int [edges.size()];
	for (int i = 0; i < edges.size(); i++) {
		for (int j = 0; j < edges.size(); j++) {
			if (edges[i].rever == edges[j].id) {
				tmp[i] = j;
				break;
			}
		}
	}

	for (int i = 0; i < edges.size(); i++) {
		edges[i].rever = tmp[i];
	}
}

// MPI send data to specific process (dest) of size 1
bool mpiSend1(int data, int dest) {
	int size = 1;
	return MPI_Send(&data, size, MPI_INT, dest, TAG, MPI_COMM_WORLD);
}

// MPI send more data to specific process (dest) of any size
bool mpiSendN(int *data, int size, int dest) {
	return MPI_Send(data, size, MPI_INT, dest, TAG, MPI_COMM_WORLD);
}

// MPI recieve data from sender of size 1
bool mpiRecv1(int *buf, int sender) {
	int size = 1;
	MPI_Status status;
	return MPI_Recv(buf, size, MPI_INT, sender, TAG, MPI_COMM_WORLD, &status);
}

// MPI recieve data from sender of any size
bool mpiRecvN(int *buf, int size, int sender) {
	MPI_Status status;
	return MPI_Recv(buf, size, MPI_INT, sender, TAG, MPI_COMM_WORLD, &status);
}

// Send data to proc0, recieve size of data and then the data
void sendToProc0(int rank, int size, int val, int *buf) {
	MPI_Status status;
	// send to proc0
	mpiSend1(val, PROC0);

	int recvDataSize;
	double log = log2(rank);
	int sender = rank - pow(2, trunc(log));
	mpiRecv1(&recvDataSize, sender);
	mpiRecvN(buf, recvDataSize, sender);

	for (int i = pow(2, floor(log + 1)); i + rank < size; i *= 2) {
		mpiSend1(size, rank + i);
		mpiSendN(buf, size, rank + i);
	}
}

// Receive data from all and store them into buffer, send size of data and the data
void recvFromAll(int rank, int size, int val, int *buf) {
	MPI_Status status;
	int recvVal;
	buf[rank] = val;

	for (int i = PROC0 + 1; i < size; i++) {
		mpiRecv1(&recvVal, i);
		buf[i] = recvVal;
	}

	for (int i = PROC0 + 1; i < size; i *= 2) {
		mpiSend1(size, i);
		mpiSendN(buf, size, i);
	}
}

// Print output, convert back to chars
void printOutput(int* output, int len) {
	for (int i = 0; i < len; i++) {
		cout << (char) output[i];
	}
	cout << endl;
}


int main(int argc, char *argv[]) {
	string inputTree = argv[1];

	// input is only one char (root)
	if (inputTree.length() == 1) {
		cout << inputTree << endl;
		return 0;
	}

/************************* MPI INIT *************************/
	
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int nextProc = rank + 1;
	int predProc = rank - 1;
	int lastProc = size - 1;

	//cout << rank << endl;

/************************* ADJACENCY LIST *************************/

	vector <Edge> edges;
	int len = inputTree.length();
	int next = NIL;
	bool hasNext = false;
	int *adjList = new int[len];

	calcAdjList(edges, len, adjList);

	/*
	if (rank == 0) {
		cout << "edges:" << endl;
		for (Edge i: edges) {
			cout << i.id << ": " << i.nd1 << " " << i.nd2 << " " << i.rever << endl;
		}
	}
	*/

/************************* EULER TOUR  *************************/

	int etour;
	int reverseEdge = edges[rank].rever;
	int *etourRoute = new int [size];

	// In parallel
	// if next(e_reverse) <> nil then Etour(e) = next(e_reverse)
	if (edges[reverseEdge].hasNext) {
		etour = edges[reverseEdge].next;
	}
	
	// else Etour(e) = AdjList(v)
	else {
		etour = adjList[edges[rank].nd2]; 
	}
	
	// root node
	if (etour == ROOT) {
		etour = rank;
	}

	// send data to proc0
	if (rank != PROC0) {
		sendToProc0(rank, size, etour, etourRoute);
	}
	
	// proc0 recieve data from the rest
	if (rank == PROC0) {
		recvFromAll(rank, size, etour, etourRoute);
	}

/************************* SUFFIX SUM *************************/

	int iterationCount = trunc(log2(size)) + 1;
	int *weight = new int [size];
	int etourTmp;
	// forward: sum = 1
	// reverse: sum = 0 
	int sum = edges[rank].forward;

	for (int i = 0; i < iterationCount; i++) {
		// send and receive
		if (rank != PROC0) {
			sendToProc0(rank, size, sum, weight);
		}
		if (rank == PROC0) {
			recvFromAll(rank, size, sum, weight);
		}

		sum += weight[etourRoute[rank]];
		etourTmp = etourRoute[etourRoute[rank]];

		// send and receive
		if (rank != PROC0) {
			sendToProc0(rank, size, etourTmp, etourRoute);
		}
		if (rank == PROC0) {
			recvFromAll(rank, size, etourTmp, etourRoute);
		}
	}

	// DESC order -> sub from size
	// n - weight(e) + 1
	sum = size/2 - sum + 1;

	// send and receive
	if (rank != PROC0) {
		sendToProc0(rank, size, sum, weight);
	}

	if (rank == PROC0) {
		recvFromAll(rank, size, sum, weight);
	}

/************************* PREORDER *************************/

	int *output = new int [len];
	int recvVal;

	// send and receive
	if (rank != PROC0) {
		mpiRecv1(&recvVal, predProc);
		mpiRecvN(output, recvVal, predProc);
		// only forward edges
		if (edges[rank].forward == 1) {
			// sort input to preorder by weights
			// convert to int so it can be send as INT TYPE
			output[weight[rank]] = (int) inputTree[edges[rank].nd2];
		}
		if (rank == lastProc) {
			// last proc: send to first
			mpiSend1(len, PROC0);
			mpiSendN(output, len, PROC0);
		} else {
			// send to next
			mpiSend1(len, nextProc);
			mpiSendN(output, len, nextProc);
		}
	}

	if (rank == PROC0) {
		// sort input to preorder by weights
		output[weight[rank]] = (int) inputTree[edges[rank].nd2];
		
		// insert root
		output[ROOT] = (int) inputTree[ROOT];

		// send to next
		mpiSend1(len, nextProc);
		mpiSendN(output, len, nextProc);

		// recieve from last
		mpiRecv1(&recvVal, lastProc);
		mpiRecvN(output, recvVal, lastProc);
	}

/************************* PRINT *************************/

	if (rank == PROC0) {
		printOutput(output, len);
	}

	MPI_Finalize();
	return 0;
}