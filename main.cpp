// ConsoleApplication8.cpp : Defines the entry point for the console application.
//

// 130 dsh.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include <stdio.h>
#include <process.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <time.h>
#include <bitset>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
//using namespace std;

#include <iomanip>
#include <sstream>

#define D 4 //number of dimensions 53
#define datasize 11 //581012
#define querysize 11 //number of queries 5???
#define K 4 //k nearest neighbors 20
#define familysize 200 //Size of hash family
#define L 40
#define Lused 25
#define Alter 10
int forcin;

#define M 15  //Number of hash functions selected
#define bucketnum 2000 //number of buckets 2000

int totalcheck[Alter] = {};
double R[Alter] = { 1,1.2,1.44,1.73,2.07,2.49,2.99,3.58,4.3,5.19 };
double data[datasize][D] = {};
double query[querysize][D] = {};
int querygroundtruth[querysize][K] = {};
int queryresult[querysize][K] = {};
int decision[datasize] = {};
//double knndistanceave[querysize] = {};

int dshrandom(int max) // generate a random num from 0 to max-1 (max may be any unsigned int)
{
	long long int longresult;
	int i = 1;
	longresult = rand();
	while (longresult < 10 * max)
	{
		i = 2 * i;
		longresult = 2 * longresult + rand();
	}
	int result = longresult%max;
	return result;
}

// generate a random variable follows Gaussian distribution where mean = 0 varaince = 1
double rand_single_gaussian()
{
	double sum = 0;
	for (int i = 0; i < 256; i++)
	{
		sum += (double)(2 * (rand() % 2) - 1) * 0.0625;
	}
	// this part is to avoid
	int test = rand() % 17;
	int sleep;
	for (int i = 0; i < test; i++) sleep = rand();
	//cout<<sum<<endl;
	//int forcin;
	//cin>>forcin;
	return sum;
}

// generate a random variable follows uniform distribution from -1 to 1
double rand_single_uniform()
{
	double sum = -1;
	double temp = 1;
	for (int i = 0; i < 20; i++)
	{
		sum += (rand() % 2) * temp;
		temp = temp / 2;
	}
	return sum;
}

//generate a simple multi gaussian distribution n dimensions variance matrix is I
void rand_multi_gaussian(double array[], int n)
{
	for (int i = 0; i < n; i++) array[i] = rand_single_gaussian();
	return;
}

//generate a n dimension random univector
void rand_multi_univector(double array[], int n)
{
	rand_multi_gaussian(array, n);
	double sum = 0;
	for (int i = 0; i < n; i++) sum += array[i] * array[i];
	sum = sqrt(sum);
	for (int i = 0; i < n; i++) array[i] = array[i] / sum;
	return;
}

int knnlist[K] = {};
double distlist[K] = {};
double bound = 0;
int tochange;

//compute the l2 square distance of two points
double distancel2sq(double id1[], double id2[], double bound)// compute the distance between id1 and id2
{
	// to be optimized: use bound to filt
	double result = 0;
	for (int i = 0; i < D; i++)
	{
		result += (id1[i] - id2[i])*(id1[i] - id2[i]);
	}
	//result = sqrt(result);
	if (result < 0)std::cout << "error negative distance" << std::endl;
	return result;
}

//compute the bound : current largest knn distance and its index: tochange
void computebound()
{
	bound = distlist[0];
	tochange = 0;
	for (int i = 1; i < K; i++)
	{
		if (distlist[i] > bound)
		{
			bound = distlist[i];
			tochange = i;
		}
	}
	return;
}

// check if data forcheck is a knn of querypoint, if yes maintain a new knn list
void addvertex(int forcheck, double querypoint[])
{
	double dist;
	dist = distancel2sq(data[forcheck], querypoint, 0); //distance calculation
	//for debug
	/*cout << "towrite " << towrite << endl;
	cout << "dist " << dist << endl;
	int forcin;
	cin >> forcin;*/
	if (knnlist[K - 1] == -1)
	{
		for (int i = 0; i < K; i++)
		{
			if (knnlist[i] == -1)
			{
				knnlist[i] = forcheck;// towrite here is a label for the datapoint
				distlist[i] = dist;
				//cout << "directly added" << endl;
				if (i == K - 1) computebound();
				return;
			}
		}
	}
	if (dist >= bound) return;
	//cout << "small than others" << endl;
	//cout << "old bound " << bound << endl;
	knnlist[tochange] = forcheck;
	distlist[tochange] = dist;
	computebound();
	//cout << "new bound " << bound << endl;
	return;
}


//indexing part
double familyvector[familysize][D + 1] = {};
//bool datahashresult[datasize][familysize][Alter] = {};
//bitset <familysize*Alter> datahashresult[datasize];
unsigned int datahashresult[datasize][L];
int hashtableindex[L][M] = {};
int datahashtable[L][datasize] = {};
int hashkeyindex[L][Alter][bucketnum] = {};
int hashkeylength[L][Alter][bucketnum] = {};

void family_generator()
{
	for (int i = 0; i < familysize; i++)
	{
		rand_multi_gaussian(familyvector[i], D + 1);
		for (int j = 0; j < D; j++)familyvector[i][j] = familyvector[i][j] / sqrt(D);
	}
	return;
}

double dotproduct(double id1[], double id2[])
{
	double result = 0;
	for (int i = 0; i < D; i++)
	{
		result += id1[i] * id2[i];
	}
	return result;
}

unsigned int getkey(int familyresult[], int tableindex[], int Rrank)
{
	unsigned int result = 0;
	for (int i = 0; i < M; i++)
	{
		result ^= familyresult[tableindex[i] * Alter + Rrank] + 0x9e3779b9 + (result << 6) + (result >> 2);
	}
	//result = result%bucketnum;
	/*if(result < 0)
	{
	cout<<"false"<<endl;
	cin>>forcin;
	}*/
	return result;
}

unsigned int getkeyindex(int familyresult[], int tableindex[])
{
	unsigned int result = 0;
	for (int i = 0; i < M; i++)
	{
		result ^= familyresult[tableindex[i]] + 0x9e3779b9 + (result << 6) + (result >> 2);
	}
	//result = result%bucketnum;
	/*if(result < 0)
	{
	cout<<"false"<<endl;
	cin>>forcin;
	}*/
	return result;
}

void pointhasher(double querypoint[], unsigned int tableresult[][Alter])
{
	double product;
	int hashresult[familysize*Alter] = {};
	for (int i = 0; i < familysize; i++)
	{
		product = dotproduct(querypoint, familyvector[i]);
		for (int j = 0; j < Alter; j++)
		{
			double temp = product;
			temp = temp / R[j];
			temp += familyvector[i][D];
			hashresult[i*Alter + j] = (int)temp;
		}
	}
	for (int i = 0; i < L; i++)
	{
		for (int j = 0; j < Alter; j++)
		{
			tableresult[i][j] = getkey(hashresult, hashtableindex[i], j);
		}
	}
	return;
}

void datasethasher()
{
	for (int k = 0; k < datasize; k++)
	{
		if (k % 100000 == 0) std::cout << "current hashing data " << k << std::endl;
		//pointhasher(data[i], datahashresult[i]);
		double product;
		int hashresult[familysize] = {};
		for (int i = 0; i < familysize; i++)
		{
			product = dotproduct(data[k], familyvector[i]);
			double temp = product;
			temp = temp / R[decision[k]];
			temp += familyvector[i][D];
			hashresult[i] = (int)temp;
		}
		for (int i = 0; i < L; i++)
		{
			datahashresult[k][i] = getkeyindex(hashresult, hashtableindex[i]);
		}
	}
	/*cout<<"finished"<<endl;
	cin>>forcin;
	while(1)
	{
	cin>>forcin;
	for(int i = 0; i < L; i++)cout<<datahashresult[forcin][i]<<endl;
	}*/
	return;
}



void familysample(int result[], int size, int needsize)
{
	std::vector<int> forchoose;
	forchoose.clear();
	for (int i = 0; i < size; i++)
		forchoose.push_back(i);
	int forswap, temp;
	for (int i = 0; i < needsize; i++)
	{
		forswap = dshrandom(size - i) + i;
		temp = forchoose[i];
		forchoose[i] = forchoose[forswap];
		forchoose[forswap] = temp;
		result[i] = forchoose[i];
	}
	forchoose.clear();
	return;
}

void generate_hashtableindex()
{
	// to be optimized in next paper
	// how to generate efficient hash table using some hashing family
	for (int i = 0; i < L; i++)
	{
		familysample(hashtableindex[i], familysize, M);
	}
	return;
}



/*vector<int> hashstorage[Alter*bucketnum];

void construct_index()
{
int hashkey;
//cout<<"test signal 4"<<endl;
for(int j = 0; j < L; j++)
{
//cout<<"test signal 5"<<endl;
for(int i = 0; i < bucketnum; i++)
{
hashstorage[i].clear();
}
cout<<"test signal 6"<<endl;
int forcin;
cin>>forcin;
for(int i = 0; i < datasize; i++)
{
hashkey = getkey(datahashresult[i],hashtableindex[j],decision[i]);
//if(i%10000 == 0){cout<<i<<"  "<<hashkey<<endl;cin>>forcin;}
//if(decision[i]*bucketnum + hashkey >= Alter*bucketnum)cout<<"test signal 3"<<endl;
hashstorage[decision[i]*bucketnum + hashkey].push_back(i);
}
//cout<<"test signal 1"<<endl;
//cin>>forcin;
int q = 0;
for(int i = 0; i < Alter * bucketnum; i++)
{
hashkeyindex[j][i/bucketnum][i%bucketnum] = q;
hashkeylength[j][i/bucketnum][i%bucketnum] = hashstorage[i].size();
for(int p = 0; p < hashstorage[i].size(); p++)
{
datahashtable[j][q++] = hashstorage[i][p];
}
hashstorage[i].clear();
}
//cout<<"test signal 2"<<endl;
//cin>>forcin;
cout<<"table "<<j<<" construted"<<endl;
}
double sumchecknum = 0;
double optimized = 0;
double avechecknum;
for(int i = 0; i < L; i++)
{
for(int k = 0; k < Alter; k++)
{
for(int j = 0; j < bucketnum; j++)
{
sumchecknum += hashkeylength[i][k][j] * hashkeylength[i][k][j] / (double) datasize;
optimized += hashkeylength[i][k][j] / (double) bucketnum;
}
}
}
avechecknum = sumchecknum / L;
cout<<"bucket balance index is "<<avechecknum<<endl;
cout<<"optimized bucket balance is "<<optimized / L<<endl;
return;
}*/

struct datakey {
	int dataid;
	int key;
};

datakey temptable[datasize];

bool datakeycompare(datakey a, datakey b)
{
	if (a.key < b.key)return true;
	if (a.key == b.key)
	{
		if (a.dataid < b.dataid)return true;
	}
	return false;
}

void construct_index()
{
	//cout<<"test signal 0"<<endl;
	//cin>>forcin;
	for (int j = 0; j < L; j++)
	{
		int hashkey;
		//cout<<"test signal 1"<<endl;
		//cin>>forcin;
		for (int i = 0; i < datasize; i++)
		{
			temptable[i].dataid = i;
			hashkey = datahashresult[i][j] % bucketnum;
			temptable[i].key = decision[i] * bucketnum + hashkey;
		}
		//cout<<"test signal 2"<<endl;
		//cin>>forcin;
		std::sort(temptable, temptable + datasize, datakeycompare);
		//cout<<"test signal 3"<<endl;
		//cin>>forcin;
		int q = 0;
		hashkeyindex[j][0][0] = 0;
		hashkeylength[j][0][0] = 0;
		for (int i = 0; i < datasize; i++)
		{
			while (q != temptable[i].key)
			{
				q++;
				/*if(q != temptable[i].key)
				{
				cout<<"test signal 101"<<endl;
				cin>>forcin;
				}*/
				/*if(q%100 == 0)
				{
				cout<<q<<" "<<i<<endl;
				cin>>forcin;
				}*/
				hashkeyindex[j][q / bucketnum][q%bucketnum] = i;
				hashkeylength[j][q / bucketnum][q%bucketnum] = 0;
			}
			hashkeylength[j][q / bucketnum][q%bucketnum]++;
			datahashtable[j][i] = temptable[i].dataid;
		}
		//cout<<"test signal 4"<<endl;
		//cin>>forcin;
		int nullsum = 0;
		int normalsum = 0;
		for (int i = 0; i < Alter * bucketnum; i++)
		{
			if (hashkeylength[j][i / bucketnum][i%bucketnum] == 0) nullsum++;
			else normalsum++;
		}
		std::cout << "table " << j << " construted" << std::endl;
		std::cout << nullsum << " " << normalsum << std::endl;
		//cout<<"test signal 5"<<endl;
		//cin>>forcin;
	}
	double sumchecknum = 0;
	double optimized = 0;
	double avechecknum;
	for (int i = 0; i < L; i++)
	{
		for (int k = 0; k < Alter; k++)
		{
			for (int j = 0; j < bucketnum; j++)
			{
				sumchecknum += hashkeylength[i][k][j] * hashkeylength[i][k][j] / (double)datasize;
				optimized += hashkeylength[i][k][j] / (double)bucketnum;
			}
		}
	}
	avechecknum = sumchecknum / L;
	std::cout << "bucket balance index is " << avechecknum << std::endl;
	std::cout << "optimized bucket balance is " << optimized / L << std::endl;
	return;
}

void diskread_double(std::string filename, double arrayA[][4]) // changed version
{
	int x;
	int colA = 0;
	int rowA = 0;
	std::string lineA;
	// Read file
    std::ifstream file;
	file.open(filename.c_str()); // Open file
	while(file.good())
    {
        while(getline(file, lineA))
        {
            std::istringstream streamA(lineA);
            colA=0;
            while(streamA >>x)
            {
                arrayA[rowA][colA] = x;
                std::cout << "This is array "<<arrayA[rowA][colA] << std::endl; //showing the content of the array.
                colA++;
                }
            rowA++;
        }
    }
	file.close();
}

void diskread_double_query(std::string filename, double arrayA[][4]) // changed version
{
	int x;
	int colA = 0;
	int rowA = 0;
	std::string lineA;
	// Read file
    std::ifstream file;
	file.open(filename.c_str()); // Open file
	while(file.good())
    {
        while(getline(file, lineA))
        {
            std::istringstream streamA(lineA);
            colA=0;
            while(streamA >>x)
            {
                arrayA[rowA][colA] = x;
                std::cout << "This is query "<<arrayA[rowA][colA] << std::endl; //showing the content of the array.
                colA++;
                }
            rowA++;
        }
    }
	file.close();
}

/*void diskread_double(std::string filename, double array[], int size) //diskread_double("covtype.data", data[0], datasize*D);
{
	FILE *fp;
	fp = fopen(filename.c_str(), "rb");
	if (fp == NULL)
	{
		std::cout << "Cannot open read file!" << std::endl;;
		exit(1);
	}
	fread(array, sizeof(double), size, fp);
	fclose(fp);

	return;
}*/

void diskwrite_double(std::string filename, double array[], int size)
{
	FILE *fp;
	/*if (fopen_s(&fp, filename.c_str(), "wb"))
	{
	cout << "Cannot open file!" << endl;;
	exit(1);
	}*/
	fp = fopen(filename.c_str(), "wb");
	if (fp == NULL)
	{
		std::cout << "Cannot open file!" << std::endl;;
		exit(1);
	}
	fwrite(array, sizeof(double), size, fp);
	fclose(fp);
	return;
}

void diskread_int(std::string filename, int array[], int size)
{
	FILE *fp;
	/*if (fopen_s(&fp, filename.c_str(), "rb"))
	{
	cout << "Cannot open file!"<< endl;;
	exit(1);
	}*/
	fp = fopen(filename.c_str(), "rb");
	if (fp == NULL)
	{
		std::cout << "Cannot open file!" << std::endl;;
		exit(1);
	}
	fread(array, sizeof(int), size, fp);
	fclose(fp);
	return;
}

void diskwrite_int(std::string filename, int array[], int size)
{
	FILE *fp;
	/*if (fopen_s(&fp, filename.c_str(), "wb"))
	{
	cout << "Cannot open file!" << endl;;
	exit(1);
	}*/
	fp = fopen(filename.c_str(), "wb");
	if (fp == NULL)
	{
		std::cout << "Cannot open write file!" << std::endl;
		exit(1);
	}
	fwrite(array, sizeof(int), size, fp);
	fclose(fp);
/*
//This is to see what's inside the results array
	std::ofstream myfile ("myResults.txt");
	if(myfile.is_open())
        myfile.close();

	std::ofstream myfile1 ("myResults.txt",std::ios_base::app | std::ios_base::out);
    if (myfile1.is_open())
      {
    myfile1 << "This is a line.\n";
    myfile1 << "This is another line.\n";
    for(int count = 0; count < size; count ++){
      myfile1 << array[count] << " " ;
    }
    myfile1 << "\n " ;
    myfile1.close();
  }
  else std::cout << "Unable to open file"; */
	return;
}


void index_module()
{
	std::cout << "indexing module begin" << std::endl;
	std::cout << "generating hashing family & tableindex" << std::endl;
	family_generator();
	generate_hashtableindex();
	std::cout << "hashing data" << std::endl;
	datasethasher();
	std::cout << "constructing index" << std::endl;
	construct_index();
	std::cout << "writing index to disk" << std::endl;
	diskwrite_double("index1.dat", familyvector[0], familysize*(D + 1));
	diskwrite_int("index2.dat", hashtableindex[0], L*M);
	diskwrite_int("index3.dat", datahashtable[0], L*datasize);
	diskwrite_int("index4.dat", hashkeyindex[0][0], L * Alter * bucketnum);
	diskwrite_int("index5.dat", hashkeylength[0][0], L * Alter * bucketnum);
	std::cout << "finished" << std::endl;
	return;
}

//querying part
int queryid[datasize];

void pointquery(double querypoint[], int result[], int id)
{
	//bitset <familysize*Alter> familyresult;
	unsigned int querytableresult[L][Alter];
	pointhasher(querypoint, querytableresult);
	for (int i = 0; i <= K; i++)
	{
		knnlist[i] = -1;
		distlist[i] = -1;
	}
	int hashkey;
	int bucketindex;
	int bucketlength;
	int tocheck;
	bound = 1000000;
	for (int n = 0; n < Alter; n++)
	{
		if (bound < 1.5*R[n])break;
		for (int i = 0; i < L; i++)
		{
			//hashkey = getkey(familyresult,hashtableindex[i],n);
			hashkey = querytableresult[i][n] % bucketnum;
			bucketindex = hashkeyindex[i][n][hashkey];
			bucketlength = hashkeylength[i][n][hashkey];
			for (int j = 0; j < bucketlength; j++)
			{
				tocheck = datahashtable[i][bucketindex + j];
				if (datahashresult[tocheck][i] != querytableresult[i][n])continue;
				if (queryid[tocheck] == id) continue;
				queryid[tocheck] = id;
				totalcheck[n]++;
				addvertex(tocheck, querypoint);
			}
		}
	}
	for (int i = 0; i < K; i++)result[i] = knnlist[i];
	return;
}

void batchquery()
{
	for (int i = 0; i < querysize; i++)
	{
		pointquery(query[i], queryresult[i], i);
		std::cout << "doing query id: " << i << " " << bound << std::endl;  // what is bound?
	}
	return;
}

// Driver function to sort the vector elements
// by second element of pairs
bool sortbysec(const std::pair<int,double> &a, const std::pair<int,double> &b)
{
	return (a.second < b.second);
}


//Computes the euclidean distance between two data points (one query and one data)
double getEuclideanDist(int queryIndex, int dataPointIndex){
  double sum=0;
  for(int i=0;i<D;i++){
    sum=(query[queryIndex][i]-data[dataPointIndex][i])*(query[queryIndex][i]-data[dataPointIndex][i]);
      }
  return sqrt(sum);
}




//run one K-nn query with euclidean distance against each data point
void runEuclidean(int queryIndex){
  //std::map<double, int> euclideanDist; //distance, index
  // declaring vector of pairs with index and euclidean distance
  std::vector< std::pair <int, double> > vect;

  for(int i=0;i<datasize;i++){
    std::cout << ".............................................." << std::endl; //added
    std::cout << "Here is the data on row "<<i<<" "<<std::endl;
    for(int j=0;j<sizeof(data[i])/8; j++)
        std::cout <<data[i][j] <<" "<< std::endl;

    std::cout << "Here is the query"<<std::endl;
    for(int j=0;j<sizeof(query[queryIndex])/8; j++)
        std::cout <<query[queryIndex][j] <<" "<< std::endl;

    std::cout << ".............................................." << std::endl; //added
    //std::cout << "One example for the query[1][1] is " << query[1][1] << std::endl;
    double dist = sqrt(distancel2sq(data[i], query[queryIndex], 100));
     std::cout << "The distance is ------------------------------------------------------------>"<< dist << std::endl;
    vect.push_back(std::make_pair(i,dist));
    //double dist = getEuclideanDist(queryIndex,i);
    //euclideanDist[dist]=i;
   // std::cout<<i<<"="<<dist<<" ";
  }
  // Using sort() function to sort by 2nd element
  // of pair
  sort(vect.begin(), vect.end(), sortbysec);
 // std::map<double, int>::iterator it=euclideanDist.begin();
  for(int i=0;i<K;i++){
     query[queryIndex][i]=vect[i].first;
    //querygroundtruth[queryIndex][i]=(it->second);
    std::cout<<vect[i].first<<" ";
  }
    std::cout<<std::endl;
}
//build the groundtruth.dat file
void buildGroundTruth(){
  std::cout << "groundtruth write to disk" << std::endl;
  for(int i=0; i<querysize; i++){
    runEuclidean(i);
  }

  diskwrite_int("groundtruth.dat", querygroundtruth[0], querysize*K);

}

void query_initialize()
{
	for (int i = 0; i < datasize; i++)
	{
		queryid[i] = -1;
	}
	return;
}



/*void query_module()
{
	diskread_double("query.dat", query[0][100]); //changed version
	std::cout << "query read from disk" << std::endl;
	std::cout << "Query size" << sizeof(query)<< std::endl; //query should be query[11][4], but why size is 352????
	query_initialize(); // initialize the query environment
	buildGroundTruth();
    //query_load(); // load the query points into memory
	diskread_double("index1.dat", familyvector[0], familysize*(D + 1));
	diskread_int("index2.dat", hashtableindex[0], L*M);
	diskread_int("index3.dat", datahashtable[0], L*datasize);
	diskread_int("index4.dat", hashkeyindex[0][0], L *Alter* bucketnum);
	diskread_int("index5.dat", hashkeylength[0][0], L *Alter* bucketnum);
	for (int k = 0; k < Alter; k++)
	{
		double sumchecknum = 0;
		double optimized = 0;
		for (int i = 0; i < L; i++)
		{
			for (int j = 0; j < bucketnum; j++)
			{
				sumchecknum += hashkeylength[i][k][j] * hashkeylength[i][k][j] / (double)datasize;
				optimized += hashkeylength[i][k][j] / (double)bucketnum;
			}
		}
		double avechecknum = sumchecknum / L;
		std::cout << k << std::endl;
		std::cout << "bucket balance index is " << avechecknum << std::endl;
		std::cout << "optimized bucket balance is " << optimized / L << std::endl;
	}
	int forcin;
	std::cin >> forcin;
	std::cout << "starting query" << std::endl;
	int start, finish;
	start = clock();
	batchquery(); // run the query
	finish = clock();
	std::cout << "query finished" << std::endl;
	std::cout << "time used: " << finish - start << std::endl;
	//int forcin;
	std::cin >> forcin;
	diskwrite_int("result.dat", queryresult[0], querysize*K);
	std::cout << "result write to disk" << std::endl;
	return;
}*/

void statistic_module()
{
	diskread_int("groundtruth.dat", querygroundtruth[0], querysize*K);
	std::cout << "groundtruth read from disk" << std::endl;
	diskread_int("result.dat", queryresult[0], querysize*K);
	std::cout << "result read from disk" << std::endl;

/*
	for(int i = 0; i < querysize; i++)
	{
	for(int j = 0; j < K; j++)
	{
	for(int p = 0; p < K; p++)
	{
	if(queryresult[i][j] == querygroundtruth[i][p]) sumrecall++;
	}
	}
	}
	double recall = (double)sumrecall/ (double)(querysize*K);*/

	std::cout<<"DSH results: "<<std::endl;
	for(int i=0;i<querysize;i++){
        for(int j=0;j<K;j++){
            std::cout<<queryresult[i][j]<<" ";
        }
        std::cout<<std::endl;
	}

	std::cout<<"Ground Truth results: "<<std::endl;
	for(int i=0;i<querysize;i++){
        for(int j=0;j<K;j++){
            std::cout<<querygroundtruth[i][j]<<" ";
        }
        std::cout<<std::endl;
	}
/*
	std::cout << "Printing data" << std::endl;
	for(int i=0;i<datasize;i++){
        for(int j=0;j<D;j++){
            std::cout<<data[i][j]<<" ";
        }
        std::cout<<std::endl;
	}

*/
	double sumerrorrate = 0;
	double aveerrorate;
	int sumrecall = 0;
	double dshdist[K], gtdist[K];
	int localrecall = 0;
	std::cout<<"Printing gtdist after: "<<std::endl;
	for (int i = 0; i < querysize; i++)
	{
		for (int j = 0; j < K; j++)
		{
			dshdist[j] = distancel2sq(data[queryresult[i][j]], query[i], 100);
			dshdist[j] = sqrt(dshdist[j]);
			gtdist[j] = distancel2sq(data[querygroundtruth[i][j]], query[i], 100);
			gtdist[j] = sqrt(gtdist[j]);
			std::cout<<gtdist[j]<<" ";
			std::cout<<dshdist[j]<<" ";
		}
		std::cout<<std::endl;


		std::sort(dshdist, dshdist + K);
		std::sort(gtdist, gtdist + K);
		for (int j = 0; j < K; j++)
		{
			if (dshdist[j] <= gtdist[K - 1] + 0.000001) { sumrecall++; localrecall++; } //0.001
			if (dshdist[j] - gtdist[j] < 0.000001) continue;                            //0.001
			double temp = (dshdist[j] - gtdist[j]) / gtdist[j];
			if (temp > 4) temp = 4;
			if (temp < -0.01) std::cout << "error perform better than optimal!" << std::endl;
			sumerrorrate += temp;
		}
		//cout<<i<<" "<<localrecall<<endl;
		//if (i % 50 == 0)
		//{
		//	std::cout << localrecall << std::endl;
		//	localrecall = 0;
		//}
	}
	double recall = (double)sumrecall / (double)(querysize*K);
	aveerrorate = sumerrorrate / (querysize*K);
	//std::cout<<"dsh running time is "<<(double)dshtime/(double)1000<<" ms"<<endl;
	//std::cout<<"naive running time is "<<(double)naivetime/(double)1000<<" ms"<<endl;
	std::cout << "recall is " << recall << std::endl;
	std::cout << "error rate is " << aveerrorate << std::endl;
	for (int i = 0; i < Alter; i++)
	{
		std::cout << i << std::endl;
		std::cout << "total check number is " << totalcheck[i] << std::endl;
		std::cout << "check rate is " << ((double)totalcheck[i]) / (1000 * datasize) << std::endl;
	}
	return;
}

int main()
{
    diskread_double("covtype.data", data);//changed version
    diskread_double_query("query.dat", query);//added function to read query
	//diskread_double("covtype.data", data[0], datasize*D);
	diskread_int("decision.dat", decision, datasize);
	//diskread_double("R.dat", R, Alter);
	for (int i = 0; i < Alter; i++)R[i] = 1.2*R[i];
	std::cout << "data read from disk" << std::endl;
	index_module();

    buildGroundTruth();
	//query_module();
	statistic_module();
	return 0;
}

