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
#include <sstream>

#define D 53
#define datasize 581012
#define familysize 80
#define L 40
#define Alter 10
#define M 15
#define bucketnum 9973
#define thresholdpoint 500
#define thresholdtable 15
using namespace std;

int forcin;
int totalcheck[Alter] = {};
double R[Alter] = {1,1.2,1.44,1.73,2.07,2.49,2.99,3.58,4.3,5.19};
double data[datasize][D] = {};
int decision[datasize] = {};
int decisionsignal[datasize] = {};
//double knndistanceave[querysize] = {};

int dshrandom(int max) // generate a random num from 0 to max-1 (max may be any unsigned int)
{
	long long int longresult;
	int i = 1;
	longresult = rand();
	while (longresult < 10*max)
	{
		i = 2 * i;
		longresult = 2*longresult + rand();
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
	int test = rand()%17;
	int sleep;
	for(int i = 0; i < test; i++) sleep = rand();
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



//indexing part
double familyvector[familysize][D+1] ={};
unsigned int datahashresult[datasize][L];
float dataproduct[datasize][familysize] = {};
int hashtableindex[L][M] = {};
int hashkeylength[L][bucketnum] ={};
int sumcount[datasize] = {};

void family_generator()
{
    for(int i = 0; i < familysize; i++)
    {
        rand_multi_gaussian(familyvector[i], D+1);
        for (int j = 0; j < D; j++)familyvector[i][j] =  familyvector[i][j]/sqrt(D);
    }
    return;
}

double dotproduct(double id1[], double id2[])
{
    double result = 0;
	for (int i = 0; i < D; i++)
	{
		result += id1[i]*id2[i];
	}
	return result;
}

unsigned int getkeyindex(int familyresult[], int tableindex[])
{
    unsigned int result = 0;
    for(int i = 0; i < M; i++)
    {
            result ^= familyresult[tableindex[i]] + 0x9e3779b9 + (result << 6) + (result >> 2);
    }
    //result = result%bucketnum;
    /*if(result < 0)
    {
        cout<<"false"<<endl;
        cin>>forcin;
    }*/
    return result%bucketnum;
}

void productcomputer()
{
    for(int i = 0; i < datasize; i++)
    {
        if(i%100000 == 0) cout<<"current hashing data "<<i<<endl;
        for(int j = 0; j < familysize; j++)
        {
            dataproduct[i][j] = dotproduct(data[i],familyvector[j]);
        }
    }
    return;
}

void datasethasher(int Rrank)
{
     for(int i = 0; i < L; i++)for(int j = 0; j < bucketnum; j++)hashkeylength[i][j] = 0;
     for(int i = 0; i < datasize; i++)sumcount[i] = 0;
     for(int k = 0; k < datasize; k++)
     {
            if(k%100000 == 0) cout<<"current hashing data "<<k<<endl;
            //pointhasher(data[i], datahashresult[i]);
            //double product;
            int hashresult[familysize] = {};
            for(int i = 0; i < familysize; i++)
            {
                //product =  dotproduct(data[k],familyvector[i]);
                double temp = dataproduct[k][i];
                temp =  temp/R[Rrank];
                temp += familyvector[i][D];
                hashresult[i] = (int)temp;
            }
            for(int i = 0; i < L; i++)
            {
                datahashresult[k][i] = getkeyindex(hashresult,hashtableindex[i]);
                hashkeylength[i][datahashresult[k][i]]++;
            }
     }
     int sum = 0;
     for(int i = 0; i < datasize; i++)
     {
         for(int j = 0; j < L; j++)
         {
             if(hashkeylength[j][datahashresult[i][j]] >= thresholdpoint)sumcount[i]++;
         }
         if(sumcount[i] >= thresholdtable)
         {
             if(decisionsignal[i] == 0)
             {
                 sum++;
                 decisionsignal[i] = 1;
                 decision[i] = Rrank;
             }
         }
     }
     cout<<sum<<" points qualified for Round"<<Rrank<<endl;
     return;
}

void familysample(int result[], int size, int needsize)
{
     vector<int> forchoose;
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
     for(int i = 0; i < L; i++)
     {
             familysample(hashtableindex[i], familysize, M);
     }
     return;
}


void diskread_double(string filename, double array[], int size)
{
    FILE *fp;
	fp = fopen(filename.c_str(),"rb");
	if(fp == NULL)
	{
	    cout << "Cannot open read file!" << endl;;
		exit(1);
	}
	fread(array, sizeof(double), size, fp);
	fclose(fp);
    return;
}

void diskwrite_double(string filename, double array[], int size)
{
    FILE *fp;
	/*if (fopen_s(&fp, filename.c_str(), "wb"))
	{
		cout << "Cannot open file!" << endl;;
		exit(1);
	}*/
	fp = fopen(filename.c_str(),"wb");
	if(fp == NULL)
	{
	    cout << "Cannot open file!" << endl;;
		exit(1);
	}
	fwrite(array, sizeof(double), size, fp);
	fclose(fp);
    return;
}

void diskread_int(string filename, int array[], int size)
{
    FILE *fp;
	/*if (fopen_s(&fp, filename.c_str(), "rb"))
	{
		cout << "Cannot open file!"<< endl;;
		exit(1);
	}*/
	fp = fopen(filename.c_str(),"rb");
	if(fp == NULL)
	{
	    cout << "Cannot open file!" << endl;;
		exit(1);
	}
	fread(array, sizeof(int), size, fp);
	fclose(fp);
    return;
}

void diskwrite_int(string filename, int array[], int size)
{
	FILE *fp;
	fp = fopen(filename.c_str(),"wb");
	if(fp == NULL)
	{
	    cout << "Cannot open write file!" << endl;
		exit(1);
	}
	fwrite(array, sizeof(int), size, fp);
	fclose(fp);
    return;
}


void index_module()
{
	cout<<"indexing module begin"<<endl;
	cout<<"generating hashing family & tableindex"<<endl;
	family_generator();
	generate_hashtableindex();
	productcomputer();
	for(int i = 0; i < datasize; i++)decisionsignal[i] = 0;
	for(int i = 0; i < Alter - 1; i++)
	{
    cout<<"Round: "<<i<<endl;
	cout<<"hashing data"<<endl;
	datasethasher(i);
	}
	for(int i = 0; i < datasize; i++)
	{
	    if(decisionsignal[i] == 0)
             {
                 decisionsignal[i] = 1;
                 decision[i] = Alter - 1;
             }
	}
	cout<<"finished"<<endl;
	int decisioncount[Alter] = {};
	for(int i = 0; i < 10; i++)
	{
	    for(int j = 0; j < Alter; j++)decisioncount[j] = 0;
	    for(int j = 0; j < datasize/10; j++)
	    {
	        decisioncount[decision[i*datasize/10+j]]++;
	    }
	    cout<<"Case "<<i<<": "<<endl;
	    for(int j = 0; j < Alter; j++)cout<<decisioncount[j]<<" ";
	    cout<<endl;
	}
	cin>>forcin;
	return;
}

void read_csv(string datapath){
   //float data[datasize][D];
    std::ifstream file(datapath.c_str());

    for(int row = 0; row < datasize; ++row)
    {
        std::string line;
        std::getline(file, line);
        if ( !file.good() )
            break;

        std::stringstream iss(line);

        for (int col = 0; col < D; ++col)
        {
            std::string val;
            std::getline(iss, val, ',');
            if ( !iss.good() )
                break;

            std::stringstream convertor(val);
            convertor >> data[row][col];
        }
    }

}


int main()
{
    string datapath ="C:\\covtype.data";
    read_csv(datapath);
    //diskread_double(datapath, data[0], datasize*D);
    diskwrite_double("covtype.data", data[0], datasize*D);
    //diskread_int("decision.dat",decision, datasize);
    //diskread_double("R.dat", R, Alter);
    for(int i = 0; i < Alter; i++)R[i] = 1.2*R[i];
    cout<<"data read from disk"<<endl;
    index_module();
    diskwrite_int("decision.dat",decision, datasize);
    //query_module();
    //statics_module();
	return 0;
}

