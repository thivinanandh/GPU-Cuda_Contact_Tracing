#include<iostream>
#include<cstring>
#include<fstream>
#include<vector>
#include<sstream>
#include<cstdlib>
#include<algorithm>
#include<omp.h>
#include<random>
#include<ctime>
#include <iomanip>      // std::setw
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>


static int N_THREADS;

int static fileShifter(int startLine, const char* filename, int threadNo)
{
    std::string file(filename);
    std::string command ;
    std::string fileout;
    std::string fileoutrm;
    fileout = "b_" + std::to_string(threadNo);
    fileoutrm = "rm " + fileout;
    command  = "head -n " + std::to_string(startLine) + " " + file + " | wc -c >> " + fileout;
    // std::cout << command << std::endl;
    
    system(command.c_str());
    int position;
    std::ifstream file1(fileout.c_str());
    file1 >> position;
    file1.close();
     system(fileoutrm.c_str());
    return position;
}

#include"sparse_parallel.h"

using namespace std;


__global__ void GPU_contactTracer(int* rowPtrDevice,int* colPtrDevice,int* infectedListDevice,int* quarantinedFlagDevice,
                                 int rowSize, int N_threadCycle , int infectedListSize , int quarantinedFlagSize , int day)
{
    const unsigned int ThreadId = blockIdx.x*blockDim.x + threadIdx.x;

    // if(!ThreadId)
    // {
    //     printf(" day : %d\n", day);

    //     for(int k = 0 ; k < infectedListSize ; k++)
    //     {
    //         printf("%d ", infectedListDevice[k]);
    //     }
    // }

    for (int ThreadCycle = 0 ; ThreadCycle < N_threadCycle  ; ThreadCycle++)
    {
        unsigned long int ThreadId_new = ThreadId + ThreadCycle*blockDim.x;             

        if(ThreadId_new < infectedListSize)
        {
            int searchVertex    = infectedListDevice[ThreadId_new];
            int begin           = rowPtrDevice[searchVertex];
            int end             = rowPtrDevice[searchVertex+1];
            int size            = end - begin;

            for (int k = begin; k < end; k++)
            {
                int vertNo = colPtrDevice[k];
                quarantinedFlagDevice[vertNo] = day+1;
            }
        }
    }

}


int main(int argc , char** argv)
{
    if(argc < 3)
    {
        cout << " ERROR : Insufficient input variables  " <<endl;
        cout << " INFO  :  arg-1 - InputFile Name "  <<endl;
        exit(0);   // EXIT Statement
    }
    N_THREADS = stoi(argv[2]);
    srand(time(NULL));

    // ---------------- MAIN FILE VARIABLES ----------------------- //
    //Variable to Store all the filenames  ---1st one - input data , remaining files
    std::vector<string> fileNames;
    std::vector<int> infectedList;          // infected People list
    std::vector<int> quarantinedList;
    std::vector<int> quarantinedFlag;
    int N_inputFiles;

    const int N_STAGE_QUARANTINE = 1;
    const double PROBABILITY_INFECTION = 0.5;
    // --------------------------- END - MAIN FILE VARIABLE --------//
    ifstream file(argv[1]);
    
    //Temp variables 
    std:string temp;
    int t_n ;
    

    // ------ Read the Filenames of the Graph from the Input File ----- //
    while(getline(file,temp))
    {
       t_n = stoi(temp);
        break;    
    }

    fileNames.resize(t_n); 
    N_inputFiles = t_n - 1;

    file.clear(); file.seekg(0);

    int t_nLines = 0;
    while(getline(file,temp))
    {
        if(t_nLines)
            fileNames[t_nLines-1] = temp;

        t_nLines++;
    }

    file.close();
    
    // ------ END Read the Filenames of the Graph from the Input File ----- //

    //------ read input infected file ------- //
    file.open(fileNames[0]);

    while(getline(file,temp))
        infectedList.push_back(stoi(temp));
    
    file.close();
    //------ END read input infected file ------- //

    // create an object for SparseMatrix based on N input files
    SparseMatrix** Graph = new SparseMatrix*[N_inputFiles];

    int maxRow = 0;
    for (int i = 0; i < N_inputFiles ; i++)
    {
        Graph[i] = new SparseMatrix;
        const char* file = fileNames[i + 1].c_str();
        Graph[i]->SparseMatrixRead1(file);

        if (Graph[i]->m_Nrow > maxRow) maxRow = Graph[i]->m_Nrow;

    }

    struct timespec start,end, TotalStart, TotalEnd;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&start);
    // Graph[1]->printArray();
    // fill the infected list flag array
    int vertices =  maxRow;
    quarantinedFlag.resize(vertices);

    std::vector<int> quarantinedFlag_sol(vertices,0);

    for (int k = 0; k < infectedList.size(); k++)
        quarantinedFlag[infectedList[k]] = 1;


    //Pragma omp shared Variables
    std::vector<int> o_N_infected(N_THREADS,0);
    std::vector<int> o_N_quarantined(N_THREADS,0);
    std::vector<int> o_N_safe(N_THREADS,0);

    
    const int* rowPtr;
    const int* colPtr;
    const double* values;

    int* quarantinedFlagDevice;
    int* infectedListDevice ; 
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&TotalStart);

    cudaMalloc((void **) &quarantinedFlagDevice, sizeof(int)*vertices);
    cudaMalloc((void **) &infectedListDevice, sizeof(int)*vertices);


    cudaMemcpy(quarantinedFlagDevice,quarantinedFlag.data(),sizeof(int)*vertices,cudaMemcpyHostToDevice);

    int InfectednewSize  = infectedList.size();
    int prevSize = 0;
    omp_set_num_threads(N_THREADS);


    double TotalTime = 0;
    
    for (int day = 0; day < N_inputFiles; day++)
    {
        double t1 = omp_get_wtime();


        colPtr = Graph[day]->ColPtr.data();
        rowPtr = Graph[day]->RowPtr.data();
        values = Graph[day]->values.data();

        int StartPos = 0;
        int Max_contact         = 0;
        int Max_contact_vertex  = 0;

        for (int stage = 0; stage < N_STAGE_QUARANTINE ; stage++)
        {
            int* rowPtrDevice ;
            int* colPtrDevice ;
            int sizeRowPtr      = Graph[day]->RowPtr.size();
            int sizeColPtr      = Graph[day]->ColPtr.size();
            int sizeInfected   = infectedList.size();

            cudaMalloc((void **) &rowPtrDevice, sizeof(int)* sizeRowPtr);
            cudaMalloc((void **) &colPtrDevice, sizeof(int)*sizeColPtr);
            cout<< "- New size  : " << InfectednewSize << " Prev size : " << prevSize <<endl;
            cudaMemcpy(rowPtrDevice,rowPtr,sizeof(int)* sizeRowPtr,cudaMemcpyHostToDevice);
            cudaMemcpy(colPtrDevice,colPtr,sizeof(int)*sizeColPtr,cudaMemcpyHostToDevice);
            cudaMemcpy(infectedListDevice + prevSize,infectedList.data() + prevSize,sizeof(int)*InfectednewSize,cudaMemcpyHostToDevice);

            prevSize = sizeInfected;
            
            int numBlocks = 1;
            
            int numThreadsPerBlock;


            if(sizeInfected < 1024) numThreadsPerBlock = sizeInfected;
            else numThreadsPerBlock = 1024;

            int N_threadCycle =  std::ceil( double(sizeInfected)/double(numThreadsPerBlock));
            dim3 dimGrid(numBlocks);
            dim3 dimBlock(numThreadsPerBlock);

            GPU_contactTracer<<<dimGrid,numThreadsPerBlock>>>(rowPtrDevice,colPtrDevice,infectedListDevice,quarantinedFlagDevice,
                                                              sizeRowPtr,N_threadCycle,sizeInfected,vertices , day);
            // cudaDeviceSynchronize();

            cudaMemcpy(quarantinedFlag_sol.data(),quarantinedFlagDevice,sizeof(int)*vertices,cudaMemcpyDeviceToHost);
            
            cudaDeviceSynchronize();
            

            InfectednewSize = 0;
            for(int j = 0 ; j < quarantinedFlag_sol.size(); j++)
            {
                if(quarantinedFlag[j] == -99 ) 
                {
                    quarantinedFlag_sol[j] = -99;
                    continue;
                }
                if((quarantinedFlag_sol[j] ==  day+1)  && !abs(quarantinedFlag[j]) )
                {

                    infectedList.push_back(j);
                    InfectednewSize += 1;
                    double randNo = double (rand())/RAND_MAX;
                    if(randNo > PROBABILITY_INFECTION)
                        quarantinedFlag_sol[j]  = -99;
               
                }

                quarantinedFlag[j] = quarantinedFlag_sol[j];
            }
            cudaFree(rowPtrDevice);
            cudaFree(colPtrDevice);
            // cudaFree(infectedListDevice);
            
        }
        

        int N_inf = 0 ;
        int N_quar = 0;
        int N_safe = 0;
        for (int i = 0; i < quarantinedFlag_sol.size(); i++)
        {
            if(quarantinedFlag_sol[i] == -99) N_inf++;
            else if  (quarantinedFlag_sol[i] >= 1) N_quar++;
            else N_safe++;
        }



        double t2 = omp_get_wtime();
        TotalTime += t2 - t1;

        
        int total = quarantinedFlag_sol.size();
        std::cout  << std::endl;
        std::cout << "*************** DAY :  " << N_inputFiles -  day << " *************** " <<std::endl;
        std::cout << "TOTAL PEOPLE  :   "  << vertices << std::endl;
        std::cout << "Max Contact Vertex: " <<Max_contact_vertex << "  Contact Count : " << Max_contact << std::endl;
        std::cout << "---------- QUARANTINED INFORMATION ------------ " << std::endl;
        std::cout <<std::setw(28)<< "INFECTED IN QUARANTINE : " <<std::setw(8) << N_inf <<std::setw(35)<< "Percentage : " << double(N_inf)/total* 100 << " %"<< std::endl;
        std::cout <<std::setw(28)<< "NOT INFECTED IN QUARANTINE : " <<std::setw(8) << N_quar<<std::setw(35)<< "Percentage : " << (double(N_quar)/total) * 100 << " %" << std::endl;
        std::cout <<std::setw(28)<< "TOTAL QUARANTINE : "  <<std::setw(8)<< N_quar + N_inf<<std::setw(35)<< "Total INfected Percentage : " << double((N_quar + N_inf))/total * 100 << " %"    << std::endl;
        std::cout <<std::setw(28)<< "TOTAL SAFE : " <<std::setw(8)<< N_safe  <<std::setw(35)<< "Total Safe Percentage : " << double(N_safe)/total * 100 << " %" << std::endl;
        std::cout <<std::setw(28)<< " ------------ TIME FOR DAY  : " <<  (end.tv_sec - start.tv_sec) + (end.tv_nsec -  start.tv_nsec)*1e-9 <<std::endl;

    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&TotalEnd);
    cout << " Time for Iteration : " << TotalTime <<endl;




    return 0;
}
