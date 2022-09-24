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

    // Graph[0]->printArray();
    // fill the infected list flag array
    int vertices =  maxRow;
    quarantinedFlag.resize(vertices);

    for (int k = 0; k < infectedList.size(); k++)
        quarantinedFlag[infectedList[k]] = 1;

    
     const int* rowPtr;
    const int* colPtr;
    const double* values;

    double startTime  = omp_get_wtime();
    double tTotal     = 0;
    
    for (int day = 0; day < N_inputFiles; day++)
    {

        colPtr = Graph[day]->ColPtr.data();
        rowPtr = Graph[day]->RowPtr.data();
        values = Graph[day]->values.data();

        int StartPos = 0;
        int Max_contact         = 0;
        int Max_contact_vertex  = 0;


        double t1 = omp_get_wtime(); 

        //Pragma omp shared Variables
        std::vector<int> o_N_infected(N_THREADS,0);
        std::vector<int> o_N_quarantined(N_THREADS,0);
        std::vector<int> o_N_safe(N_THREADS,0);
        std::vector<int>  o_N_MaxContact(N_THREADS,0);
        std::vector<int>  o_N_MaxContactVert(N_THREADS,0);




        for (int stage = 0; stage < N_STAGE_QUARANTINE ; stage++)
        {
            std::cout << "Stage : " << stage << std::endl;
            int N_infected = infectedList.size();
            #pragma omp parallel shared(quarantinedFlag,StartPos,infectedList,rowPtr, colPtr)
            {
                int o_totalSearch       = N_infected - StartPos;
                int o_searchPerThread   = o_totalSearch/N_THREADS;
                int threadId = omp_get_thread_num();
                int o_start = threadId*o_searchPerThread;
                o_start = StartPos + o_start;
                int o_end   = (threadId+1)*o_searchPerThread;
                if(threadId == N_THREADS-1) o_end = N_infected;

                for (int i = o_start; i < o_end; i++)
                {
                    int searchVertex =  infectedList[i];
                    int newSize = 0;

                    int begin = rowPtr[searchVertex];
                    int end   = rowPtr[searchVertex+1];
                    int size = end - begin;
                    // std::cout << "size : " << size <<" maxCon : " <<Max_contact <<" vert : "<< searchVertex << std::endl;
                   
                    // if(size > o_N_MaxContact[threadId] )
                    // {
                    //     o_N_MaxContact[threadId] = size;
                    //     o_N_MaxContactVert[threadId] = searchVertex;
                    // }

                    for (int k = begin; k < end; k++)
                    {
                        int vertNo  = colPtr[k];
                        // std::cout << "vert : " <<  searchVertex << "  vertno["<<k<<"] : " << vertNo 
                        //             << " Flag : " << quarantinedFlag[vertNo] << std::endl;
                        if(!abs(quarantinedFlag[vertNo]))
                        {   
                            #pragma omp critical
                            {
                                infectedList.emplace_back(vertNo);
                            }
                            quarantinedFlag[vertNo]++;
                            double randNo = double (rand())/RAND_MAX;
                            if(randNo > PROBABILITY_INFECTION)
                                quarantinedFlag[vertNo]  = -99;
                        }   
                    }
                }
            }
            StartPos = N_infected; 

            // int distance = std::distance(o_N_MaxContact.data(), std::max_element(o_N_MaxContact.data(),o_N_MaxContact.data()+N_THREADS))  ;
            // Max_contact         = o_N_MaxContact[distance];
            // Max_contact_vertex  = o_N_MaxContactVert[distance];

        }

        int N_inf = 0 ;
        int N_quar = 0;
        int N_safe = 0;

        #pragma omp parallel shared(quarantinedFlag)
        {
            int threadId        =   omp_get_thread_num();
            int o_total         =   quarantinedFlag.size();
            int o_vertPerThread =   o_total/N_THREADS;
            int o_start         =   o_vertPerThread * threadId;
            int o_end           =   o_vertPerThread * (threadId + 1);
            if(threadId == N_THREADS-1) o_end = o_total;
            for (int i = o_start; i < o_end; i++)
            {
                if(quarantinedFlag[i] == -99) o_N_infected[threadId]++;
                else if  (quarantinedFlag[i] ==1) o_N_quarantined[threadId]++;
                else o_N_safe[threadId]++;
            }
        }
        
        N_inf   =   std::accumulate(o_N_infected.begin(),o_N_infected.end(),0);
        N_quar  =   std::accumulate(o_N_quarantined.begin(),o_N_quarantined.end(),0);
        N_safe  =   std::accumulate(o_N_safe.begin(),o_N_safe.end(),0);


        double t2 = omp_get_wtime();
        tTotal += t2 - t1;

        int total = quarantinedFlag.size();
        std::cout  << std::endl;
        std::cout << "*************** DAY :  " << N_inputFiles -  day << " *************** " <<std::endl;
        std::cout << "TOTAL PEOPLE  :   "  << vertices << std::endl;
        std::cout << "Max Contact Vertex: " <<Max_contact_vertex << "  Contact Count : " << Max_contact << std::endl;
        std::cout << "---------- QUARANTINED INFORMATION ------------ " << std::endl;
        std::cout <<std::setw(28)<< "INFECTED IN QUARANTINE : " <<std::setw(8) << N_inf <<std::setw(35)<< "Percentage : " << double(N_inf)/total* 100 << " %"<< std::endl;
        std::cout <<std::setw(28)<< "NOT INFECTED IN QUARANTINE : " <<std::setw(8) << N_quar<<std::setw(35)<< "Percentage : " << (double(N_quar)/total) * 100 << " %" << std::endl;
        std::cout <<std::setw(28)<< "TOTAL QUARANTINE : "  <<std::setw(8)<< N_quar + N_inf<<std::setw(35)<< "Total INfected Percentage : " << double((N_quar + N_inf))/total * 100 << " %"    << std::endl;
        std::cout <<std::setw(28)<< "TOTAL SAFE : " <<std::setw(8)<< N_safe  <<std::setw(35)<< "Total Safe Percentage : " << double(N_safe)/total * 100 << " %" << std::endl;
        std::cout <<std::setw(28)<< "TIME TAKEN: " <<std::setw(8)<< t2 - t1 <<" s" << std::endl;

   }

    double endTime  = omp_get_wtime();

    std::cout << " EXECUTION TIME : " << tTotal<< std::endl;





    return 0;
}
