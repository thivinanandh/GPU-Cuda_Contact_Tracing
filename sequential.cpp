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

#include"sparse.h"

using namespace std;

int main(int argc , char** argv)
{
    if(argc < 2)
    {
        cout << " ERROR : Insufficient input variables  " <<endl;
        cout << " INFO  :  arg-1 - InputFile Name "  <<endl;
        exit(0);   // EXIT Statement
    }
    
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
    

    //------ read input infected file ------- //
    file.open(fileNames[0]);

    while(getline(file,temp))
        infectedList.push_back(stoi(temp));
    
    file.close();
    //------ END read input infected file ------- //


    // create an object for SparseMatrix
    SparseMatrix** Graph = new SparseMatrix*[N_inputFiles];

    for (int i = 0; i < N_inputFiles ; i++)
    {
        Graph[i] = new SparseMatrix;
        const char* file = fileNames[i + 1].c_str();
        Graph[i]->SparseMatrixRead(file);
    }

    // Graph[0]->printArray();
    // fill the infected list flag array
    int vertices =  Graph[0]->m_Nrow;
    quarantinedFlag.resize(vertices);

    // std::cout << "Resize Done " << std::endl;


    for (int k = 0; k < infectedList.size(); k++)
        quarantinedFlag[infectedList[k]] = 1;
    
   
    const int* rowPtr;
    const int* colPtr;
    const double* values;

    double t1 = omp_get_wtime();
    double tStart = t1;

    for (int day = 0; day < N_inputFiles; day++)
    {

        colPtr = Graph[day]->ColPtr.data();
        rowPtr = Graph[day]->RowPtr.data();
        values = Graph[day]->values.data();

        int StartPos = 0;
        int Max_contact = 0;
        int Max_contact_vertex = 9;

        for (int stage = 0; stage < N_STAGE_QUARANTINE ; stage++)
        {
            // std::cout << "Stage : " << stage << std::endl;
            int N_infected = infectedList.size();
            for (int i = StartPos; i < N_infected; i++)
            {
                int searchVertex =  infectedList[i];
                int newSize = 0;

                int begin = rowPtr[searchVertex];
                int end   = rowPtr[searchVertex+1];
                int size = end - begin;
                

                //size Control ( if Size > 1000 ) , restrict Size to 1000
                // if(size > 50)  {end = begin + 50;  size = 50;   }
                // std::cout << "size : " << size <<" maxCon : " <<Max_contact <<" vert : "<< searchVertex << std::endl;
                if(size > Max_contact )
                {
                    Max_contact = size;
                    Max_contact_vertex = searchVertex;
                }

                for (int k = begin; k < end; k++)
                {
                    int vertNo  = colPtr[k];
                    // std::cout << "vert : " <<  searchVertex << "  vertno["<<k<<"] : " << vertNo 
                    //             << " Flag : " << quarantinedFlag[vertNo] << std::endl;
                    if(!quarantinedFlag[vertNo])
                    {   
                        infectedList.emplace_back(vertNo);
                        quarantinedFlag[vertNo]++;
                        double randNo = double (rand())/RAND_MAX;
                        if(randNo > PROBABILITY_INFECTION)
                            quarantinedFlag[vertNo]++;
                    }   
                }
            }

            StartPos = N_infected;       
        }
        int N_inf = 0 ;
        int N_quar = 0;
        int N_safe = 0;
        for (int i = 0; i < quarantinedFlag.size(); i++)
        {
            if(quarantinedFlag[i] ==2) N_inf++;
            else if  (quarantinedFlag[i] ==1) N_quar++;
            else N_safe++;
        }
        double t2 = omp_get_wtime();
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
        std::cout << std::setw(28)<<" --- --------    TIME TAKEN PER DAY : " << t2 - t1 <<std::endl;

   }
   double tEnd = omp_get_wtime();

   std::cout << "################# FINAL TIME TAKEN : " << tEnd - tStart << " #########################" << std::endl;

    return 0;
}
