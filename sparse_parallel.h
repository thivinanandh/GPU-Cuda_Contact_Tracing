#include<iostream>
#include<cstring>
#include<fstream>
#include<vector>
#include<sstream>
#include<cstdlib>
#include<algorithm>
#include <omp.h>



using namespace std;

class SparseMatrix
{
    public:
        std::vector<int> RowPtr;
        std::vector<int> ColPtr;
        std::vector<double> values;

        //Variables
        int m_NNZ;
        int m_Nrow;
    
        SparseMatrix()
        {
            
        }


    public:
        void printArray();

        void SparseMatrixRead(const char* filename);

        void SparseMatrixRead1(const char* filename);


};


void SparseMatrix::printArray()
{
    cout << "Row Pointer "<<endl;

    for (int i = 0; i < m_Nrow + 1; i++)
        cout << RowPtr[i]<<"\t";
    cout << endl;

    cout << "Col Pointer "<<endl;
    for (int i = 0; i < m_NNZ; i++)
        cout << ColPtr[i]<<"\t";
    cout << endl;

    cout << "Values"<<endl;
    for (int i = 0; i < m_NNZ; i++)
        cout << values[i]<<"\t";
    cout << endl;

    
}




// Sequential File read 
void SparseMatrix::SparseMatrixRead(const char* filename)
{

    int N_linesFile = 0;
    int N_Vertices = 0;


    cout << " ---------------- REAAD START --------------------------------- " <<endl;
    cout << " ----- FILE NAM E: " << filename << "  -------------------------" <<endl;


    std::string command;
    std::string t_file(filename);
    command = "cat ";
    command += t_file + " |  wc -l >> a";
    cout << " cpmmand : " << command <<endl;
    system(command.c_str());
    ifstream file2;
    file2.open("a");
    int N_LINE_TOTAL = 0;
    file2 >> N_LINE_TOTAL;
    file2.close();
    system("rm a");
    // N_LINE_TOTAL++;

    // #pragma omp parallel private(file,t_val,t_max1,temp,temp1,ss,GraphRowData,N_linesFile) 
    // Graph Parameters 
    ifstream file(filename);

    //Temp Variables
    int t_prevRow = -123234345;  // Random NUmber
    int t_val[2];
    int t_max1 = -1000000000;
    int t_max2 = -100000000;
    std::string temp;
    std::string temp1;
    std::vector<std::pair<int,int> >  GraphRowData;  // ROw number and the prefix sum of all the rows

    double t1 = omp_get_wtime();

    GraphRowData.resize(2*N_LINE_TOTAL);
    cout << " NLINE TOAL : " << N_LINE_TOTAL <<endl;
    cout << " Resize done" <<endl;
    N_linesFile = 0;
   
    while(getline(file,temp))
    {
        int len = temp.length();
        size_t found = temp.find(" "); 
        int a = stoi(temp.substr(0,found));
        int b = stoi(temp.substr(found+1,len-found));
        
        if(N_linesFile >= N_LINE_TOTAL*2) cout<< " DANGER "  <<endl;
        GraphRowData[N_linesFile] = make_pair(a, b);
        GraphRowData[N_linesFile + 1] = make_pair(b, a);    
        N_linesFile += 2;
        
        if(t_max1 < a) t_max1 = a;
        if(t_max1 < b) t_max1 = b;

    }
    // Close the file pointer
    file.close();

    N_Vertices = t_max1 + 1;

    double t2 = omp_get_wtime();

    std::cout << "time Read: " <<  t2 - t1 << std::endl;


    t1 = omp_get_wtime();
    // Sort the Array
    sort(GraphRowData.begin(),GraphRowData.end());

    t2 = omp_get_wtime();

    std::cout << "time Sort: " <<  t2 - t1 << std::endl;


    cout << " N Lines in file  : " << N_linesFile <<endl;
    cout << " N Vertex in file  : " << GraphRowData.size() <<endl;
    cout << " m_Nrows : " << N_Vertices <<endl;


    values.resize(N_linesFile);
    ColPtr.resize(N_linesFile);
    m_Nrow = N_Vertices;
    RowPtr.resize(m_Nrow + 1);
    m_NNZ = N_linesFile;


    std::fill(values.begin(),values.end(),1);
    std::fill(RowPtr.begin(),RowPtr.end(),0);

    RowPtr[0] = 0;

    int* rowPtr = RowPtr.data();
    int* colPtr  = ColPtr.data();

    t1 =  omp_get_wtime();
    
    for (int i = 0; i < N_linesFile; i++)
    {
        colPtr[i] = GraphRowData[i].second;
        rowPtr[GraphRowData[i].first + 1]++;
    }
    t2 =  omp_get_wtime(); 

    std::cout << "time : " << t2 - t1 << std::endl;
 
    for (int i = 0; i < m_Nrow; i++)
        rowPtr[i+1] += rowPtr[i];

    GraphRowData.clear();
}


// ------ Parallel File Read using OpenMP ------- //

void SparseMatrix::SparseMatrixRead1(const char* filename)
{

    int N_Vertices = 0;

    double t1_start = omp_get_wtime();

    std::string command;
    std::string t_file(filename);
    command = "cat ";
    command += t_file + " |  wc -l >> a";
    system(command.c_str());
    ifstream file2;
    file2.open("a");
    int N_LINE_TOTAL = 0;
    file2 >> N_LINE_TOTAL;
    file2.close();
    N_LINE_TOTAL++;
    double t2 = omp_get_wtime();
    system("rm a");
    cout << " LINES IN IFLE : " << N_LINE_TOTAL  << " TIme "  <<  (t2 - t1_start)<<endl;

    // Declare a GLobal Vector
    std::vector<std::pair<int,int>>  GlobalgraphRowData;
    GlobalgraphRowData.resize(2*N_LINE_TOTAL);

    // Declare Variables for each Thread
    std::vector<int> maxRowValue(N_THREADS);

    // Public Variables 
    ifstream file;
    std::ifstream file1;
    std::string temp;
    std::string temp1;
    std::vector<std::pair<int,int>>  GraphRowData;

    omp_set_num_threads(N_THREADS);
    double t1 = omp_get_wtime();
    #pragma omp parallel shared(GlobalgraphRowData) private(file,file1,temp, temp1, GraphRowData)
    { 
        double t01 = omp_get_wtime();
        int threadId = omp_get_thread_num();

        file.open(filename);

        int o_linesPerThread =  N_LINE_TOTAL/N_THREADS;
        if(threadId == N_THREADS-1) o_linesPerThread += (N_LINE_TOTAL % N_THREADS );

        int startPosition;
        startPosition =  (N_LINE_TOTAL/N_THREADS) * threadId ; 

        // int position = fileShifter(startPosition,filename,threadId);

        int position;

        // ---------  Calculate file shift position for each threads  ------- // 
        std::string fileN(filename);
        std::string command ;
        std::string fileout;
        std::string fileoutrm;
        fileout = "b_" + std::to_string(threadId);
        fileoutrm = "rm " + fileout;
        command  = "head -n " + std::to_string(startPosition) + " " + fileN + " | wc -c >> " + fileout;
        
        system(command.c_str());
        file1.open(fileout.c_str());
        file1 >> position;
        file1.close();
        system(fileoutrm.c_str());
        // ----- end of file shift

        file.seekg(position);
        double t1111 = omp_get_wtime();

        //Temp Variables
        int t_val[2];
        int t_max1 = -1000000000;
          // ROw number and the prefix sum of all the rows
        GraphRowData.resize(o_linesPerThread*2);
        int N_linesFile = 0;
        int N_pair = 0;

        // cout << std::setprecision(15) << " Thread : " <<threadId << "  Start time : " << omp_get_wtime() <<endl;
        while (getline(file,temp) && (N_linesFile < o_linesPerThread))
        {
            int t_counter = 0;
            int len = temp.length();
            size_t found = temp.find(" "); 
            int a = stoi(temp.substr(0,found));
            int b = stoi(temp.substr(found+1,len-found));

            GraphRowData[N_pair] = make_pair(a, b);
            GraphRowData[N_pair+1] = make_pair(b, a);    
            N_linesFile++;
            N_pair += 2;

            if(t_max1 < a) t_max1 = a;
            if(t_max1 < b) t_max1 = b;
        }
        file.close();
        double t02 = omp_get_wtime();

        // -- Local sorting
        maxRowValue[threadId] = t_max1;
        std::sort(GraphRowData.begin(),GraphRowData.end());
        

        // Copy the Individual Arrays to the GLobal 
        std::copy(GraphRowData.begin(),GraphRowData.end(),GlobalgraphRowData.begin() + ((N_LINE_TOTAL/N_THREADS) * threadId )*2 );

        GraphRowData.clear();
    }  
    t2 = omp_get_wtime();
    std::cout << " FILE IN : " << t2 - t1 <<endl;

    t1 = omp_get_wtime();

    //  -- Hierarchial Merge Sort of local Arrays ----- // 

    #pragma omp parallel shared(GlobalgraphRowData)
    {
        int threadId = omp_get_thread_num();
        int N_Stages = log2(N_THREADS);

        for (int stage = 0; stage < N_Stages ; stage++)
        {
            double tstart = omp_get_wtime();
            if(threadId % (1 << (stage + 1)) == 0 )
            {   
                int jump    = 1 << stage ; // pow( 2 ^ stage )
                // std::cout << "Stage : "<< stage  << " Thred : " << threadId << " jump : " << jump<< std::endl;

                int start   = ((N_LINE_TOTAL/N_THREADS) * threadId )*2;
                int middle  = ((N_LINE_TOTAL/N_THREADS) * (threadId+jump) )*2;
                int end;
                if(threadId + (jump * 2)== N_THREADS) end = N_LINE_TOTAL*2;
                else                               end = ((N_LINE_TOTAL/N_THREADS) * (threadId+(jump*2)) )*2;

                // std::cout << "Tid : "<< threadId << "Start : " << start << " midd : " << middle << " end : "<<end << std::endl;
                std::inplace_merge(GlobalgraphRowData.begin() + start ,
                                    GlobalgraphRowData.begin() + middle, 
                                    GlobalgraphRowData.begin() + end);
                
            } 
        }
    
    }

    t2 = omp_get_wtime();
    std::cout << " FILE SORT : " << t2 - t1 <<endl;

    // std::for_each(maxRowValue.begin(),maxRowValue.end(),
    //                 [](auto e){std::cout << e<< std::endl;});

    int* n = std::max_element(maxRowValue.data(),maxRowValue.data() + N_THREADS);
    N_Vertices = *n;
    N_Vertices++;
      

    std::cout << " RowPtr : " << N_Vertices << std::endl;
    std::cout << " N Lines in file  : " << N_LINE_TOTAL <<endl;
    std::cout << " N_Vertices : " << N_Vertices <<endl;


    // --- Create Adjacency Matrix -------- //
    values.resize(N_LINE_TOTAL*2);
    ColPtr.resize(N_LINE_TOTAL*2);
    m_Nrow = N_Vertices;
    RowPtr.resize(m_Nrow + 1);
    m_NNZ = N_LINE_TOTAL*2;


    std::fill(values.begin(),values.end(),1);
    std::fill(RowPtr.begin(),RowPtr.end(),0);

    RowPtr[0] = 0;

    int* rowPtr = RowPtr.data();
    int* colPtr  = ColPtr.data();

    #pragma omp parallel for 
    for (int i = 0; i < m_NNZ; i++)
    {
        colPtr[i] = GlobalgraphRowData[i].second;
        rowPtr[GlobalgraphRowData[i].first + 1]++;
    }

 
    for (int i = 0; i < m_Nrow; i++)
        rowPtr[i+1] += rowPtr[i];

    t2 = omp_get_wtime(); 
    std::cout << "----------------- FILE LOAD TIME : " <<  t2 - t1_start << std::endl; 

    GlobalgraphRowData.clear();

}



