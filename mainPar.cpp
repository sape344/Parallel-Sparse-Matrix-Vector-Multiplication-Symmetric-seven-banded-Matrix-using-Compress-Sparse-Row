#include <iostream>
#include <random>
#include <chrono>
#include <sys/time.h>
#include "mpi.h"


#define NNZ_SIZE 7*(M_SIZE-6)+30

typedef double deftype;
template<typename A>
void SpMV(A *Vector, A* ValVector, int * RowVector, int* ColVector,int size, A*& Result);

template<typename A,int size>
void CreateSymmetricSevenBandedMatrix(A ** &Matrix);

template<typename A>
void CreateVector(A * &Vector,int size);

template<typename A,int size>
void PrintMatrix(A **Matrix);

template<typename A>
void PrintVector( const A *Vector,int size);

template<typename A,int size>
void PrintVector(const A *Vector);

template<typename A, int size>
void ConvertToCSRFormat(A **Matrix,A * &ValVector, int * &RowVector, int* &ColVector);

template<typename A>
void ConvertToCSRFormat(A * &ValVector, int * &RowVector, int* &ColVector, int size);
template<typename A>
A* SpMV(A *Vector, A* ValVector, int * RowVector, int* ColVector,int size);

template<typename A,int size>
A* SpMV(A *Vector, A** Matrix);

template<typename A>
A* SpMVParellel(A *ValAVector, A* BVector,int rank,int localsize,int nproc);

template<typename A,int size>
void DeleteMatrix(A **Matrix);

template<typename A>
void DeleteVector(A *Vector);
int** ArraysForScatterv(int N,int nproc);


int nproc,rank;

int main(int argc, char **argv){
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    unsigned long long SIZE=std::atoi(argv[1])*1.25e6;
    int** Matrix=ArraysForScatterv(SIZE,nproc); //Needed Matrix to use Scatterv and gatterv  
    int local_size=Matrix[5][rank];

    //std::chrono::_V2::system_clock::time_point start,stop;
    double startTime,stopTime;
    double scatterDiffTime;
    double gatherDiffTime;
    double compDiffTime;

    //Pointers for produce Random matrix and vector
    deftype *AValVector;
    deftype *BVector;
    deftype *CResult;
    int *RowVector;
    int *ColVector;

    //Pointers for Local values 
    deftype *LocAValVector;
    deftype *LocBVector;
    deftype *LocCResult;
    LocAValVector = new deftype[Matrix[1][rank]];
    LocBVector= new deftype[Matrix[3][rank]];

    MPI_Request requestA,requestB;
    MPI_Status status;
    if(rank ==0 ){

        deftype *Result= new deftype[SIZE];
        CreateVector<deftype>(CResult,SIZE);
        CreateVector<deftype>(BVector,SIZE);
        ConvertToCSRFormat<deftype>(AValVector, RowVector, ColVector,SIZE);
       
    if(nproc==1){

        startTime=MPI_Wtime();
        auto CResult2=SpMV<deftype>(BVector,AValVector,RowVector,ColVector,SIZE);
        stopTime=MPI_Wtime();
        auto diffSerial=stopTime-startTime;

        std::cout<<"Simple Matrix multiplication WallTime(s): " <<stopTime-startTime<<" Size: "<<SIZE<<" nprocs: "<<nproc << std::endl;
        DeleteMatrix<int,6>(Matrix);
        DeleteVector(AValVector);
        DeleteVector(CResult);
        DeleteVector(CResult2);
        DeleteVector(BVector);
        DeleteVector(ColVector);
        DeleteVector(RowVector); 
        MPI_Finalize();
        return 0;
    }

    if (rank==0)
    {
        startTime=MPI_Wtime();

    }
    
    MPI_Scatterv(AValVector,Matrix[1],Matrix[0],MPI_DOUBLE,LocAValVector,Matrix[1][rank],MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Scatterv(BVector,Matrix[3],Matrix[2],MPI_DOUBLE,LocBVector,Matrix[3][rank],MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0)
    {
        stopTime=MPI_Wtime();
        scatterDiffTime=stopTime-startTime;
        startTime=MPI_Wtime();
    }
    


    
    LocCResult=SpMVParellel(LocAValVector, LocBVector, rank, local_size, nproc);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0)
    {
        
        stopTime=MPI_Wtime();
        compDiffTime=stopTime-startTime;
        startTime=MPI_Wtime();
    }
    
    MPI_Gatherv(LocCResult,local_size,MPI_DOUBLE,CResult,Matrix[5],Matrix[4],MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);


    if (rank==0)
    {
       // MPI_Wait(&requestB,&status);

        stopTime=MPI_Wtime();
        gatherDiffTime= stopTime-startTime;


    }


    if (rank ==0 )
    {   

      //  stop=std::chrono::high_resolution_clock::now();
//        std::cout << "elapsed time for parelel : " << std::chrono::duration<double> (stop-start).count() << "\n";
        std::cout<<"Parelel Matrix multiplication WallTime(s): " <<gatherDiffTime+compDiffTime+scatterDiffTime<<
        " GatherTime(s): "<<gatherDiffTime<<" CompTime(s): "<<compDiffTime<<" Scattertime(s): "<<scatterDiffTime<<" Size: "<<SIZE<<" nprocs: "<<nproc << std::endl;
        //start=std::chrono::high_resolution_clock::now();
        /*startTime=MPI_Wtime();

        auto CResult2=SpMV<deftype>(BVector,AValVector,RowVector,ColVector,SIZE);
        stopTime=MPI_Wtime();
        auto diffSerial=stopTime-startTime;

        //stop=std::chrono::high_resolution_clock::now();
        //std::cout << "elapsed time for Clasic : " << std::chrono::duration<double> (stop-start).count() << "\n";
        std::cout<<"Simple Matrix multiplication WallTime(s): " <<stopTime-startTime<<" Size: "<<SIZE << std::endl;
        std::cout<<"Diff rate Serial/Parelel ="<<diffSerial/compDiffTime<<"\n";
        
        //PrintVector(CResult,SIZE);
        //PrintVector(CResult2,SIZE);
        for (int i = 0; i < SIZE; i++)
        {
            if(CResult2[i]!=CResult[i]){
                std::cout<<"Hatalı! Size"<<SIZE<<"\n";
                break;
            }
        }
        DeleteVector(CResult2);*/
        


    }    



    if (rank==0)
    {
    DeleteVector(AValVector);
    DeleteVector(CResult);
    DeleteVector(BVector);
    DeleteVector(ColVector);
    DeleteVector(RowVector);    

    }
    DeleteVector(LocCResult);
    DeleteVector(LocBVector);
    DeleteVector(LocAValVector);
    DeleteMatrix<int,6>(Matrix);
    MPI_Finalize();
    return 0;

    



}
}

template<typename A,int size>
void CreateSymmetricSevenBandedMatrix(A ** &Matrix){
    Matrix= new A*[size];
    for (int i = 0; i < size; i++)
    {
        Matrix[i]=new A[size];
    }
    std::uniform_real_distribution<double> unif(1,1000);
    std::default_random_engine re;
    for (int i = 0,counter=0; i < size; i++)
    {
       for (int j = 0; j < size; j++)
       {
           //if (std::abs(i-j)<4 )
           if ((i-j<4 )&&(-1<i-j))
           
           {
              Matrix[i][j]=static_cast<A>(unif(re));
              // Matrix[i][j]=counter;

               Matrix[j][i]=Matrix[i][j];
               counter+=1;
           }
           else{
               Matrix[i][j]=0;
           }
           
       }
       
    }
}



template<typename A,int size>
void PrintMatrix(A **Matrix){
    for (int i = 0; i < size; i++)
    {
       for (int j = 0; j < size; j++)
       {

           std::cout<<Matrix[i][j]<<"\t";
           
       }
           std::cout<<"\n";
       
    }


}

template<typename A>
void PrintVector( const A *Vector,int size){
     std::cout<< "\n------------------------------------\n";
           for (int j = 0; j < size; j++)
       {

           std::cout<<Vector[j]<<"\t";
           
       }
        std::cout<< "\n------------------------------------\n";
}

template<typename A>
void ConvertToCSRFormat(A **Matrix,A * &ValVector, int * &RowVector, int* &ColVector, int size){
    int nonzero_size=6*5+(size-6)*7;
    ValVector =new A[nonzero_size];
    ColVector =new int[nonzero_size];
    RowVector =new int[size+1];


    int count=0;
    for (int i = 0; i < size; i++)
    {
        RowVector[i]=count;
       for (int j = 0; j < size; j++)
       {

           if(Matrix[i][j]!=0){
                
                ValVector[count] =Matrix[i][j];
                ColVector[count]=j;
                count+=1;
           }
       }
       
    }
    RowVector[size]=count;

}


template<typename A,int size>
void DeleteMatrix(A **Matrix){
for (int i = 0; i < size; i++)
{
    delete[] Matrix[i];
}
delete[] Matrix;

}

template<typename A>
void DeleteVector(A *Vector){
    if (Vector!=nullptr)
    {
           delete[] Vector;

    }
    

}
template<typename A>
A* SpMV(A *Vector, A* ValVector, int * RowVector, int* ColVector,int size){
    A *Result= new A[size];



    int in_row=0;
    for (int i = 0; i < size; i++)
    {        
        Result[i]=0;
     for (int j  = RowVector[i]; j < RowVector[i+1]; j++)
     {

        Result[i]+= Vector[ColVector[j]]*ValVector[j];
     }
     

    }
    
    return Result;
}


template<typename A>
void SpMV(A *Vector, A* ValVector, int * RowVector, int* ColVector,int size, A*& Result){
    //A *Result= new A[size];



    int in_row=0;
    for (int i = 0; i < size; i++)
    {        
        Result[i]=0;
     for (int j  = RowVector[i]; j < RowVector[i+1]; j++)
     {

        Result[i]+= Vector[ColVector[j]]*ValVector[j];
     }
     

    }
    
}

template<typename A>
void CreateVector(A * &Vector,int size){
     Vector= new A[size];
    std::uniform_real_distribution<double>  unif(1,1000);
    std::default_random_engine re;

    for (int i = 0; i < size; i++)
    {
       Vector[i]=static_cast<A>(unif(re));
    }
    
}

template<typename A,int size>
A* SpMV(A *Vector, A** Matrix){
     A* Result = new A[size];
    for (size_t i = 0; i < size; i++)
    {
        Result[i]=0;
        for (int j = 0; j < size; j++)
        {
            Result[i]+=Matrix[i][j]*Vector[j];
        }
        
    }
    
    return Result;


}

template<typename A>
void ConvertToCSRFormat(A * &ValVector, int * &RowVector, int* &ColVector, int size){
    int nonzero_size=6*5+(size-6)*7;
    ValVector =new A[nonzero_size];
    ColVector =new int[nonzero_size];
    RowVector =new int[size+1];
    RowVector[0]=0;
    std::uniform_real_distribution<double> unif(1,1000);
    std::default_random_engine re;

       for (int i = 0,N=0,M=4,counter=0,jUpper; i < size; i++)
    {
        jUpper=M+N;
        if (M+N>size)
        {
            jUpper=size;
        }
        
        for (int j = N; j < jUpper ; j++)
        {
             ColVector[counter]=j;
             counter+=1;
        }
        RowVector[i+1]=counter;

        if (M<7)
        {
            M+=1;
        }
        else{
             N+=1;
        }      
        
    }

    for (int i = 0,counter=0; i < size; i++)
    {
      

        for (int j = RowVector[i],indJ; j < RowVector[i+1]; j++)
        {

            
            if (ColVector[j]==i)
            {

               ValVector[j]=static_cast<double>(unif(re));

                indJ =ColVector[j];
                
            }
            else if(ColVector[j]>i){
               ValVector[j]=static_cast<double>(unif(re));
                
                int abs= ColVector[j]-indJ;


                if(i+abs< size && indJ<size){
                ValVector[RowVector[i+abs]+indJ-ColVector[RowVector[i+abs]]]=ValVector[j];

                }


            }
            
        }
        
    }



}

int** ArraysForScatterv(int N,int nproc){

    //N matrix size
    //nprocs kullanılan cpu

    int localN= N/nproc;  //local size
    int rem=  N-localN*nproc;   //remander



    int **Matrix= new int*[6];
    for (int i = 0; i < 6; i++)
    {
        Matrix[i]=new int[nproc];
    }
    /*
    Matrix[0] idxA
    Matrix[1] countSendA
    Matrix[2] idxB
    Matrix[3] countSendB
    Matrix[4] idxC
    Matrix[5] countSendC- local size
    */

    for (int i = 0; i < nproc; i++)
    {
        if (i==0)
        {
            Matrix[0][i]=0;
            Matrix[2][i]=0;
            Matrix[4][i]=0;
            Matrix[1][i]=7*localN-6;
            Matrix[3][i]=localN+3;
            Matrix[5][i]=localN;


        }
        else if(i==nproc-1){
            Matrix[0][i]=Matrix[0][i-1]+Matrix[1][i-1];
            Matrix[2][i]=Matrix[2][i-1]+Matrix[3][i-1]-6;
            Matrix[4][i]=Matrix[4][i-1]+Matrix[5][i-1];
            Matrix[1][i]=7*localN-6;
            Matrix[3][i]=localN+3;
            Matrix[5][i]=localN;

        }
        else{
            Matrix[0][i]=Matrix[0][i-1]+Matrix[1][i-1];
            Matrix[2][i]=Matrix[2][i-1]+Matrix[3][i-1]-6;
            Matrix[4][i]=Matrix[4][i-1]+Matrix[5][i-1];
            Matrix[1][i]=7*localN;
            Matrix[3][i]=localN+6;
            Matrix[5][i]=localN;  
        }
        if (rem>0)
        {
            rem-=1;
            Matrix[1][i]+=7;
            Matrix[3][i]+=1;
            Matrix[5][i]+=1;  

        }
    }

    
    return Matrix;
}

template<typename A>
A* SpMVParellel(A *ValAVector, A* BVector,int rank,int localsize,int nproc){
    A* C= new deftype[localsize];
    if (rank==0)
    {   int j=0;
        C[0]=0;
        for (int i=0 ; i < 4; i++,j++)
        {
            C[0]+=ValAVector[j]*BVector[i];
        }
        C[1]=0;
        for (int i=0  ; i < 5; i++,j++)
        {
            C[1]+=ValAVector[j]*BVector[i];
        }
        C[2]=0;
        for (int i=0  ; i < 6; i++,j++)
        {
            C[2]+=ValAVector[j]*BVector[i];
        }
        for (int k = 3; k < localsize; k++)
        {
            C[k]=0;
            for (int i = k-3; i < k+4; i++,j++)
            {
                C[k]+=ValAVector[j]*BVector[i];
            }            
        }              
    }
    else if (nproc-1==rank)
    {   int j=0;

        for (int k = 0; k < localsize-3; k++)
        {
            C[k]=0;
            for (int i = k; i < k+7; i++,j++)
            {
                C[k]+=ValAVector[j]*BVector[i];
            }            
        }
        C[localsize-3]=0;
        for (int i=localsize-3  ; i < localsize+3; i++,j++)
        {
            C[localsize-3]+=ValAVector[j]*BVector[i];
        }
        C[localsize-2]=0;
        for (int i=localsize-2  ; i < localsize+3; i++,j++)
        {
            C[localsize-2]+=ValAVector[j]*BVector[i];
        }
        C[localsize-1]=0;
        for (int i=localsize-1  ; i < localsize+3; i++,j++)
        {
            C[localsize-1]+=ValAVector[j]*BVector[i];
        }

              
    }

    else{
        int j=0;
        for (int k = 0; k < localsize; k++)
        {
            C[k]=0;
            for (int i = k; i < k+7; i++,j++)
            {
                C[k]+=ValAVector[j]*BVector[i];
            }            
        }
    }
    return C;

}