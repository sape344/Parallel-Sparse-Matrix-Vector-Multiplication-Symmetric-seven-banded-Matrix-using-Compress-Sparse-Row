#include <iostream>
#include <random>
#include <chrono>

#define M_SIZE 3'000'000
#define NNZ_SIZE 7*(M_SIZE-6)+30

template<typename A,int size>
void CreateSymmetricSevenBandedMatrix(A ** &Matrix);

template<typename A,int size>
void CreateVector(A * &Vector);

template<typename A,int size>
void PrintMatrix(A **Matrix);

template<typename A,int size>
void PrintVector(const A *Vector);

template<typename A, int size>
void ConvertToCSRFormat(A **Matrix,A * &ValVector, int * &RowVector, int* &ColVector);

template<typename A, int size>
void ConvertToCSRFormat(A * &ValVector, int * &RowVector, int* &ColVector);
template<typename A,int size>
A* SpMV(A *Vector, A* ValVector, int * RowVector, int* ColVector);

template<typename A,int size>
A* SpMV(A *Vector, A** Matrix);

template<typename A,int size>
void DeleteMatrix(A **Matrix);

template<typename A,int size>
void DeleteVector(A *Vector);

int main(){
   // double  **Matrix;
    double *ValVector;
    double *Vector;

    int *RowVector;
    int *ColVector;
   // CreateSymmetricSevenBandedMatrix<double,M_SIZE>(Matrix);
    CreateVector<double,M_SIZE>(Vector);

    //ConvertToCSRFormat<double,M_SIZE>(Matrix,ValVector, RowVector, ColVector);
    ConvertToCSRFormat<double,M_SIZE>(ValVector, RowVector, ColVector);

    if (ValVector == nullptr ||RowVector == nullptr ||ColVector == nullptr   )
    {
        std::cout<<"Malloc error mk!!";
    }
    
    
   /* PrintMatrix<double,M_SIZE>(Matrix);
    PrintVector<double,M_SIZE>(Vector);
    PrintVector<double,NNZ_SIZE>(ValVector);
    PrintVector<int,NNZ_SIZE>(ColVector);
    PrintVector<int,M_SIZE>(RowVector);*/
    auto start=std::chrono::high_resolution_clock::now();

    auto Result=SpMV<double,M_SIZE>(Vector,ValVector,RowVector,ColVector);

    auto stop=std::chrono::high_resolution_clock::now();

    
    std::cout << "elapsed time for SpMV : " << std::chrono::duration<double> (stop-start).count() << "\n";

    start=std::chrono::high_resolution_clock::now();
    auto Result2=SpMV<double,M_SIZE>(Vector,ValVector,RowVector,ColVector);
    stop=std::chrono::high_resolution_clock::now();

    std::cout << "elapsed time for Clasic : " << std::chrono::duration<double> (stop-start).count() << "\n";



    int count =0;
    for (int i = 0; i < M_SIZE; i++)
    {
        if(Result[i]==Result2[i]){
            count+=1;
        }
    }
    std::cout<<"Correct rate: "<<count/M_SIZE<<"\n";
    


    //PrintVector<double,M_SIZE>(Result);






   // DeleteMatrix<double,M_SIZE>(Matrix);
    DeleteVector<double,M_SIZE>(ValVector);
    DeleteVector<double,M_SIZE>(Result);
   // DeleteVector<double,M_SIZE>(Result2);


    DeleteVector<double,M_SIZE>(Vector);
    DeleteVector<int,M_SIZE>(RowVector);
    DeleteVector<int,M_SIZE>(ColVector);



    return 0;
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

template<typename A,int size>
void PrintVector( const A *Vector){
     std::cout<< "\n------------------------------------\n";
           for (int j = 0; j < size; j++)
       {

           std::cout<<Vector[j]<<"\t";
           
       }
        std::cout<< "\n------------------------------------\n";
}

template<typename A, int size>
void ConvertToCSRFormat(A **Matrix,A * &ValVector, int * &RowVector, int* &ColVector){
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

template<typename A,int size>
void DeleteVector(A *Vector){
    delete[] Vector;

}
template<typename A,int size>
A* SpMV(A *Vector, A* ValVector, int * RowVector, int* ColVector){
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

template<typename A,int size>
void CreateVector(A *&Vector){
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

template<typename A, int size>
void ConvertToCSRFormat(A * &ValVector, int * &RowVector, int* &ColVector){
    int nonzero_size=6*5+(size-6)*7;
    std::cout<<"nonzero_size= "<<nonzero_size<<"\n";
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

                /*int row=RowVector[i+abs];
                int simetrik=RowVector[i+abs]+indJ-ColVector[row];
                std::cout<<"simetrik = "<<simetrik<<"\n";
                std::cout<<"row = "<<row<<"\n";
                std::cout<<"indJ = "<<indJ<<"\n";
                std::cout<<"First Row = "<<i<<"\tColumn=  "<<ColVector[j]<<"\n";
                std::cout<<"Second Row = "<<i+abs<<"\tColumn=  "<<indJ<<"\n\n\n";*/

                ValVector[RowVector[i+abs]+indJ-ColVector[RowVector[i+abs]]]=ValVector[j];

                }


            }
            
        }
        
    }



}
/*
template<typename A, int size>
void ConvertToCSRFormat(A * &ValVector, int * &RowVector, int* &ColVector){
    int nonzero_size=6*5+(size-6)*7;
    ValVector =new A[nonzero_size];
    ColVector =new int[nonzero_size];
    RowVector =new int[size+1];
    RowVector[0]=0;
    std::uniform_real_distribution<double> unif(0,100);
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
            // ValVector[counter]=static_cast<A>(unif(re));
             ValVector[counter]=counter;

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



}*/