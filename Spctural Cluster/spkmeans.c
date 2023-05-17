#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "spkmeans.h"

int main(int argc,char** argv) {
   
    int i;
    int j;
    char c;
    double **vectors;
    int vDim;
    double charc;
    FILE *file;
    int n;
    double** result;
     if(argc <0 || argc >3){
    exit(1);
    }
    file = fopen(argv[2],"r");
    if(file==NULL)
    {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    n=0;
    vDim=0;
    while((c=fgetc(file))!=EOF)
    {
        if(c=='\n')
        {
            n+=1;
            vDim+=1;
        }
        if(c==',')
        {
            vDim+=1;
        }
    }
    rewind(file);
    vectors=(double**)malloc(n*sizeof(double*));
    if(vectors==NULL)
      { printf("An Error Has Occurred\n");
        exit(1);}
    vDim=vDim/n;
    for(i=0;i<n;i++)
    {
        vectors[i]=(double*)malloc(vDim*sizeof(double));
        if(vectors[i]==NULL){
              printf("An Error Has Occurred\n");
              exit(1);}
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<vDim;j++)
        {
            fscanf(file,"%lf",&charc);
            vectors[i][j]=(double) charc;
            fgetc(file);
        }
    }
    fclose(file);
    result= resultMatrix(vectors,argv[1],n,vDim);
    freeMatrix(vectors,n);
    if(strcmp(argv[1],"jacobi")==0){
        print_matrix(result,n+1,n);
        freeMatrix(result,n+1);
        
    }
    else
    {
        print_matrix(result,n,n);
        freeMatrix(result,n);
        
    }
  return 1;

}
double** Matrix_init(int rows_num , int cul_num){
    double** mat;
    int i,j;
    mat= (double**) malloc(sizeof(double*)*rows_num);
    if(mat == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i=0;i<rows_num;i++){
        mat[i]= (double*) malloc(sizeof(double)*cul_num);
        if(mat[i] == NULL) {
            printf("An Error Has Occurred\n");
            exit(1);
        }
        for (j=0;j<cul_num;j++){
            mat[i][j]=0;
        }
    }
    return mat;
}

void freeMatrix(double** mat, int len){
    int i;
    for (i = 0; i < len; i++) {
        free(mat[i]);
    }
    free(mat);
}
double norm(double *x1,double *x2,int vDim){
    double result;
    int i;
    result=0;
    for(i=0;i<vDim;i++)
    {
        result+=pow(x1[i]-x2[i],2);
    }
    result=sqrt(result);
    return result;
}

double** adjancyMatrix(double **x,int n,int vDim){
    double**W;
    int i,j;
    W=Matrix_init(n,n);
    for( i=0;i<n;i++)
    {
        for( j=0;j<n;j++)
        {
            if(i!=j){
              W[i][j]=exp(-((norm(x[i],x[j],vDim)/2)));
                
            }
            
        }
    }
    return W;
}


double** diagonalMatrix(double** W,int n) {   
    double**D;
    int i,j;
    D=Matrix_init(n,n);
    for ( i = 0; i < n; i++)
    {
        for ( j = 0; j < n; j++)
        {
            D[i][i] += W[i][j] ;
        }
    }
    return D;
}

double ** multiplyMatrix(double** D,double** W,int n){
    double**DW,**DWD;
    int i,j,k;
    DW = Matrix_init(n,n);
    DWD = Matrix_init(n,n);
    for ( i = 0; i < n; i++)
    {
        D[i][i]=pow(D[i][i],-0.5);
    }
    for (i = 0; i < n; i++)
    {
        for ( j = 0; j < n; j++)
        {
            for( k=0;k<n;k++)
            {
                DW[i][j]+=D[i][k]*W[k][j];
            }
        }
    }
    for ( i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for( k=0;k<n;k++)
            {
                DWD[i][j]+=DW[i][k]*D[k][j];
            }
        }
    }
    freeMatrix(DW,n);
    return DWD;
}
double** multiplyMatrix2(double** mat1,double** mat2,int n){
    double **result ;
    int i,j,k;
    result=Matrix_init(n,n);
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            for(k=0;k<n;k++){
                result[i][j]+=mat1[i][k]*mat2[k][j];
            }
        }
    }
    return result;
}
double** LaplacianMatrix(double** D,double** W,int n){
    double**L;
    int i;
    int j;
    L=multiplyMatrix(D,W,n);
    for(i=0;i<n;i++)
    {
        L[i][i]=1-L[i][i];
    }
    for(i=0;i<n;i++){
      for(j=0;j<n;j++){
          if(i!=j)
          {L[i][j]=-L[i][j];}
        }
    }
    return L;
}

void Rotation_Matrix(double **A, double **P,double **Pt,int n)
{
    int i,j,sign;
    int row,col;
    double max,c,theta,t,s;
    max=A[0][1];
    row=0;
    col=1;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if(fabs(A[i][j])>fabs(max) && i!=j)
            {
                max=A[i][j];
                row=i;
                col=j;
            }
            if(i==j)
            {
                P[i][j]=1;
                Pt[i][j]=1;
            }
        }
    }
    theta=(A[col][col]-A[row][row])/(2*(A[row][col]));
    if(theta>=0){sign=1;}
    else{sign=-1;}
    t=sign/(fabs(theta)+sqrt(pow(theta,2)+1));
    c=1/sqrt(pow(t,2)+1);
    s=t*c;
    P[row][row]=c;
    P[row][col]=s;
    P[col][row]=-s;
    P[col][col]=c;
    Pt[row][row]=c;
    Pt[row][col]=-s;
    Pt[col][row]=s;
    Pt[col][col]=c;
    
}
void Rotate(double **A,double** P,double**Pt,double **A2,int n)
{   int i,j;
    double** tmp,**tmp2;
    tmp=multiplyMatrix2(Pt,A,n);   
    tmp2=multiplyMatrix2(tmp,P,n);
    for(i=0;i<n;i++){
      for(j=0;j<n;j++){
        A2[i][j]=tmp2[i][j];
      }
    }
    freeMatrix(tmp2,n);
    freeMatrix(tmp,n);
    
}  
int convergence(double **A,double **A2,int n)
{
    int i,j;
    double sum1,sum2;
    sum1=0;
    sum2=0;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if(i!=j)
            {
                sum1+=pow(A[i][j],2);
                sum2+=pow(A2[i][j],2);
            }
        }
    }
    if(sum1-sum2<=pow(10,-5)){return 1;}
    return 0;

}
double **Jacobi(double **lnorm,int n,double *eignvalues)
{
    double **P,**Pt,**A,**A2,**eignVectors,**tmp,**tmp2;
    int i,j,iter;
    iter=0;
    A=lnorm;
    eignVectors=Matrix_init(n,n);
    
    for(i=0;i<n;i++)
    {
        eignVectors[i][i]=1;
    }
    P=Matrix_init(n,n);
    Pt=Matrix_init(n,n);
        while(iter<100)
    {
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++){
                P[i][j]=0;
                Pt[i][j]=0;
            }
        }
        Rotation_Matrix(A,P,Pt,n);
        tmp=multiplyMatrix2(eignVectors,P,n);
        freeMatrix(eignVectors,n);
        eignVectors=tmp;   
        A2=Matrix_init(n,n);
        Rotate(A,P,Pt,A2,n);    
        if(convergence(A,A2,n)==1)
        {
            freeMatrix(A,n);
            A=A2;
            break;
        }
        freeMatrix(A,n);
        A=A2;
        iter++;
    }    
    
    for(i=0;i<n;i++)
    {
        eignvalues[i]=A[i][i];
    }
    freeMatrix(A,n);
    freeMatrix(P,n);
    freeMatrix(Pt,n);
    tmp2=merge(eignVectors,eignvalues,n);
    freeMatrix(eignVectors,n);
    return tmp2;
}
double** merge(double** eignVectors,double* eignValues,int n)
{
    int i,j;
    double** tmp;
    tmp=Matrix_init(n+1,n);
    for(i=0;i<=n;i++){
        for(j=0;j<n;j++) {
            if (i == 0) {
                tmp[i][j] = eignValues[j];
            } else {
                tmp[i][j] = eignVectors[i - 1][j];
            }
        }
    }
    return tmp;
}
double** resultMatrix(double **vectors,char* goal,int n,int vDim)
{
    double **W,**D,**L,**eignVectors;
    double* eignvalues;
    W=adjancyMatrix(vectors,n,vDim);
    if(strcmp(goal,"wam")==0){
        return W;
    }
    D=diagonalMatrix(W,n);
    if(strcmp(goal,"ddg")==0){
        freeMatrix(W,n);
        return D;

    }
    L=LaplacianMatrix(D,W,n);
    freeMatrix(W,n);
    freeMatrix(D,n);
    
    if(strcmp(goal,"lnorm")==0){
        return L;
    }
    else{
        eignvalues=(double*)malloc(n*sizeof(double));
        if(eignvalues==NULL){
              printf("An Error Has Occurred\n");
              exit(1);
        }
        eignVectors=Jacobi(vectors,n,eignvalues);
        freeMatrix(L,n);
        free(eignvalues);
        return eignVectors;
    }
}
void swap(int i,int j,double* tmp,int * ind){
    double tmp1;
    int tmp2;
    tmp1=tmp[i];
    tmp[i]=tmp[j];
    tmp[j]=tmp1;
    tmp2=ind[i];
    ind[i]=ind[j];
    ind[j]=tmp2;
}
void sorteigen(double*tmp,int n, int *ind){
    int i,j;
    
    for(i=0;i<n;i++)
    {
      for (j=i+1;j<n;j++)
      {
        if(tmp[i]<tmp[j])
        {
          swap(i,j,tmp,ind);
        }
      }
    }
     
}
int eigen(double* eigenvalues,int n,int*ind){
    double* tmp; 
    int i;
    int index;
    double max;
    tmp=(double*)malloc(n*sizeof(double));
    for( i=0;i<n;i++)
    {
        tmp[i]=eigenvalues[i];
    }
    sorteigen(tmp,n,ind);
    index=0;
    max=fabs(tmp[0]-tmp[1]);
    for(i=0;i<floor(n/2);i++)
    {
      if (fabs(tmp[i]-tmp[i+1])>max)
      {
          max=fabs(tmp[i]-tmp[i+1]);
          index=i;
      }
    }

    return index+1;
}
double** UMatrix(int k,double** sliced, int* ind,int n){
    int i,j;
    double** U2;
    U2=Matrix_init(n,k);
    for(i=0;i<n;i++)
    {
        for(j=0;j<k;j++)
        {
            U2[i][j]=sliced[i][ind[j]];
        }
    } 
    free(ind);   
    return U2;
}
double getNorm(double** U,int i,int k)
{
    double sum;
    int j;
    sum=0;
    for(j=0;j<k;j++){
        sum+=pow(U[i][j],2);
    }
    return sum;
}
double** slice(double** jacobi,int n){
    int i,j;
    double** result;  
    result= Matrix_init(n,n);
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++){
            result[i][j]=jacobi[i+1][j];
       }
    }
    return result;
}

double** TMatrix(int k,double** vectors,int n,int vDim){
    double** D,**L,**W,**J;
    double** T,**U;
    double* eignvalues2;
    int *ind;
    int i,j;
    ind=(int*) malloc(n*sizeof (int));
    if(ind==NULL){  
        printf("An Error Has Occurred\n");
        exit(1);}
     for(i=0;i<n;i++)
    {
        ind[i]=i;
    }    
    eignvalues2=(double *) malloc(n*sizeof (double));
    if(eignvalues2==NULL){  
        printf("An Error Has Occurred\n");
        exit(1);}
    W= adjancyMatrix(vectors,n,vDim);
    D= diagonalMatrix(W,n);
    L= LaplacianMatrix(D,W,n);
    J= Jacobi(L,n,eignvalues2);
    if(k==0){
        k= eigen(J[0],n,ind);
        K=k;
    }
    else{
    sorteigen(J[0],n,ind);
    }
    U=UMatrix(k,slice(J,n),ind,n);
    free(eignvalues2);
    T= Matrix_init(n,k);
  for(i=0;i<n;i++){
        for(j=0;j<k;j++)
        {
          if(getNorm(U,i,k)!=0)
            { T[i][j]=((U[i][j])/(sqrt(getNorm(U,i,k))));}
          else{T[i][j]=U[i][j];}
        }
    }
    
    
    return T;

}

int closest(double **centroids2, double *vector1,int size1,int size2){
    double close;
    int result;
    int i,j;
    close=0;
    result=0;
    for( i=0;i<size1;i++ ){
        close+=(centroids2[0][i]-vector1[i])*(centroids2[0][i]-vector1[i]);
    }
    for( i=0;i<size2;i++){
        double tmp=0;
        for( j=0;j<size1;j++){
            tmp+=(centroids2[i][j]-vector1[j])*(centroids2[i][j]-vector1[j]);
        }
        if(tmp<close){
            close=tmp;
            result=i;
        }
    }
    return result;
}
double* divide(double *arr,int size){
    double *result;
    int i;
    result=(double*)malloc((size)*sizeof(double));
    if(result == NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(i=0 ;i<size;i++){
        result[i]=arr[i]/arr[size];
    }
    return result;
}
int check(double **centroids1,double **centroids2,int size1,int size2){
    double sum;
    int i;
    int j;
    for ( i = 0; i < size1; i++)
    {
        sum=0.0;

        for (j = 0; j < size2; j++)
        {
            sum+=pow((centroids1[i][j]-centroids2[i][j]),2);
        }
        if((sum*sum)>0){
            return -1;
        }
    }

    return 1;
}

void KMEANS(int k,double** T,int n,double** centroids3)
{
    int flag,q,close;
    double** matrix;
    double**centroids2;
    int i,j;
    flag=-1;
    q=0;
    centroids= (double**)malloc(k*sizeof (double*));
    for(i=0;i<k;i++)
    {
        centroids[i]=(double*) malloc(k*sizeof (double ));
    }
    centroids2= (double**)malloc(k*sizeof (double*));
    for(i=0;i<k;i++)
    {
        centroids2[i]=(double*) malloc(k*sizeof (double ));
    }
    for(i=0;i<k;i++)
        {
            for(j=0;j<k;j++){
                centroids[i][j]=centroids3[i][j];
            }
        }
    freeMatrix(centroids3,k);
    matrix=Matrix_init(k,k+1);
    while(q<300 && flag==-1)
    {
        for(i=0;i<n;i++)
        {
            close= closest(centroids,T[i],k,k);
            for(j=0;j<k;j++)
            {
                matrix[close][j]+=T[i][j];
                if(j==0){
                    matrix[close][k]+=1;
                }
            }
        }
        for(i=0;i<k;i++)
        {
            for(j=0;j<k;j++){
                centroids2[i][j]=centroids[i][j];
            }
        }
        for(i=0;i<k;i++)
        {
            free(centroids[i]);
            centroids[i]= divide(matrix[i],k);
        }
        flag= check(centroids,centroids2,k,k);
        q+=1;
    }
    freeMatrix(T,n);
    freeMatrix(centroids2,k);
    freeMatrix(matrix,k);
    


}
void print_matrix(double **M, int row_number, int col_number){
    int i,j;
    for(i=0;i<row_number;i++){
        for(j=0;j<col_number;j++){
            printf("%.4f",M[i][j]);
            if(j!=col_number-1) printf(",");
        }
        printf("\n");
    }
}