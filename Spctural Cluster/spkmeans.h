

#ifndef FINALPROJECT1_SPKMEANS_H
#define FINALPROJECT1_SPKMEANS_H
int K;
double** centroids;
/**
 * @param rows_num
 * @param cul_num
 * @return initializing matrix
 */
double** Matrix_init(int rows_num , int cul_num);
/**
 * @param mat
 * @param len
 */
void freeMatrix(double** mat, int len);
/**
 * @param x1
 * @param x2
 * @param n
 * @return normalized matrix
 */
double norm(double x1[],double x2[],int vDim);
/**
 *
 * @param x
 * @param n
 * @return adjancy matrix
 */
double** adjancyMatrix(double **x,int n,int vDim);
/**
 *
 * @param W
 * @param n
 * @return diagonal matrix
 */
double** diagonalMatrix(double** W,int n);
/**
 *
 * @param D
 * @param W
 * @param n
 * @return multiply 3 matrix(D^-0.5WD^-0.5)
 */
double ** multiplyMatrix(double** D,double** W,int n);
/**
 *
 * @param mat1
 * @param mat2
 * @param n
 * @return multiply matrix
 */
double** multiplyMatrix2(double** mat1,double** mat2,int n);
/**
 *
 * @param D
 * @param W
 * @param n
 * @return laplacian matrix
 */
double** LaplacianMatrix(double** D,double** W,int n);
/**
 *
 * @param A
 * @param P
 * @param Pt
 * @param n
 *calculates P,Pt matricies 
 */
void Rotation_Matrix(double **A, double **P,double **Pt,int n);
/**
 *
 * @param A
 * @param P
 * @param Pt
 * @param A2
 * @param n
 *Rotating the matricies as described in the algorithm
 */
void Rotate(double **A,double** P,double**Pt,double **A2,int n);
/**
 *
 * @param A
 * @param A2
 * @param n
 * @return checking the convergence for jacobi  algorithms
 */
int convergence(double **A,double **A2,int n);
/**
 *
 * @param eignVectors
 * @param eignValues
 * @param n
 * @return eighvectors matrix with eingvalues in first row
 */
double** merge(double** eignVectors,double* eignValues,int n);
/**
 *
 * @param lnorm
 * @param n
 * @param eignvalues
 * @return jacobi matrix after running the jacobi algorithms
 */
double **Jacobi(double **lnorm,int n,double *eignvalues);
/**
 *
 * @param vectors
 * @param goal
 * @param n
 *@param vDim
 * @return the matrix depend om the goal
 */
double** resultMatrix(double **vectors,char* goal,int n,int vDim);
/**
 *
 * @param i
 * @param j
 * @param tmp
 * @param ind
 * swap the elements in tmp and ind according to i,j
 */
void swap(int i,int j,double* tmp,int * ind);
/**
 *
 * @param eigneValues
 * @param first
 * @param last
 * @param ind
 * @return  sort array of eignValues(sorted by quicksort)
 */
void sorteigen(double*eigneValues,int n, int *ind);
/**
 *
 * @param eigenvalues
 * @param n
 * @param ind
 * @return the max eigen to determine the k
 */
int eigen(double* eigenvalues,int n,int*ind);
/**
 *
 * @param eignvalues
 * @param eignvectors
 * @param ind
 * @param n
 * @return the matrix that contain the vectors of the biggest k eigenvalues
 */
double** UMatrix(int k,double** sliced, int* ind,int n);
/**
 *
 * @param UMatrix
 * @param i
 * @param k
 * @return the norm of the ith vector of UMatrix
 */
double getNorm(double** UMatrix,int i,int k);
/*
  *@param jacobi
  *@param n
  *@return the eignvectors from the jacobi matrix
  */
double** slice(double** jacobi,int n);
/**
 *
 * @param n
 * @param vectors
 * @param k
 * @return normalized the UMatrix
 */
 
double** TMatrix(int k,double** vectors,int n,int vDim);
/**
 *
 * @param centroids
 * @param vector1
 * @param size1
 * @param size2
 * @return the distance of closest point to the centroids
 */
int closest(double **centroids, double *vector1,int size1,int size2);
/**
 *
 * @param arr
 * @param size
 * @return divide the values of the array by the size
 */
double* divide(double *arr,int size);
/**
 *
 * @param centroids1
 * @param centroids2
 * @param size1
 * @param size2
 * @return either to stop the kmeans algorithm or not
 */
int check(double **centroids1,double **centroids2,int size1,int size2);
/**
 *
 * @param k
 * @param TMatrix
 * @param n
 * @param centroids
 *runs the kmean algorithm
 */
void KMEANS(int k,double** TMatrix,int n, double** centroids);
/**
 *
 * @param M
 * @param row_number
 * @param col_number
 * printing the matrix as wanted 
 */
void print_matrix(double **M, int row_number, int col_number);


#endif 
