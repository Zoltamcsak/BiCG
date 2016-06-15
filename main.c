#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>
 
double* matrixVectorMultiplication(double** matrix, double* vector, int n) {
    double* result = (double*)malloc(n*sizeof(double));
    int i, j;
    for (i = 0; i<n; i++) {
        for (j=0; j<n; j++) {
            result[i] += matrix[i][j]*vector[j];
        }
    }
    return result;
}
 
double* matvecprod(double* A, double* v, int N) {
    double alpha = 1.0, beta = 0.0;
    char no = 'N', tr = 'T';
    double* u = (double*)malloc(N*sizeof(double));
    int m = N, n = N, lda = N, incx = 1, incy = 1, i;
    dgemv_(&tr, &m, &n, &alpha, A, &lda, v, &incx, &beta, u, &incy);
    return u;
}
 
double* vectorMatrixMultiplication(double *vector, double **matrix, int n) {
    double *result = (double*) malloc(n*sizeof(double));
    int i,j;
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            result[i] += vector[j]*matrix[j][i];
        }
    }
    return result;
}
 
double* vecmatprod(double* v, double* A, int N) {
    double alpha= 1.0, beta= 0.0;
    char no= 'N', tr= 'T';
    int m= N, n= 1, k= N, lda= N, ldb= N, ldc= N;
    double* u = (double*)malloc(N*sizeof(double));
    dgemm_(&no,&no,&m,&n,&k,&alpha,A,&lda,v,&ldb,&beta,u,&ldc);
    return u;
}
 
double* vectorSubtraction(double* vec1, double* vec2, int n) {
    int i;
    double *result = (double*)malloc(n*sizeof(double));
    for (i=0; i<n; i++) {
        result[i] = vec1[i] - vec2[i];
    }
    return result;
}
 
double* vectorAddition(double* vec1, double* vec2, int n) {
    int i;
    double *result = (double*)malloc(n*sizeof(double));
    for (i=0; i<n; i++) {
        result[i] = vec1[i] + vec2[i];
    }
    return result;
}
 
double vectorVectorMultiplication(double* vec1, double* vec2, int n) {
    int i;
    double result = 0;
    for (i=0; i<n; i++) {
        result += vec1[i] * vec2[i];
    }
    return result;
}
 
double vecvecprod(double* v, double* A, int N) {
    double alpha = 1.0, beta = 0.0;
    char no = 'N', tr = 'T';
    double* u = (double*)malloc(N*sizeof(double));
    int m = N, n = 1, lda = N, incx = 1, incy = 1, i;
    dgemv_(&tr, &m, &n, &alpha, A, &lda, v, &incx, &beta, u, &incy);
    return u[0];
}
 
double* vectorNumMultiplication(double* vec, double num, int n) {
    int i;
    double* resultVector = (double*)malloc(n*sizeof(double));
    for (i=0; i<n; i++) {
 
        resultVector[i] = vec[i] * num;
    }
    return resultVector;
}
 
double* transposeMatrix(double* mat, int n) {
    double* transpose = (double*)malloc(n*n*sizeof(double *));
    int i, j;
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            transpose[i*n+j] = mat[j*n+i];
        }
    }
    return transpose;
}
 
double* inverseMatrix(double* mat, int n) {
    int *PIV = (int*) malloc(n*sizeof(int));
    int LWORK = n*n;
    double *work = (double*) malloc(n*sizeof(double));;
    int info;
    dgetrf_(&n, &n, mat, &n, PIV, &info);
    dgetri_(&n, mat, &n, PIV, work, &LWORK, &info);
    return mat;
}
 
void luFactorization(int n, double* a, double* b) {
    int k,i,j;
    double* c = (double*)malloc(n*sizeof(double));
    double* m = (double*)malloc(n*n*sizeof(double *));
    for (i=0; i<n; i++) {
        m[i*n+i] = 1;
    }
    for (k=0; k<n-1; k++){
        if (a[k*n+k] == 0) {
            break;
        }
 
        for (i=k+1; i<n; i++) {
          m[i*n+k] = (a[i*n+k]/a[k*n+k]);
          b[i] = b[i] - (m[i*n+k]*b[k]);
            for (j=k+1; j<n; j++) {
                a[i*n+j] = a[i*n+j] - (m[i*n+k]*a[k*n+j]);
            }
        }
    }
}
 
double* backwardSubstitution(double *upperMatrix, double *y, int n) {
    int i, j;
    double* resultVector = (double*)malloc(n*sizeof(double));
    for (i=n-1; i>=0; i--) {
        resultVector[i] = y[i];
        for (j=i+1; j<n; j++) {
            resultVector[i] -= upperMatrix[i*n+j] * resultVector[j];
        }
        resultVector[i] /= upperMatrix[i*n+i];
    }
    return resultVector;
}
 
double *copyMatrix(double* A, int n) {
    int i,j;
    double* copyRes = (double*)malloc(n*n*sizeof(double *));
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            copyRes[i*n+j] = A[i*n+j];
        }
    }
    return copyRes;
}
 
double *cop(double *A, int n) {
    int i, j, x = 0;
    double *result = (double*)malloc(n*sizeof(double));
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            result[x] = A[i*n+j];
            x++;
        }
    }
    return result;
}
 
double **backcop(double *A, int n) {
    int i,j,x = 0;
    double** copyRes = (double**)malloc(n*sizeof(double *));
    for (i=0; i<n; i++) {
        copyRes[i] = (double*)malloc(n*sizeof(double));
        for (j=0; j<n; j++) {
            copyRes[i][j] = A[x];
            x++;
        }
    }
    return copyRes;
}
 
void printMat(int n, double* a) {
    int i,j;
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            printf("%f ", a[i*n+j]);
        }
        printf("\n");
    }
}
 
void printVector(double* a, int n) {
    int i;
    for (i=0; i<n; i++) {
        printf("%f \n", a[i]);
    }
}
 
void BiCG(double *A, double *b, int n) {
    double *r0, *x0, *r0Star, *z0, *z0Star, alphaK, *xKplus1, *xKplus1Star, *rKplus1, *rKplus1Star, *pKplus1, *pKplus1Star, *x1;
    double rho1, rho2, beta, *p0, *p0Star, *q, *qStar;
    double *M, *Minverse;
    int i;
    luFactorization(n, copyMatrix(A, n), b);
    x0 = backwardSubstitution(A, b, n);
    r0 = vectorSubtraction(b, matvecprod(A, x0, n), n);
    r0Star = r0;
    M = cop(A, n);
    Minverse = inverseMatrix(M, n);
    Minverse[0] = A[0];
    //printMat(n, Minverse);
    for (i=1; i<5; i++) {
        z0 = matvecprod(Minverse, r0, n);
 
        z0Star = vecmatprod(r0Star, Minverse, n);
 
        rho1 = vecvecprod(z0, r0Star, n);
 
        if (rho1 == 0) return;
        if (i == 1) {
            p0 = z0;
            p0Star = z0Star;
        } else {
 
            beta = rho1/rho2;
 
            p0 = vectorAddition(z0, vectorNumMultiplication(p0, beta, n), n);
 
            p0Star = vectorAddition(z0Star, vectorNumMultiplication(p0Star, beta, n), n);
 
        }
 
 
        q = matvecprod(A, p0, n);
        qStar = matvecprod(transposeMatrix(A, n), p0Star, n);
 
        alphaK = rho1/vecvecprod(p0Star, q, n);
 
        x0 = vectorAddition(x0, vectorNumMultiplication(p0, alphaK, n), n);
 
        printf("szor u: \n");
        printVector(vectorNumMultiplication(q, alphaK, n), n);
        printf("\n");
        printf("r0 u: \n");
        printVector(r0, n);
        printf("\n");
        r0 = vectorSubtraction(r0, vectorNumMultiplication(q, alphaK, n), n);
 
        r0Star = vectorSubtraction(r0Star, vectorNumMultiplication(qStar, alphaK, n), n);
        rho2 = rho1;
    }
}
 
double* generateRandomMatrix(int n) {
    int i, j;
    double* resultMatrix = (double*)malloc(n*n*sizeof(double *));
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            resultMatrix[i*n+j] = -1+2*((double)rand())/RAND_MAX;
        }
    }
    return resultMatrix;
}
 
double* calculateBForMatrix(double *lowerMatrix, int n) {
    int i,j;
    double* resultVector = (double*)malloc(n*sizeof(double));
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            //matrix is multiplied by unit matrix
            resultVector[i] += lowerMatrix[i*n+j] * 1;
        }
    }
    return resultVector;
}
 
int main()
{
    int n = 2, i, j;
    double *A = (double*)malloc(n*n*sizeof(double *));
    double *b;
    A = generateRandomMatrix(n);
    b = calculateBForMatrix(A, n);
    BiCG(A, b, n);
 
    return 0;
}