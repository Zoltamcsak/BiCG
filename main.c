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

double* vectorSubtraction(double* vec1, double* vec2, int n) {
    int i;
    for (i=0; i<n; i++) {
        vec1[i] -= vec2[2];
    }
    return vec1;
}

double* vectorAddition(double* vec1, double* vec2, int n) {
    int i;
    for (i=0; i<n; i++) {
        vec1[i] += vec2[i];
    }
    return vec1;
}

double vectorVectorMultiplication(double* vec1, double* vec2, int n) {
    int i;
    double result = 0;
    for (i=0; i<n; i++) {
        result += vec1[i] * vec2[i];
    }
}

double* vectorNumMultiplication(double* vec, double num, int n) {
    int i;
    double* resultVector = (double*)malloc(n*sizeof(double));
    for (i=0; i<n; i++) {
        resultVector[i] = resultVector[i] * num;
    }
    return resultVector;
}

double** transposeMatrix(double** mat, int n) {
    double** transpose = (double**)malloc(n*sizeof(double *));
    int i, j;
    for (i=0; i<n; i++) {
        transpose[i] = (double*)malloc(n*sizeof(double));
        for (j=0; j<n; j++) {
            transpose[i][j] = mat[j][i];
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

void luFactorization(int n, double** a, double* b) {
    int k,i,j;
    double* c = (double*)malloc(n*sizeof(double));
    double** m = (double**)malloc(n*sizeof(double *));
    for (i=0; i<n; i++) {
        m[i] = (double*)malloc(n*sizeof(double));
        m[i][i] = 1;
    }
    for (k=0; k<n-1; k++){
        if (a[k][k] == 0) {
            break;
        }

        for (i=k+1; i<n; i++) {
          m[i][k] = (a[i][k]/a[k][k]);
          b[i] = b[i] - (m[i][k]*b[k]);
            for (j=k+1; j<n; j++) {
                a[i][j] = a[i][j] - (m[i][k]*a[k][j]);
            }
        }
    }
}

double* backwardSubstitution(double **upperMatrix, double *y, int n) {
    int i, j;
    double* resultVector = (double*)malloc(n*sizeof(double));
    for (i=n-1; i>=0; i--) {
        resultVector[i] = y[i];
        for (j=i+1; j<n; j++) {
            resultVector[i] -= upperMatrix[i][j] * resultVector[j];
        }
        resultVector[i] /= upperMatrix[i][i];
    }
    return resultVector;
}

double **copyMatrix(double** A, int n) {
    int i,j;
    double** copyRes = (double**)malloc(n*sizeof(double *));
    for (i=0; i<n; i++) {
        copyRes[i] = (double*)malloc(n*sizeof(double));
        for (j=0; j<n; j++) {
            copyRes[i][j] = A[i][j];
        }
    }
    return copyRes;
}

double *cop(double **A, int n) {
    int i, j, x = 0;
    double *result = (double*)malloc(n*sizeof(double));
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            result[x] = A[i][j];
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

void printMat(int n, double** a) {
    int i,j;
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            printf("%f ", a[i][j]);
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

void BiCG(double **A, double *b, int n) {
    double *r0, *x0, *r0Star, *z0, *z0Star, alphaK, *xKplus1, *xKplus1Star, *rKplus1, *rKplus1Star, *pKplus1, *pKplus1Star;
    double rho1, rho2, beta, *p0, *p0Star, *q, *qStar;
    double *M, **Minverse;
    int i;
    luFactorization(n, copyMatrix(A, n), b);
    x0 = backwardSubstitution(A, b, n);
    r0 = vectorSubtraction(b, matrixVectorMultiplication(A, x0, n), n);
    r0Star = r0;
    M = cop(A, n);
    for (i=1; i<n; i++) {
        Minverse = backcop(inverseMatrix(M, n), n);
        z0 = matrixVectorMultiplication(Minverse, r0, n);
        z0Star =vectorMatrixMultiplication(r0Star, Minverse, n);
        rho1 = vectorVectorMultiplication(z0, r0Star, n);
        if (rho1 == 0) return;
        if (i=1) {
            p0 = z0;
            p0Star = z0Star;
        } else {
            beta = rho1/rho2;
            p0 = vectorAddition(z0, vectorNumMultiplication(p0, beta, n), n);
            p0Star = vectorAddition(z0Star, vectorNumMultiplication(p0Star, beta, n), n);
        }
        q = matrixVectorMultiplication(A, p0, n);
        qStar = matrixVectorMultiplication(transposeMatrix(A, n), p0Star, n);
        alphaK = rho1/vectorVectorMultiplication(p0Star, q, n);
        x0 = vectorAddition(x0, vectorNumMultiplication(p0, alphaK, n), n);
        r0 = vectorSubtraction(r0, vectorNumMultiplication(q, alphaK, n), n);
        r0Star = vectorSubtraction(r0Star, vectorNumMultiplication(qStar, alphaK, n), n);
    }
}

double** generateRandomMatrix(int n) {
    int i, j;
    double** resultMatrix = (double**)malloc(n*sizeof(double *));
    for (i=0; i<n; i++) {
        resultMatrix[i] = (double*)malloc(n*sizeof(double));
        for (j=0; j<n; j++) {
            resultMatrix[i][j] = -1+2*((double)rand())/RAND_MAX;
        }
    }
    return resultMatrix;
}

double* calculateBForMatrix(double **lowerMatrix, int n) {
    int i,j;
    double* resultVector = (double*)malloc(n*sizeof(double));
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            //matrix is multiplied by unit matrix
            resultVector[i] += lowerMatrix[i][j] * 1;
        }
    }
    return resultVector;
}

int main()
{
    int n = 2, i;
    double **A = (double**)malloc(n*sizeof(double *));
    double *b;
    for(i=0; i<n; i++) {
        A[i] = (double*)malloc(n*sizeof(double));
    }
    A = generateRandomMatrix(n);
    b = calculateBForMatrix(A, n);
    BiCG(A, b, n);
    return 0;
}
