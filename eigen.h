
/*
computePMatrix() function is modified from the code of Simon Whelan,
all the rest comes from Ziheng Yang's software package paml3.14.
*/

#ifndef EIGEN_H
#define EIGEN_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>

class Eigen
{
public:
    Eigen();
    ~Eigen();

    int getpi_sqrt (double pi[], double pi_sqrt[], int n, int *npi0);
    int eigenQREV (double Q[], double pi[], double pi_sqrt[], int n, int npi0, double Root[], double U[], double V[]);
    int eigenRealSym(double A[], int n, double Root[], double work[]);
    void computePMatrix(int n, double* pMat, double* U, double* V, double* Root, double time);

    void HouseholderRealSym(double a[], int n, double d[], double e[]);
    int EigenTridagQLImplicit(double d[], double e[], int n, double z[]);
    void EigenSort(double d[], double U[], int n);
    
};

#endif
