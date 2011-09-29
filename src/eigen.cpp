
/*
 * This stuff comes from Ziheng Yang's PAML
 */

#include "eigen.h"
extern int NOISE;

Eigen::Eigen()
{
}

Eigen::~Eigen()
{
}

int Eigen::getpi_sqrt (double pi[], double pi_sqrt[], int n, int *npi0)
{
    int j;
    for (j=0,*npi0=0; j<n; j++)
        if (pi[j]) pi_sqrt[(*npi0)++]=sqrt(pi[j]);
    *npi0 = n - *npi0;
    return(0);
}

int Eigen::eigenQREV (double Q[], double pi[], double pi_sqrt[], int n, int npi0,
                      double Root[], double U[], double V[])
{
    /*
       This finds the eigen solution of the rate matrix Q for a time-reversible
       Markov process, using the algorithm for a real symmetric matrix.
       Rate matrix Q = S * diag{pi} = U * diag{Root} * V,
       where S is symmetrical, all elements of pi are positive, and U*V = I.
       pi_sqrt[n-npi0] has to be calculated before calling this routine.

       [U 0] [Q_0 0] [U^-1 0]    [Root  0]
       [0 I] [0   0] [0    I]  = [0     0]

       Ziheng Yang, 25 December 2001 (ref is CME/eigenQ.pdf)
    */
    int i,j, inew,jnew, nnew=n-npi0, status;

    /* store in U the symmetrical matrix S = sqrt(D) * Q * sqrt(-D) */
    if (pi_sqrt==NULL) { printf("\nError: pi_sqrt should be calculated before.\n"); exit(-1); }

    if (npi0==0) {
        for (i=0; i<n; i++)
            for (j=0,U[i*n+i] = Q[i*n+i]; j<i; j++)
                U[i*n+j] = U[j*n+i] = (Q[i*n+j] * pi_sqrt[i]/pi_sqrt[j]);
        status=eigenRealSym(U, n, Root, V);
        for (i=0;i<n;i++) for (j=0;j<n;j++)  V[i*n+j] = U[j*n+i] * pi_sqrt[j];
        for (i=0;i<n;i++) for (j=0;j<n;j++)  U[i*n+j] /= pi_sqrt[i];
    }
    else {
        for (i=0,inew=0; i<n; i++) {
            if (pi[i]) {
                for (j=0,jnew=0; j<i; j++)
                    if (pi[j]) {
                        U[inew*nnew+jnew] = U[jnew*nnew+inew]
                                            = Q[i*n+j] * pi_sqrt[inew]/pi_sqrt[jnew];
                        jnew++;
                    }
                U[inew*nnew+inew] = Q[i*n+i];
                inew++;
            }
        }
        status=eigenRealSym(U, nnew, Root, V);

        for (i=n-1,inew=nnew-1; i>=0; i--)   /* construct Root */
            Root[i] = (pi[i] ? Root[inew--] : 0);
        for (i=n-1,inew=nnew-1; i>=0; i--) {  /* construct V */
            if (pi[i]) {
                for (j=n-1,jnew=nnew-1; j>=0; j--)
                    if (pi[j]) {
                        V[i*n+j] = U[jnew*nnew+inew]*pi_sqrt[jnew];
                        jnew--;
                    }
                    else
                        V[i*n+j] = (i==j);
                inew--;
            }
            else
                for (j=0; j<n; j++)  V[i*n+j] = (i==j);
        }
        for (i=n-1,inew=nnew-1; i>=0; i--) {  /* construct U */
            if (pi[i]) {
                for (j=n-1,jnew=nnew-1;j>=0;j--)
                    if (pi[j]) {
                        U[i*n+j] = U[inew*nnew+jnew]/pi_sqrt[inew];
                        jnew--;
                    }
                    else
                        U[i*n+j] = (i==j);
                inew--;
            }
            else
                for (j=0;j<n;j++)
                    U[i*n+j] = (i==j);
        }
    }

    if (fabs(Root[0])>1e-10 && NOISE>0) printf("Root[0] = %.5e\n",Root[0]);
    Root[0]=0;
    return(status);
}


/* eigen solution for real symmetric matrix */


int Eigen::eigenRealSym(double A[], int n, double Root[], double work[])
{
    /* This finds the eigen solution of a real symmetrical matrix A[n*n].  In return,
       A has the right vectors and Root has the eigenvalues. work[n] is the working space.
       The matrix is first reduced to a tridiagonal matrix using HouseholderRealSym(),
       and then using the QL algorithm with implicit shifts.

       Adapted from routine tqli in Numerical Recipes in C, with reference to LAPACK
       Ziheng Yang, 23 May 2001
    */
    int status=0;
    HouseholderRealSym(A, n, Root, work);
    status=EigenTridagQLImplicit(Root, work, n, A);
    EigenSort(Root, A, n);

    return(status);
}


void Eigen::EigenSort(double d[], double U[], int n)
{
    /* this sorts the eigen values d[] and rearrange the (right) eigen vectors U[]
    */
    int k,j,i;
    double p;

    for (i=0;i<n-1;i++) {
        p=d[k=i];
        for (j=i+1;j<n;j++)
            if (d[j] >= p) p=d[k=j];
        if (k != i) {
            d[k]=d[i];
            d[i]=p;
            for (j=0;j<n;j++) {
                p=U[j*n+i];
                U[j*n+i]=U[j*n+k];
                U[j*n+k]=p;
            }
        }
    }
}

void Eigen::HouseholderRealSym(double a[], int n, double d[], double e[])
{
    /* This uses HouseholderRealSym transformation to reduce a real symmetrical matrix
       a[n*n] into a tridiagonal matrix represented by d and e.
       d[] is the diagonal (eigends), and e[] the off-diagonal.
    */
    int m,k,j,i;
    double scale,hh,h,g,f;

    for (i=n-1;i>=1;i--) {
        m=i-1;
        h=scale=0;
        if (m > 0) {
            for (k=0;k<=m;k++)
                scale += fabs(a[i*n+k]);
            if (scale == 0)
                e[i]=a[i*n+m];
            else {
                for (k=0;k<=m;k++) {
                    a[i*n+k] /= scale;
                    h += a[i*n+k]*a[i*n+k];
                }
                f=a[i*n+m];
                g=(f >= 0 ? -sqrt(h) : sqrt(h));
                e[i]=scale*g;
                h -= f*g;
                a[i*n+m]=f-g;
                f=0;
                for (j=0;j<=m;j++) {
                    a[j*n+i]=a[i*n+j]/h;
                    g=0;
                    for (k=0;k<=j;k++)
                        g += a[j*n+k]*a[i*n+k];
                    for (k=j+1;k<=m;k++)
                        g += a[k*n+j]*a[i*n+k];
                    e[j]=g/h;
                    f += e[j]*a[i*n+j];
                }
                hh=f/(h*2);
                for (j=0;j<=m;j++) {
                    f=a[i*n+j];
                    e[j]=g=e[j]-hh*f;
                    for (k=0;k<=j;k++)
                        a[j*n+k] -= (f*e[k]+g*a[i*n+k]);
                }
            }
        }
        else
            e[i]=a[i*n+m];
        d[i]=h;
    }
    d[0]=e[0]=0;

    /* Get eigenvectors */
    for (i=0;i<n;i++) {
        m=i-1;
        if (d[i]) {
            for (j=0;j<=m;j++) {
                g=0;
                for (k=0;k<=m;k++)
                    g += a[i*n+k]*a[k*n+j];
                for (k=0;k<=m;k++)
                    a[k*n+j] -= g*a[k*n+i];
            }
        }
        d[i]=a[i*n+i];
        a[i*n+i]=1;
        for (j=0;j<=m;j++) a[j*n+i]=a[i*n+j]=0;
    }
}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

int Eigen::EigenTridagQLImplicit(double d[], double e[], int n, double z[])
{
    /* This finds the eigen solution of a tridiagonal matrix represented by d and e.
       d[] is the diagonal (eigenvalues), e[] is the off-diagonal
       z[n*n]: as input should have the identity matrix to get the eigen solution of the
       tridiagonal matrix, or the output from HouseholderRealSym() to get the
       eigen solution to the original real symmetric matrix.
       z[n*n]: has the orthogonal matrix as output

       Adapted from routine tqli in Numerical Recipes in C, with reference to
       LAPACK fortran code.
       Ziheng Yang, May 2001
    */
    int m,j,iter,niter=30, status=0, i,k;
    double s,r,p,g,f,dd,c,b, aa,bb;

    for (i=1;i<n;i++) e[i-1]=e[i];  e[n-1]=0;
    for (j=0;j<n;j++) {
        iter=0;
        do {
            for (m=j;m<n-1;m++) {
                dd=fabs(d[m])+fabs(d[m+1]);
                if (fabs(e[m])+dd == dd) break;  /* ??? */
            }
            if (m != j) {
                if (iter++ == niter) {
                    status=-1;
                    break;
                }
                g=(d[j+1]-d[j])/(2*e[j]);

                /* r=pythag(g,1); */

                if ((aa=fabs(g))>1)  r=aa*sqrt(1+1/(g*g));
                else                r=sqrt(1+g*g);

                g=d[m]-d[j]+e[j]/(g+SIGN(r,g));
                s=c=1;
                p=0;
                for (i=m-1;i>=j;i--) {
                    f=s*e[i];
                    b=c*e[i];

                    /*  r=pythag(f,g);  */
                    aa=fabs(f); bb=fabs(g);
                    if (aa>bb)       { bb/=aa;  r=aa*sqrt(1+bb*bb); }
                    else if (bb==0)             r=0;
                else            { aa/=bb;  r=bb*sqrt(1+aa*aa); }

                    e[i+1]=r;
                    if (r == 0) {
                        d[i+1] -= p;
                        e[m]=0;
                        break;
                    }
                    s=f/r;
                    c=g/r;
                    g=d[i+1]-p;
                    r=(d[i]-g)*s+2*c*b;
                    d[i+1]=g+(p=s*r);
                    g=c*r-b;
                    for (k=0;k<n;k++) {
                        f=z[k*n+i+1];
                        z[k*n+i+1]=s*z[k*n+i]+c*f;
                        z[k*n+i]=c*z[k*n+i]-s*f;
                    }
                }
                if (r == 0 && i >= j) continue;
                d[j]-=p; e[j]=g; e[m]=0;
            }
        } while (m != j);
    }
    return(status);
}

#undef SIGN

void Eigen::computePMatrix(int n, double* pMat, double* U, double* V, double* Root, double time) {

    // Create the P(T) matrix
    double *P = pMat;
    for (int i=0;i<n*n;i++) {
        *(P++) = 0;
    }

    double *pdV,*pdU;
    double e1,e2;

    for (int k=0;k<n;k++) {
        P = pMat;
        pdU = &U[k];  // Set P for pointer arithmatic

        e1 = exp(time*Root[k]);

        for (int i=0;i<n;i++) {

            e2 = *pdU*e1;
            pdV = &V[k*n];
            pdU += n;

            for (int j=0;j<n;j++) {
                *P++ += (e2 * *(pdV++));
            }
        }
    }
}
