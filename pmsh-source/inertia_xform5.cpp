#include "defns5.h"

void     nerror(char *)   ;
void     eig_jacobi(double **,double *,double **,int,int *) ;
void     free_dvector(double *) ;
void     eigsrt(double *,double **,int) ;
void     free_dmatrix(double **,int) ; 

void     inertia_xform(int,pSortedKey) ; 
void     xform_points(pSortedKey,int,double **) ; 
void     comp_inertia_matrix(pSortedKey,int,double **) ; 
void     free_matrix(double **,int) ; 


double  *dvector(int)     ;
double **matrix(int,int)  ;


/*
 *    Transforms the mesh to 'principal axis' coordinate system.
 */
void inertia_xform(
int        n,
pSortedKey keys)
{
   int      nrot ;
   double  *D,
          **M,
          **R ; 

   M = matrix(3,3) ;
   R = matrix(3,3) ;
   D = dvector(3)  ;
   
   /* 1. Compute mesh   centroid and the moment of inertia matrix M  */
   comp_inertia_matrix(keys,n,M) ; 
   
   /* 2. find the eigenvectors of the inertia matrix 
    *    i.e.  R^(-1)*M*R = D      where D = diag(lambda1,lambda2,lambda3)
    */
    eig_jacobi(M,D,R,3,&nrot) ;   /* by jacobi method                    */
    eigsrt(D,R,3) ;           /* sort eigenvalues in ascending order */
   
   /* 3. transform the coordinates to new coordinate system by:
    *      [x' y' z'] = [x y z]*R 
    */
   xform_points(keys,n,R) ;

   free_dvector(D) ; 
   free_matrix(M,3)  ;
   free_matrix(R,3)  ;
}

void comp_inertia_matrix(
pSortedKey keys,
int        n,
double   **M)
{
    int    i, j              ;
    double centr[4]          ; 
    double globM[6], locM[6] ;

    for(i=0 ; i < 4 ; i++) {
        centr[i] =  0.0 ;
    }
    for(i=0 ; i < n ; i++) {
      for(j=0 ; j < 3 ; j++)      
        centr[j] += keys[i].xyz[j] ; 
    }

    for(j=0 ; j < 3 ; j++)
      centr[j] /= n ; 
    
    for(i=0 ; i < n ; i++) {

      for(j=0 ; j < 3 ; j++)
        keys[i].xyz[j] -= centr[j] ;

      /* moments of inertia */
      M[1][1] += keys[i].xyz[1]*keys[i].xyz[1]+keys[i].xyz[2]*keys[i].xyz[2] ; 
      M[2][2] += keys[i].xyz[0]*keys[i].xyz[0]+keys[i].xyz[2]*keys[i].xyz[2] ; 
      M[3][3] += keys[i].xyz[0]*keys[i].xyz[0]+keys[i].xyz[1]*keys[i].xyz[1] ; 

      /* product of inertia terms */
      M[1][2] = M[2][1] += - keys[i].xyz[0]*keys[i].xyz[1] ; 
      M[1][3] = M[3][1] += - keys[i].xyz[0]*keys[i].xyz[2] ;
      M[2][3] = M[3][2] += - keys[i].xyz[1]*keys[i].xyz[2] ;

    }
}


void xform_points(
pSortedKey keys,
int       n,
double   **R)
{
    int i ; 
    int  j,k ;
    double xyz[3] ; 
    
    for(i=0 ; i < n ; i++) {
       for(j=0 ; j < 3 ; j++) {
         xyz[j] = 0.0 ; 
         for(k=0 ; k < 3 ; k++) 
           xyz[j] += keys[i].xyz[k]*R[k+1][j+1] ; 
       }
       for(j=0 ; j < 3 ; j++) 
          keys[i].xyz[j] = xyz[j] ; 
    }
}
