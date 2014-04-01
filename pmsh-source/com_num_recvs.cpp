#include "defns5.h"

int com_num_recvs(
MPI_Comm  comm,
int       numprocs,
int       mypid,
int       num_s,
int      *target)
{
    int *sends_arry ;   
    int *recvs_arry ;   
    int  i,ret      ;

    MYCALLOC(sends_arry,int *,numprocs,sizeof(int)) ; 
    MYCALLOC(recvs_arry,int *,numprocs,sizeof(int)) ; 

    for(i=0 ; i < num_s ; i++) {
       sends_arry[target[i]] = 1 ;
    }

    MPI_Allreduce(sends_arry,recvs_arry,numprocs,
                  MPI_INT,MPI_SUM,comm) ;

    ret = recvs_arry[mypid] ; 

    free(sends_arry) ; 
    free(recvs_arry) ; 

    return(ret) ; 
}