#include "defns5.h"

int com_sr_datatype(
MPI_Comm      comm,
int           num_s,
int           num_r,
int          *dest,
int          *src,
int          *s_length,
int          *r_length,
Barycentric  **s_data,
Barycentric  **r_data,
MPI_Datatype  datatype,
int           mypid)
{
    int          i ; 
    MPI_Status  *stat  ;
    MPI_Request *req   ;
    int          rc    ;

    if ( (num_s + num_r) > 0 ) {
      MYCALLOC(req,MPI_Request *,(num_s + num_r),sizeof(MPI_Request) ) ;
      MYCALLOC(stat,MPI_Status *,(num_s+num_r),sizeof(MPI_Status) ) ;
    }
    else {
      return(MPI_SUCCESS) ; 
    }

    for(i=0 ; i < num_s ; i++) {
         rc = MPI_Isend(s_data[i],s_length[i],datatype,dest[i],mypid,comm,
                      &(req[i]) ) ;
         if (rc !=  MPI_SUCCESS) { printf("Error in mpi\n") ; exit(1) ; }
    }

    for(i= 0 ; i < num_r ; i++) {
         rc=MPI_Irecv(r_data[i],r_length[i],datatype,src[i],src[i],comm,
                      &(req[num_s+i]) ) ;
        if (rc !=  MPI_SUCCESS) { printf("Error in mpi\n") ; exit(1) ; }
    }

    if (num_s + num_r) {
        rc=MPI_Waitall(num_s + num_r,req,stat) ;
        if (rc !=  MPI_SUCCESS) { printf("Error in mpi\n") ; exit(1) ; }
    }

    free(req) ; 
    free(stat) ; 
  
    return(rc) ; 
}
