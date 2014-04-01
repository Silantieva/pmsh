#include "defns5.h"

void com_length_from(
MPI_Comm comm,
int      num_s,
int      num_r,
int     *dest,
int     *src,
int     *length_s,
int     *length_r,
int      tag)
{
   int          i     ; 
   MPI_Status  *stat  ;
   MPI_Request *req   ;
   int          rc    ;

   if (num_s + num_r) {
      MYCALLOC(req,MPI_Request *,(num_s + num_r),sizeof(MPI_Request) ) ;
      MYCALLOC(stat,MPI_Status *,(num_s + num_r),sizeof(MPI_Status) ) ;
   }
   
   for(i=0 ; i < num_s ; i++) {
       rc = MPI_Isend(&(length_s[i]),1,MPI_INT,dest[i],tag,
                       comm,&(req[i])) ;
       if (rc !=  MPI_SUCCESS) { printf("Error in mpi\n") ; abort() ; }
   }

   for(i= 0 ; i < num_r ; i++) {
       rc=MPI_Irecv(&(length_r[i]),1,MPI_INT,MPI_ANY_SOURCE,tag,
                    comm,&(req[num_s+i])) ;
       if (rc !=  MPI_SUCCESS) { printf("Error in mpi\n") ; abort() ; }
   }

   if (num_s + num_r) {
        rc=MPI_Waitall(num_s + num_r,req,stat) ;
        if (rc !=  MPI_SUCCESS) { printf("Error in mpi\n") ; abort() ; }
   }

   for(i=0 ; i < num_r ; i++) {
      src[i] = stat[num_s+i].MPI_SOURCE ;
   }

   if (num_s + num_r) {
      free(req)  ; 
      free(stat) ; 
   }

   MPI_Barrier(comm) ;
}
