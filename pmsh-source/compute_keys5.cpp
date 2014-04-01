#include "defns5.h"

void compute_keys(
nglib::Ng_Mesh        *mesh,
int            *num_keys,
pSortedKey     *keys)
{
   int     i,k,n     ;
   double  xyz[3]    ; 
   int     entity[4] ;

   *keys = NULL ; 
   *num_keys = Ng_GetNE(mesh) ; 
   if ( *num_keys > 0) {
     MYMALLOC(*keys,pSortedKey,(*num_keys)*sizeof(SortedKey)) ;
   }
   
   for(i=0 ; i < *num_keys ; i++) {
      Ng_GetVolumeElement(mesh,i+1,entity) ;
      (*keys)[i].xyz[0] = 0.0 ;
      (*keys)[i].xyz[1] = 0.0 ;
      (*keys)[i].xyz[2] = 0.0 ;
      for(int k = 0 ; k < 4 ; k++) { 
          Ng_GetPoint (mesh,entity[k],xyz) ; 
          (*keys)[i].xyz[0] += xyz[0] ; 
          (*keys)[i].xyz[1] += xyz[1] ; 
          (*keys)[i].xyz[2] += xyz[2] ;
      } 
      (*keys)[i].xyz[0] += xyz[0]/4.0 ; 
      (*keys)[i].xyz[1] += xyz[1]/4.0 ; 
      (*keys)[i].xyz[2] += xyz[2]/4.0 ;
      (*keys)[i].id = i+1 ;
   }
   
   
   //for(i=0 ; i < *num_keys ; i++) {
   //  printf("%3d\n",(*keys)[i].id) ;  
   //}

}
