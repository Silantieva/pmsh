#include "defns5.h"

extern idx_t *eptr   ; 
extern idx_t *eind   ;
extern idx_t *edest  ;  
extern idx_t *ndest  ;

void partwmetis(
nglib::Ng_Mesh  *mesh,
int nparts)
{
   idx_t  ne          ;  // number of elements in the mesh 
   idx_t  nn          ;  // number of nodes in the mesh
   idx_t  ncommon = 3 ; 
   idx_t  objval      ; 
   int    vrts[4]     ;
   int    rc          ; // return code 
   int    i           ;
   int    elemno      ; 
   idx_t options[METIS_NOPTIONS];
   
   METIS_SetDefaultOptions(options); 
   options[METIS_OPTION_CONTIG] = 1 ; 
   nn = Ng_GetNP(mesh) ; 
   ne = Ng_GetNE(mesh) ;
   printf("MY: (nn,ne) (%d %d)\n",nn,ne) ; 
   MYMALLOC(eptr,idx_t *,((ne+1)*sizeof(idx_t))) ;
   MYMALLOC(edest,idx_t *,(ne*sizeof(idx_t))) ;
   MYMALLOC(eind,idx_t *, (4*ne*sizeof(idx_t)) ) ;
   MYMALLOC(ndest,idx_t *,(nn*sizeof(idx_t))) ;
   
   eptr[0] = 0 ; 
   i = 0 ; 
   for(elemno=1 ; elemno <= ne ; elemno++) {
      Ng_GetVolumeElement(mesh,elemno,vrts) ;
      eind[i++] = vrts[0]-1 ;
      eind[i++] = vrts[1]-1 ;
      eind[i++] = vrts[2]-1 ;
      eind[i++] = vrts[3]-1 ;
      eptr[elemno] = eptr[elemno-1] + 4 ;  
   }
   
   /****
   printf("MY: (nn,ne) (%d %d)\n",nn,ne) ;
   for(elemno=0 ; elemno < ne ; elemno++) {
       for(int j = eptr[elemno] ; j < eptr[elemno+1] ; j++) {
          printf(" %d ",eind[j]) ; 
       }
       printf("\n") ; 
   } 
   exit(0) ; 
   **/
   
   rc = METIS_PartMeshDual(&ne,&nn,eptr,eind,NULL,NULL,&ncommon,&nparts,NULL,options, 
        &objval,edest,ndest) ;
        
    printf("MY: metis objval = %d\n",objval) ;  
}
