#include "defns5.h"

extern std::map< int, FaceInfo > facemap ;

void verifysubmesh(
nglib::Ng_Mesh   *mesh,
nglib::Ng_Mesh   *submesh,
int        mypid)
{
   map< int, FaceInfo >::iterator itf ;
   int facecount ; 
   int fid ; 
   FaceInfo finfo ;
   
   
   facecount = 0 ; 

   for (itf=facemap.begin(); itf != facemap.end(); ++itf) {
      fid   = itf->first  ; 
      finfo = itf->second ;
      if ( finfo.procids[0] == mypid ) {
         facecount++ ;  
      }   
      if ( finfo.procids[1] == mypid ) {
        facecount++ ;  
      }      
   } 

   printf("MY: verify (mesh,submesh %d faces): %d %d\n",mypid,facecount,
                                                    Ng_GetNSE(submesh)); 
}
