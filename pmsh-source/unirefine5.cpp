#include "defns5.h"

extern std::map< int, FaceInfo > facemap ;
extern std::map< IntPair, int , IntPairCompare > edgemap ;  
extern int  vrtxcount ;

extern vector<Vertex> newvertices ; 
extern list<Face>     newfaces ; 


void unirefine(
nglib::Ng_Mesh  *submesh,
int numlevels)
{
   int  v[3] ;
   int  u[3] ; 
   list<Face> ::iterator li ; 
   std::map< IntPair, int >::iterator ei ;
   Face f ;
   double vxyz[3][3] ;
   double uxyz[3][3] ;
   int listsize ; 
   IntPair pr ;
   
   li=newfaces.begin() ;
   for(int l=0 ; l < numlevels ; l++) {
      listsize = newfaces.size() ; 
      //li=newfaces.begin() ;
      for(int t = 0 ; t <  listsize ; t++) { 
         for(int k=0 ; k < 3 ; k++) {
            v[k] = (*li).lsvrtx[k] ;
            Ng_GetPoint(submesh,v[k],vxyz[k]) ; 
         } 
         for(int k=0 ; k < 3 ; k++) {
            if ( v[k] < v[(k+1)%3]) {
              pr.x = v[k] ;
              pr.y = v[(k+1)%3] ; 
            }
            else {
              pr.y = v[k] ;
              pr.x = v[(k+1)%3] ;
            }

            //printf("(%d,%d) => %llu ", v[k],v[(k+1)%3],pr) ; 
            ei = edgemap.find(pr);
            if( ei == edgemap.end() ) {
               //printf("no\n") ;
               uxyz[k][0] = 0.5*(vxyz[k][0] + vxyz[(k+1)%3][0]) ;
               uxyz[k][1] = 0.5*(vxyz[k][1] + vxyz[(k+1)%3][1]) ;
               uxyz[k][2] = 0.5*(vxyz[k][2] + vxyz[(k+1)%3][2]) ;
               Ng_AddPoint(submesh,uxyz[k]);
               vrtxcount++ ; 
               u[k] = vrtxcount ;
               edgemap[pr] = vrtxcount ; 
            }
            else {
               u[k] = ei->second; 
               //printf("yes %d\n",u[k]) ; 
            }
         } 
         
         // outer triangles 
         for(int k=0 ; k < 3 ; k++) {
            f.lsvrtx[0] = v[k];
            f.lsvrtx[1] = u[k] ;
            f.lsvrtx[2] = u[(k+2)%3];
            newfaces.push_back(f) ; 
         }
         // interior triangle 
         f.lsvrtx[0] = u[0]; 
         f.lsvrtx[1] = u[1]; 
         f.lsvrtx[2] = u[2];
         newfaces.push_back(f) ;
         
         li=newfaces.erase(li) ;
      }
	printf("MY: refining level %d\n",l+1);
   }
   
   for(li=newfaces.begin() ; li != newfaces.end() ; ++li) {
      Ng_AddSurfaceElement(submesh,nglib::NG_TRIG,(*li).lsvrtx); 
   }
	printf("MY: refinement ended\n");
}
