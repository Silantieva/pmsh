#include "defns5.h"

extern std::map< int, FaceInfo > facemap ;
extern std::map< int, int > g2lvrtxmap ; 
extern int  vrtxcount ;
extern vector<Vertex> newvertices ; 
extern list<Face>     newfaces ; 

extern int filetype;
extern char geofilename[100];
extern Ng_STL_Geometry *stl_geom; // Define pointer to STL Geometry

std::map< Barycentric, int, CompBarycentric > baryc2locvrtxmap ;
std::map< IntPair, int , IntPairCompare > edgemap ;

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
   Barycentric midvrtxbaryc ; 
   Barycentric vbarycv[3] ;
   Barycentric ubarycv[3] ;
   std::map< Barycentric, int, CompBarycentric >::iterator bi ;
   vector<double> projectVertexList;
   vector<int> projectBoundaryList;
   
   li=newfaces.begin() ;
   for(int l=0 ; l < numlevels ; l++) {
      listsize = newfaces.size() ; 
      for(int t = 0 ; t <  listsize ; t++) { 
         for(int k=0 ; k < 3 ; k++) {
            v[k] = (*li).lsvrtx[k] ;
            vbarycv[k] = (*li).barycv[k] ;
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
            barycmidpoint(vbarycv[k],vbarycv[(k+1)%3],midvrtxbaryc) ;
            ei = edgemap.find(pr);
            bi = baryc2locvrtxmap.find(midvrtxbaryc);
            if( bi == baryc2locvrtxmap.end() ) {
            //if( ei == edgemap.end() ) {
               uxyz[k][0] = 0.5*(vxyz[k][0] + vxyz[(k+1)%3][0]) ;
               uxyz[k][1] = 0.5*(vxyz[k][1] + vxyz[(k+1)%3][1]) ;
               uxyz[k][2] = 0.5*(vxyz[k][2] + vxyz[(k+1)%3][2]) ;
               for (int ctr = 0; ctr < 3; ctr++) {
                 projectVertexList.push_back(uxyz[k][ctr]);
                }
              projectBoundaryList.push_back((*li).geoboundary);
               //Ng_AddPoint(submesh,uxyz[k]);
               vrtxcount++ ; 
               u[k] = vrtxcount ;
               edgemap[pr] = vrtxcount ;
               ubarycv[k]  = midvrtxbaryc ; 
               baryc2locvrtxmap[midvrtxbaryc] = vrtxcount ; 	
            }
            else {
               u[k] = ei->second ; 
               ubarycv[k]  = midvrtxbaryc ; 
            }
         } 
         
         // outer triangles 
         for(int k=0 ; k < 3 ; k++) {
            f.lsvrtx[0] = v[k];
            f.barycv[0] = vbarycv[k];
            f.lsvrtx[1] = u[k] ;
            f.barycv[1] = ubarycv[k] ;
            f.lsvrtx[2] = u[(k+2)%3];
            f.barycv[2] = ubarycv[(k+2)%3];
            f.geoboundary = (*li).geoboundary;
            newfaces.push_back(f) ; 
         }
         // interior triangle 
         f.lsvrtx[0] = u[0];
         f.barycv[0] = ubarycv[0]; 
         f.lsvrtx[1] = u[1];
         f.barycv[1] = ubarycv[1]; 
         f.lsvrtx[2] = u[2];
         f.barycv[2] = ubarycv[2];
         f.geoboundary = (*li).geoboundary;
         newfaces.push_back(f) ;
         
         li=newfaces.erase(li) ;
      }
   }
   
  if (filetype == 1) { 
    // prjoects points in projectVertexList to geometric id's in projectBoundaryList 
    // STL IMPLEMENTATION
    My_Ng_Project_Add_Points(submesh, projectVertexList, projectBoundaryList, stl_geom);
  }
  else if (filetype == 2) { 
    // prjoects points in projectVertexList to geometric id's in projectBoundaryList 
    // using geofilename, adds them to submesh
    // CSG IMPLEMENTATION
    My_Ng_Project_Add_Points(submesh, projectVertexList, projectBoundaryList, geofilename);
  }
   
   for(li=newfaces.begin() ; li != newfaces.end() ; ++li) {
      Ng_AddSurfaceElement(submesh,nglib::NG_TRIG,(*li).lsvrtx); 
   }
}


