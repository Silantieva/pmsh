#include "defns5.h"

extern std::map< Barycentric, int, CompBarycentric > baryc2locvrtxmap ; 
std::map< int, FaceInfo > facemap ;
std::map< int, int > g2lvrtxmap ;
bool netgenmeshupdate = true; 

vector<Vertex> newvertices ; 
list<Face>     newfaces ; 

int  vrtxcount ; 

int tetraor[4][3] = { {1,2,3}, {0,3,2}, {0,1,3}, {0,2,1} } ; 

surfMap_t surfMap;

idx_t *eptr        ; 
idx_t *eind        ;
idx_t *edest       ;  
idx_t *ndest       ;

void getSurfPoints(
nglib::Ng_Mesh *mesh) {
  int nse = Ng_GetNSE(mesh);
  int *surfpoints = new int[3];
  for (int i = 0; i < nse; i++) {
    Ng_GetSurfaceElement(mesh, i+1, surfpoints);
    FacePoints fp;
    sort(surfpoints, surfpoints + 3);
    fp.p[0] = surfpoints[0];
    fp.p[1] = surfpoints[1];
    fp.p[2] = surfpoints[2];
    surfMap[fp] = GetBoundaryID(mesh, i + 1);
  }
  /*
  int i = 1;
  for (mymap::const_iterator it = surfMap.begin(); it != surfMap.end(); ++it) {
    printf("SURF: %d: %d %d %d\n",i++, it->first.x, it->first.y, it->second);
  }
  */
  
}

int findSurfElem(
nglib::Ng_Mesh *mesh,
int *elem)
{
  int nse = Ng_GetNSE(mesh);
  //printf("---- SEARCHING FOR %d %d %d\n",elem[0], elem[1], elem[2]);
  sort(elem, elem + 3);
  FacePoints fp;
  fp.p[0] = elem[0];
  fp.p[1] = elem[1];
  fp.p[2] = elem[2];
  surfMap_t::const_iterator it = surfMap.find(fp);
  if (it == surfMap.end())
    return -1;
  else
    return it->second;
}

void  cre_geom_surf_mesh(
nglib::Ng_Mesh        *mesh1,
nglib::Ng_Mesh        *mesh2)
{
   
   int     np = Ng_GetNP(mesh1) ;
   int     nse = Ng_GetNSE(mesh1) ;  
   double  xyz[3] ; 
   int     entity[4] ; 
   
   // printf("%d %d\n",np,nse) ; 
   for(int i=1 ; i <= np ; i++) {
     Ng_GetPoint (mesh1,i,xyz) ;
     Ng_AddPoint (mesh2,xyz); 
     // printf("%lf %lf %lf\n",xyz[0],xyz[1],xyz[2]) ; 
   }
   
   for(int i=1 ; i <= nse ; i++) {
      Ng_GetSurfaceElement(mesh1,i,entity) ;
      Ng_AddSurfaceElement(mesh2,nglib::NG_TRIG,entity);
      // printf("%d %d %d\n",entity[0],entity[1],entity[2]) ;
   }
}

void cre_part_surf_mesh(
nglib::Ng_Mesh    *mesh1,
int         numprocs,
int         mypid,
nglib::Ng_Mesh    *submesh)
{
    int fids[4] ; 
    int vrts[4] ;
    int orient[4] ;
    int p,i,k,elemno  ; 
    FaceInfo finfo ; 
    int ne ; 
    
    getSurfPoints(mesh1);
    int nse = Ng_GetNSE(mesh1);
//  for (int i = 0; i < nse; i++) {
//    printf("SURF: %d: %d %d %d\n",i+1, surfpoints[i][0], surfpoints[i][1], surfpoints[i][2]); 
//  }
    
    ne = Ng_GetNE(mesh1) ;
    for(elemno=1 ; elemno <= ne ; elemno++) { 
      p = edest[elemno-1] ; 
         if (netgenmeshupdate) printf("MY: starting update\n") ; 
         My_Ng_GetElement_Faces(mesh1,elemno,fids,orient,netgenmeshupdate);
         //Ng_GetVolumeElement(mesh1, elemno, fids);
         //printf("YYY got element faces of %d : %d %d %d %d\n", elemno, fids[0], fids[1], fids[2], fids[3]);
         netgenmeshupdate = false; 
         //if (netgenmeshupdate) printf("MY: finishing update\n") ; 
         Ng_GetVolumeElement(mesh1,elemno,vrts) ;
         //printf("YYY got element points of %d : %d %d %d %d\n", elemno, vrts[0], vrts[1], vrts[2], vrts[3]);
         for(k=0 ; k < 4 ; k++) {
            insertetra2map(mesh1,fids[k],p,vrts) ;  
         }
    } 
    partfacecreate(mesh1,submesh,numprocs,mypid) ;
}


void cre_part_surf_mesh(
nglib::Ng_Mesh    *mesh1,
int         numprocs,
int         mypid,
int         *indx,
pSortedKey  keys,
nglib::Ng_Mesh    *submesh)
{
    int fids[4] ; 
    int vrts[4] ;
    int orient[4] ;
    int p,i,k,elemno  ; 
    FaceInfo finfo ; 
    
    for(p=0 ; p < numprocs ; p++) {
       for(i=indx[p] ; i < indx[p+1] ; i++) { 
          elemno = keys[i].id ;
          if (netgenmeshupdate) printf("MY: starting update\n") ; 
          My_Ng_GetElement_Faces(mesh1,elemno,fids,orient,netgenmeshupdate);
//        printf("
          netgenmeshupdate = false; 
          if (netgenmeshupdate) printf("MY: finishing update\n") ;
          //printf("QQ:inserttetra2map: %d %d : %d %d %d %d\n",p,elemno,fids[0],fids[1],fids[2],fids[3]) ; 
          Ng_GetVolumeElement(mesh1,elemno,vrts) ;
          for(k=0 ; k < 4 ; k++) {
             insertetra2map(mesh1,fids[k],p,vrts) ;  
          }
       } 
        
    } 
    partfacecreate(mesh1,submesh,numprocs,mypid) ;
}

int tetrafaceoutward(
int *tverts,
int *fverts)
{
    int i,j,s,t  ;
    int sumt ;
    
/*
    printf("TETRAOUT: ");
    for (i=0;i<3;i++) 
      printf("fverts[%d]: %d, ",i,fverts[i]);
    for (i=0;i<3;i++) 
      printf("tverts[%d]: %d, ",i,tverts[i]);
    printf("\n");
*/          
    
    for(i=0 ; i < 2 ; i++) {
       for(s=0 ; s < 3 ; s++) {
         if ( tverts[i] ==  fverts[s] ) break ;
       }
       if (s < 3) break ; 
    }
    if (s == 3) {
      printf("MY: error in outward face computation (1).\n") ;
      exit(1) ; 
    }
    sumt = 0 ; 
    for(int i=0 ; i < 4 ; i++) {
      t = 0 ; 
      for(int j=0 ; j < 3 ; j++) {
        if ( tverts[tetraor[i][j]] ==  fverts[(s+j)%3] ) t++ ; 
      } 
      if (t == 3) sumt++ ;      
    } 
    if (sumt == 1) {
      return(1) ; 
    }
    else if (sumt > 1) {
      printf("MY: error in outward face computation (2).\n") ;
      exit(1) ;       
    }
    return(0) ;   
}


void insertetra2map(
nglib::Ng_Mesh  *mesh,
int fid,
int procid,
int *tverts)
{
    FaceInfo finfo ;
    std::map< int, FaceInfo >::iterator it ;
    int fverts[3] ;
    int outw ; 
    
               
    My_Ng_GetFace_Vertices(mesh,fid,fverts);
    it= facemap.find(fid);
    outw = tetrafaceoutward(tverts,fverts) ; 
    if( it == facemap.end() ) {  
        // std::cout << "C: " << it->second << "\n";
        finfo.outw = outw ;
        finfo.procids[0] = procid ;
        finfo.procids[1] = -1 ;
        finfo.svrtx[0][0] = fverts[0] ; 
        if (! outw) {
          finfo.svrtx[0][1] =  fverts[1] ;
          finfo.svrtx[0][2] =  fverts[2] ;        
        }
        else {
         finfo.svrtx[0][1] =  fverts[2]  ;  
         finfo.svrtx[0][2] =  fverts[1]  ;         
        }
        facemap[fid] = finfo ;
        
        //printf("XXX fid %d, svrtx[0][0-1-2] is %d %d %d\n", fid, finfo.svrtx[0][0], finfo.svrtx[0][1], finfo.svrtx[0][2]);
    }
    else if ( (it->second).procids[0] == procid) { //internal face 
        facemap.erase(it) ; 
    } 
    else {  // partition boundary face 
        (it->second).procids[1]  = procid ; 
        (it->second).svrtx[1][0] = fverts[0] ; 
        if (! outw) {
           (it->second).svrtx[1][1] = fverts[1] ;
           (it->second).svrtx[1][2] = fverts[2] ;
        }
        else {
            (it->second).svrtx[1][1] = fverts[2] ;
            (it->second).svrtx[1][2] = fverts[1] ;     
        }                 
    }
}


void partfacecreate(
nglib::Ng_Mesh  *mesh1,
nglib::Ng_Mesh  *submesh,
int      numprocs,
int       mypid)
{
    map< int, FaceInfo >::iterator itf ;
    map< int, int > ::iterator itv ;
    int fid ; 
    FaceInfo finfo ;
    Face  f ; 
    int  lsvrtx[3] ;
    int  pid, i ; 
    double xyz[3] ; 
    Barycentric barycv ;   
    
    vrtxcount = 0 ;
    
    // add vertices to mesh 
    for (itf=facemap.begin(); itf != facemap.end(); ++itf) {
        fid   = itf->first     ; 
        finfo = itf->second    ;  
        // insert face vertices into verts map 
        for(int j=0 ; j < 2 ; j++) {    
           pid   = finfo.procids[j] ; 
           if ( pid == mypid ) {    
             for(int k=0 ; k < 3 ; k++) {
                itv = g2lvrtxmap.find(finfo.svrtx[j][k]);
                if( itv == g2lvrtxmap.end() ) { 
                   vrtxcount++ ; 
                   g2lvrtxmap[finfo.svrtx[j][k]] = vrtxcount ;
                   barycv = initbarycv(finfo.svrtx[j][k]) ;
                   baryc2locvrtxmap[barycv] = vrtxcount ;
                   Ng_GetPoint(mesh1,finfo.svrtx[j][k],xyz) ;
                  
                  Ng_AddPoint(submesh,xyz); 
                }
             }
           } 
        } 
    }

    // add surface elements to mesh 
    for (itf=facemap.begin(); itf != facemap.end(); ++itf) {
        fid   = itf->first  ; 
        finfo = itf->second ;
        for(int j=0 ; j < 2 ; j++) {
          pid = finfo.procids[j] ;
          if ( pid == mypid ) { 
            f.outw = finfo.outw ;
            //printf("GETTING BOUNDARY ID OF %d, j %d pid %d\n", fid,j,pid);
            //int boundaryID = GetBoundaryID(mesh1, fid);
            //printf("BOUNDARY ID OF %d is %d\n", fid, boundaryID);
            f.geoboundary = (finfo.procids[1] == -1) ? 1 : 0 ; 
            //f.geoboundary = (finfo.procids[(j+1)%2] != -1) ?  1 : 0 ; 
            f.lsvrtx[0] = g2lvrtxmap[finfo.svrtx[j][0]] ;
            f.lsvrtx[1] = g2lvrtxmap[finfo.svrtx[j][1]] ;
            f.lsvrtx[2] = g2lvrtxmap[finfo.svrtx[j][2]] ;
            f.barycv[0] = initbarycv(finfo.svrtx[j][0]) ;
            f.barycv[1] = initbarycv(finfo.svrtx[j][1]) ;
            f.barycv[2] = initbarycv(finfo.svrtx[j][2]) ;
            //printf("MY: addsurfel: %d : %d %d %d\n",pid,f.lsvrtx[0],f.lsvrtx[1],f.lsvrtx[2]) ;
            //if (f.geoboundary == 1) {
              int geo = findSurfElem(mesh1, finfo.svrtx[j]);
              f.geoboundary = geo;
              //printf("BOUNDARY ID OF %d, SETTING GEO TO %d, elmts %d %d %d, finfp.procids[0]: %d finfo.procids[1]: %d\n", fid, f.geoboundary, finfo.svrtx[j][0], finfo.svrtx[j][1], finfo.svrtx[j][2], finfo.procids[0], finfo.procids[1]);
            //}
            newfaces.push_back(f) ; 
          }
        }        
    }
}


void sort3int(
int *v)
{
   if (v[0] > v[1]) { v[0] ^= v[1] ; v[1] ^= v[0] ; v[0] ^= v[1] ; } 
   if (v[1] > v[2]) { v[1] ^= v[2] ; v[2] ^= v[1] ; v[1] ^= v[2] ; }
   if (v[0] > v[1]) { v[0] ^= v[1] ; v[1] ^= v[0] ; v[0] ^= v[1] ; }
}

void sort2int(
int *v)
{
   if (v[1] > v[0]) { v[0] ^= v[1] ; v[1] ^= v[0] ; v[0] ^= v[1] ; } 
}
