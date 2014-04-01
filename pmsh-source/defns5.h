#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h> 
#include <vector>
#include <map>
#include <set>
#include <list>
#include <algorithm>
#include <nginterface.h>
#include "metis.h"

#ifndef DEFNS_H
#define DEFNS_H
 
namespace nglib {
#include <nglib.h>
}

using namespace nglib;
using namespace std ; 


/* memory allocation routines */
#define  MYMALLOC(v,type,size) {  v = (type) malloc(size) ;                  \
                                  if( v == NULL) {                           \
                                    printf("Memory alloc error:%lu\n",size); \
                                    abort();                                 \
                                  }                                          \
                                }


#define  MYREALLOC(v,type,size) {  v = (type) realloc(v,size) ;              \
                                  if( v == NULL) {                           \
                                    printf("Memory alloc error:%lu\n",size); \
                                    abort();                                 \
                                  }                                          \
                                }


#define MYCALLOC(v,type,n,size) {  v = (type) calloc(n,size) ;               \
                                  if( v == NULL) {                           \
                                    abort();                                 \
                                  }                                          \
                                }

typedef struct SortedKey {
 int       id       ;
 double    xyz[3]   ;
} SortedKey         ;

typedef SortedKey *pSortedKey ; 

typedef struct FaceInfo {
 int  procids[2] ;
 int  svrtx[2][3] ;
 short outw ; 
} FaceInfo   ;

typedef struct Barycentric {
  int    gvrtx[3] ;  // in ascending order of global vertex ids (starting with 1) 
  short  coord[3] ;
  int    newgid   ;  // global data 
} Barycentric   ;

typedef std::vector<Barycentric> BarycVector ;

typedef struct Baryvrtx {
  int    gvrtx[3] ;  // in ascending order of global vertex ids (starting with 1) 
} Barycvrtx   ;

typedef struct Face {
 int         lsvrtx[3] ;
 Barycentric barycv[3] ;
 short       outw ; 
 short       geoboundary ; 
} Face   ;

struct CompBarycentric {
    bool operator()( const Barycentric & first, const Barycentric & second) const {
       for(int i=0 ; i < 3 ; i++) {
         if (first.gvrtx[i] < second.gvrtx[i]) {
            return(true) ; 
        }
        else if (first.gvrtx[i] > second.gvrtx[i]) {
            return(false) ;    
        }
       } 
       for(int i=0 ; i < 3 ; i++) {
         if (first.coord[i] < second.coord[i]) {
            return(true) ; 
         }
         else if (first.coord[i] > second.coord[i]) {
            return(false) ;  
          }  
       }
       return(false) ;  
    }      
}  ;

struct CompBarycvrtx {
    bool operator()( const Barycvrtx & first, const Barycvrtx & second) const {
       for(int i=0 ; i < 3 ; i++) {
         if (first.gvrtx[i] < second.gvrtx[i]) {
            return(true) ; 
        }
        else if (first.gvrtx[i] > second.gvrtx[i]) {
            return(false) ;    
        }
       } 
       return(false) ;  
    }      
}  ;

typedef struct Vertex {
double  xyz[3] ;
} Vertex   ;

typedef struct IntPair {
   int  x ; 
   int  y ; 
} IntPair ; 

struct IntPairCompare {
    bool operator()( const IntPair & first, const IntPair & second) const {
      return first.x < second.x || (first.x == second.x && first.y < second.y) ;
    } 
} ;

typedef std::map<int, int> intmap ;

typedef struct FacePoints {
  int p[3];
} FacePoints;

struct FacePointCompare {
    bool operator()( const FacePoints & first, const FacePoints & second) const {
      return (first.p[0] < second.p[0]) ||
        (first.p[0] == second.p[0] && first.p[1] < second.p[1]) ||
        (first.p[0] == second.p[0] && first.p[1] == second.p[1] && first.p[2] < second.p[2]);
    } 
} ;

typedef std::map<FacePoints, int, FacePointCompare> surfMap_t;


void inertia_xform(int,pSortedKey) ;
void compute_keys(nglib::Ng_Mesh *mesh,int *num_keys,pSortedKey *keys) ;
void seq_rec_subdiv(int lp,int rp,int l,int r,int *indx,pSortedKey keys) ; 
void insertetra2map(nglib::Ng_Mesh  *mesh1,int fid,int procid,int *vrts) ; 
void cre_part_surf_mesh(nglib::Ng_Mesh *mesh1,int numprocs,int mypid,int *indx,pSortedKey keys,
                        nglib::Ng_Mesh *submesh) ;
void cre_part_surf_mesh(nglib::Ng_Mesh *mesh1,int numprocs,int mypid,nglib::Ng_Mesh *submesh) ;
void cre_geom_surf_mesh(nglib::Ng_Mesh *mesh1,nglib::Ng_Mesh *mesh2) ; 
void partfacecreate(nglib::Ng_Mesh  *mesh1,nglib::Ng_Mesh *submesh,int numprocs,int mypid) ; 
void outgeomview(nglib::Ng_Mesh  *mesh,float    *rgb,char     *filename) ; 
void outgeomview(int mypd,int numprocs, nglib::Ng_Mesh *mesh, char *file) ;
int  tetrafaceoutward(int *tverts,int *fverts) ; 
void outsmy(nglib::Ng_Mesh  *mesh,char     *filename) ;
void ncolors(int n,vector<float> & rgbvals) ;
void hsv2rgb(float *rgb, float *hsv) ;
void sort3int(int *v) ;
void sort2int(int *v) ; 
void verifysubmesh(nglib::Ng_Mesh *mesh,nglib::Ng_Mesh  *submesh, int mypid) ;
void unirefine(nglib::Ng_Mesh  *submesh,int numlevels) ;
void partwmetis(nglib::Ng_Mesh  *mesh,int nparts) ;
void compvrtxadj() ; 
void insadjlist(int vrtx,int pid) ;
void outbaryc(nglib::Ng_Mesh *submesh,int mypid)  ;
void barycmidpoint(Barycentric p1,Barycentric p2, Barycentric & res) ;
Barycentric initbarycv(int v) ;
void insbaryadjlist(Barycvrtx bvrtx,int pid) ; 
int com_sr_datatype(MPI_Comm comm,int num_s,int num_r,int *dest,int *src,int *s_length,
                    int *r_length,Barycentric **s_data,Barycentric **r_data,
                    MPI_Datatype  datatype,int  mypid) ; 
void com_barycoords(nglib::Ng_Mesh  *submesh,MPI_Comm comm,int numprocs,int mypid) ; 
void computeadj(int mypid) ; 
int com_num_recvs(MPI_Comm  comm,int numprocs,int mypid,int num_s,int *target) ;
void com_length_from(MPI_Comm comm,int num_s,int num_r,int *dest,int *src,
                     int *length_s,int *length_r,int  tag) ;
                     
                     
                    
#endif

