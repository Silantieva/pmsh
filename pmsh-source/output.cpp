#include "defns5.h"

extern std::map< int, FaceInfo > facemap ;
extern std::map< int, int > g2lvrtxmap ;

std::map< int, list<int> > coarsevrtx2adjprocmap ;

std::map< int, intmap * > proc2vrtxmap ;

std::map< int, int > locvrtx2globidmap ;

void compvrtxadj()
{
    map< int, FaceInfo >::iterator itf ;
    map< int, int >::iterator itv ;
    int fid ; 
    FaceInfo finfo ;
    int  pid, i ;
    map< int, list<int> >::iterator ivrtx ; 
    list<int>::iterator li ;

    for (itf=facemap.begin(); itf != facemap.end(); ++itf) {
       fid   = itf->first     ; 
       finfo = itf->second    ;  
       for(int j=0 ; j < 2 ; j++) {    
           pid   = finfo.procids[j] ;  
           for(int k=0 ; k < 3 ; k++) {
                itv = g2lvrtxmap.find(finfo.svrtx[j][k]);
                if( itv != g2lvrtxmap.end() ) { 
                   insadjlist(itv->first,pid) ; 
                }
           } 
       } 
    }
     
    
    for (ivrtx = coarsevrtx2adjprocmap.begin() ; 
          ivrtx != coarsevrtx2adjprocmap.end(); ++ivrtx) {
       printf("%d ==> ",ivrtx->first) ;
       for(li = (ivrtx->second).begin() ; li != (ivrtx->second).end() ; ++li) {
          printf(" %d ",(*li)) ;
       } 
       printf("\n") ; 
    }
    
}


void insadjlist(
int vrtx,
int pid)
{
     map< int, list<int> >::iterator ivrtx ; 
     list<int> pidlist ; 
     list<int>::iterator li ;
     
     // printf("%d %d\n",vrtx,pid) ; 
     ivrtx = coarsevrtx2adjprocmap.find(vrtx) ;
     if( ivrtx == coarsevrtx2adjprocmap.end() ) { 
        coarsevrtx2adjprocmap[vrtx] = pidlist ; 
        coarsevrtx2adjprocmap[vrtx].push_back(pid) ;   
     }
     else {
        for(li = (ivrtx->second).begin() ; li != (ivrtx->second).end() ; ++li) {
          if ( (*li) == pid) break  ; 
        }  
        if (  li == (ivrtx->second).end() ) {
           coarsevrtx2adjprocmap[vrtx].push_back(pid) ;  
        }   
     }
}






