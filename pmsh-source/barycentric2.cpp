#include "defns5.h"

extern int maxbarycoord ;
extern std::map< int, FaceInfo > facemap ;
extern std::map< int, int > g2lvrtxmap ;
extern std::map< Barycentric, int, CompBarycentric > baryc2locvrtxmap ;

std::map< Barycvrtx, list<int>, CompBarycvrtx > barycvrtx2adjprocsmap ;
std::map< int, BarycVector * > pidmap ;
std::map< int, int > srcmap ; // new
int newglobalnocounter = 0 ; 
int *newgid ; 

MPI_Datatype  barycdatatype ;

Barycentric initbarycv(
int v)
{
   Barycentric p ;
   
   p.gvrtx[0] = 0 ;  
   p.gvrtx[1] = 0 ;
   p.gvrtx[2] = v ;
   p.coord[0] = 0 ; 
   p.coord[1] = 0 ; 
   p.coord[2] = maxbarycoord ;  
   return(p) ; 
}

void barycmidpoint(
Barycentric p1,
Barycentric p2,
Barycentric & res)
{
    set<int> s ; 
    set< int >::iterator is ;
    int k ;  
    short w1, w2 ; 
    
    /*
    printf("== %d(%hd) %d(%hd) %d(%hd) \n",p1.gvrtx[0],p1.coord[0],
                                           p1.gvrtx[1],p1.coord[1],
                                           p1.gvrtx[2],p1.coord[2]) ; 
    printf("== %d(%hd) %d(%hd) %d(%hd) \n",p2.gvrtx[0],p2.coord[0],
                                           p2.gvrtx[1],p2.coord[1],
                                           p2.gvrtx[2],p2.coord[2]) ;
    */ 
    
    for(int i=0 ; i < 3 ; i++) {
      if (p1.coord[i]) s.insert(p1.gvrtx[i]) ; 
      if (p2.coord[i]) s.insert(p2.gvrtx[i]) ;
    }
    if ( s.size() > 3)  {
      printf("Error in barycentric coordinates\n") ; 
    }
    
    k = 0 ; 
    for(is=s.begin(); is != s.end(); is++) {
      res.gvrtx[k] = *is ;
      k++ ; 
    }
    while (k < 3)  {
      res.gvrtx[k] = 0 ; 
      k++ ;  
    }
    sort3int(res.gvrtx) ; 
    
    w1 = w2 = 0 ; 
    for(int i=0 ; i < 3 ; i++) {
        for(int j=0 ; j < 3 ; j++) { 
          if ( res.gvrtx[i] == p1.gvrtx[j] ) {
             w1 = p1.coord[j] ; 
             break ; 
          }
        }
        for(int j=0 ; j < 3 ; j++) { 
          if ( res.gvrtx[i] == p2.gvrtx[j] ) {
             w2 = p2.coord[j] ; 
             break ; 
          }
        }
        res.coord[i] = (w1 + w2) / 2 ; 
    }
}

void outbaryc(
nglib::Ng_Mesh  *submesh,
int       mypid)
{
    map< Barycentric, int, CompBarycentric >::iterator itc ;
    Barycentric    c      ; 
    int         locid  ; 
    double      xyz[3] ;
    double      cxyz[3] ;
    double      ccxyz[3][3] ;
    int         i, j ; 

    for (itc = baryc2locvrtxmap.begin() ; itc != baryc2locvrtxmap.end() ; ++itc) {
       c      = itc->first     ; 
       locid  = itc->second    ;  
       Ng_GetPoint(submesh,locid,xyz) ; 
       for(i=0 ; i < 3 ; i++) 
          for(j=0 ; j < 3 ; j++)  ccxyz[i][j] = 0.0 ; 
       if (c.gvrtx[0] > 0) Ng_GetPoint(submesh,g2lvrtxmap[c.gvrtx[0]],ccxyz[0]) ;
       if (c.gvrtx[1] > 0) Ng_GetPoint(submesh,g2lvrtxmap[c.gvrtx[1]],ccxyz[1]) ;
       if (c.gvrtx[2] > 0) Ng_GetPoint(submesh,g2lvrtxmap[c.gvrtx[2]],ccxyz[2]) ;
       cxyz[0] = c.coord[0]*ccxyz[0][0] + c.coord[1]*ccxyz[1][0] + c.coord[2]*ccxyz[2][0];
       cxyz[1] = c.coord[0]*ccxyz[0][1] + c.coord[1]*ccxyz[1][1] + c.coord[2]*ccxyz[2][1];
       cxyz[2] = c.coord[0]*ccxyz[0][2] + c.coord[1]*ccxyz[1][2] + c.coord[2]*ccxyz[2][2];
       cxyz[0] = cxyz[0] / maxbarycoord ; 
       cxyz[1] = cxyz[1] / maxbarycoord ;
       cxyz[2] = cxyz[2] / maxbarycoord ;
       printf("%3d(%d) %3d(%d) %3d(%d) %lf %lf %lf  [%2d] %lf %lf %lf\n",
               c.gvrtx[0],c.coord[0],c.gvrtx[1],c.coord[1],c.gvrtx[2],c.coord[2],
               xyz[0],xyz[1],xyz[2],mypid,
               cxyz[0],cxyz[1],cxyz[2]) ; 
    }
}

void computeadj(
int mypid)
{
    map< int, FaceInfo >::iterator itf ;
    map< int, int >::iterator itv ;
    int fid ; 
    FaceInfo finfo ;
    int  pid, i ;
    map< Barycvrtx, list<int>, CompBarycvrtx >::iterator ib ; 
    list<int>::iterator li ;
    int count ; 
    Barycvrtx bvrtx ; 
    Barycvrtx countbvrtx ; 

    for (itf=facemap.begin(); itf != facemap.end(); ++itf) {
       fid   = itf->first     ; 
       finfo = itf->second    ;  
       if (  (finfo.procids[0] != mypid) && (finfo.procids[1] == -1) ) continue ;
       if (  (finfo.procids[1] != mypid) && (finfo.procids[0] == -1) ) continue ; 
       for(int j=0 ; j < 2 ; j++) {    
           pid   = finfo.procids[j] ; 
           if ( (pid == mypid) ||  (pid == -1) ) continue ;  
           count = 0 ;           
           for(int k=0 ; k < 3 ; k++) {
                countbvrtx.gvrtx[k] = 0 ; 
                itv = g2lvrtxmap.find(finfo.svrtx[j][k]);
                if( itv != g2lvrtxmap.end() ) { 
                   bvrtx.gvrtx[0]      = 0 ; 
                   bvrtx.gvrtx[1]      = 0 ; 
                   bvrtx.gvrtx[2]      = finfo.svrtx[j][k] ;
                   countbvrtx.gvrtx[k] = finfo.svrtx[j][k] ; 
                   count++ ; 
                   insbaryadjlist(bvrtx,pid) ; 
                }
           }
           sort3int(countbvrtx.gvrtx) ; 
           if (count == 2) {  // insert edge adjacency
              insbaryadjlist(countbvrtx,pid) ;
           }
           if (count == 3) {  // insert edge and face adjacencies
               bvrtx.gvrtx[0] = 0 ;
               bvrtx.gvrtx[1] = countbvrtx.gvrtx[0] ;
               bvrtx.gvrtx[2] = countbvrtx.gvrtx[1] ; 
               insbaryadjlist(bvrtx,pid) ;
               bvrtx.gvrtx[0] = 0 ;
               bvrtx.gvrtx[1] = countbvrtx.gvrtx[0] ;
               bvrtx.gvrtx[2] = countbvrtx.gvrtx[2] ;                
               insbaryadjlist(bvrtx,pid) ;
               bvrtx.gvrtx[0] = 0 ;
               bvrtx.gvrtx[1] = countbvrtx.gvrtx[1] ;
               bvrtx.gvrtx[2] = countbvrtx.gvrtx[2] ; 
               insbaryadjlist(bvrtx,pid) ;
               insbaryadjlist(countbvrtx,pid) ;   // face adjacency 
           }
       } 
    } 
    
    /**
    for (ib = barycvrtx2adjprocsmap.begin() ; 
          ib != barycvrtx2adjprocsmap.end(); ++ib) {
       printf("%d>(%2d %2d %2d) ==> ",mypid,ib->first.gvrtx[0],ib->first.gvrtx[1],ib->first.gvrtx[2]) ;
       for(li = (ib->second).begin() ; li != (ib->second).end() ; ++li) {
          printf(" %d ",(*li)) ;
       } 
       printf("\n") ; 
    } 
    **/
        
}


void insbaryadjlist(
Barycvrtx bvrtx,
int       pid)
{
     map< Barycvrtx, list<int>, CompBarycvrtx >::iterator ivrtx ; 
     list<int> pidlist ; 
     list<int>::iterator li ;
     
     // printf("%d %d\n",vrtx,pid) ; 
     ivrtx = barycvrtx2adjprocsmap.find(bvrtx) ;
     if( ivrtx == barycvrtx2adjprocsmap.end() ) { 
        barycvrtx2adjprocsmap[bvrtx] = pidlist ; 
        barycvrtx2adjprocsmap[bvrtx].push_back(pid) ;   
     }
     else {
        for(li = (ivrtx->second).begin() ; li != (ivrtx->second).end() ; ++li) {
          if ( (*li) == pid) break  ; 
        }  
        if (  li == (ivrtx->second).end() ) {
           barycvrtx2adjprocsmap[bvrtx].push_back(pid) ;  
        }   
     }
}


void com_barycoords(
nglib::Ng_Mesh  *submesh,
MPI_Comm comm,
int numprocs,
int mypid)
{
    int count = 3            ; 
    int blocklens[3]         ;
    MPI_Aint addrs[3]        ; 
    MPI_Datatype mpitypes[3] ; 
    MPI_Datatype mpibaryctype ;
    int num_s, num_r         ;    // number of sends and receives 
    int *dest, *src ;             //  destination and source processors 
    int  *s_length, *r_length ;   // sent/received data size 
    Barycentric **s_data, **r_data ;  // sent/received data 
    Barycentric brcy ;              // struct used to get struct member 
    map< Barycvrtx, list<int> >::iterator ibc ;
    list<int>::iterator li ;
    std::map< Barycentric, int, CompBarycentric >::iterator ib ;
    Baryvrtx bvrtx ; 
    map< int, BarycVector * >::iterator ipdt ;
    std::map< int, int >::iterator srcit ; // new
    int i,j,indxowner, ownerpid, locid ;
    vector<int> holders ; 
    int numverts ; 
    int *globoffsets, *globoffsetsm1 ; 
    
    // initialize new global ids array
    numverts = Ng_GetNP(submesh) ;
    MYCALLOC(newgid,int *,(numverts+1),sizeof(int)) ; // ids start with 1, initialized to 0
    
    num_s = 0 ;
    num_r = 0 ; 
    
    // construct MPI Datatype 
    blocklens[0] = 3 ;
    blocklens[1] = 3 ; 
    blocklens[2] = 1 ;
    mpitypes[0]  = MPI_INT ;
    mpitypes[1]  = MPI_SHORT ; 
    mpitypes[2]  = MPI_INT ;
    MPI_Address(&brcy.gvrtx,addrs) ;
    MPI_Address(&brcy.coord,addrs+1) ;
    MPI_Address(&brcy.newgid,addrs+2) ;
    addrs[1] = addrs[1] - addrs[0] ;
    addrs[2] = addrs[2] - addrs[0] ;
    addrs[0] = (MPI_Aint) 0 ; 
    MPI_Type_struct(count,blocklens,addrs,mpitypes,&mpibaryctype) ;
    MPI_Type_commit(&mpibaryctype) ; 
    
    newglobalnocounter = 0 ;   // initialize global number counter 
    // Loop goes over the  geometric and partition boundary vertices in the new mesh
    for (ib=baryc2locvrtxmap.begin(); ib != baryc2locvrtxmap.end() ; ++ib) {
      brcy  = ib->first ; 
      locid = ib->second ; 
      bvrtx.gvrtx[0] = brcy.gvrtx[0] ; 
      bvrtx.gvrtx[1] = brcy.gvrtx[1] ; 
      bvrtx.gvrtx[2] = brcy.gvrtx[2] ;
      ibc = barycvrtx2adjprocsmap.find(bvrtx) ; 
      if ( ibc == barycvrtx2adjprocsmap.end() ) { 
         // Not : these are entities on the geometric boundary, therefore skip them  
         //printf("%d>Error: barycentric coordinate not found in proc map: (%d,%d,%d)\n",
         //        mypid,bvrtx.gvrtx[0],bvrtx.gvrtx[1],bvrtx.gvrtx[2]) ;
         continue ; 
      }
      
      // compute owners of shared vertices (held by multiple processors - called holders) 
      holders.clear() ; 
      holders.push_back(mypid) ; // i am also a holder 
      for(li = (ibc->second).begin() ; li != (ibc->second).end() ; ++li) {
        holders.push_back((*li)) ;  // insert other holders 
      }
      std::sort(holders.begin(), holders.end() ) ;
      indxowner = (brcy.gvrtx[0] + brcy.gvrtx[1] + brcy.gvrtx[2] + 
              brcy.coord[0] + brcy.coord[1] + brcy.coord[2] ) % (ibc->second).size() ;
      ownerpid = holders[indxowner] ; 
      
      // if I am the owner, append it to the message to be sent 
      if (ownerpid == mypid) {
         newglobalnocounter++ ; 
         newgid[locid] = newglobalnocounter ;
         brcy.newgid    = newglobalnocounter ; 
         for(li = (ibc->second).begin() ; li != (ibc->second).end() ; ++li) {
           ipdt = pidmap.find((*li)) ; 
           if( ipdt == pidmap.end() ) {  
             pidmap[(*li)] = new BarycVector() ; 
           }
           pidmap[(*li)]->push_back(brcy)  ;
         }  
       }
       else {  // i am not the owner of this partition boundary vertex 
          newgid[locid] = -1 ; 
		  // new 
          srcit = srcmap.find(ownerpid) ; 
          if( srcit == srcmap.end() ) {  
             srcmap[ownerpid] = 1 ; 
          }
          else {
            srcit->second = srcit->second + 1 ;
          } 
		  // new
       }  
    }
    
    // now do global numbering of the interior vertices 
    // i.e. non-partition boundary vertices
    //  +ve ones already numbered, -1 ones owned by other processors  
    for(locid=1 ; locid <= numverts ; locid++) {
        if (newgid[locid] == 0) {
           newglobalnocounter++ ;
           newgid[locid] = newglobalnocounter ; 
        }
    }
    
    // compute pre-scan of all newglobalnocounter in array globoffsets
    MYCALLOC(globoffsetsm1,int *,(numprocs+1),sizeof(int)) ;
    globoffsets = globoffsetsm1 + 1 ; 
    globoffsets[-1] = 0 ; 
    MPI_Allgather(&newglobalnocounter,1,MPI_INT,globoffsets,1,MPI_INT,comm);
    for(i=0 ; i < numprocs ; i++) globoffsets[i] += globoffsets[i-1] ; 
    
    // add offsets to global numbers 
	// new
    for(locid=1 ; locid <= numverts ; locid++) {
        if (newgid[locid] != -1) {
           newgid[locid] += globoffsets[mypid-1] ; 
        }  
    }
	// new
    
    num_s = pidmap.size() ; 
    if (num_s > 0) MYCALLOC(s_length,int *,num_s,sizeof(int)) ; 
    if (num_s > 0) MYCALLOC(dest,int *,num_s,sizeof(int)) ;
    if (num_s > 0) MYCALLOC(s_data,Barycentric **,num_s,sizeof(Barycentric *)) ;
    i = 0 ; 
    for (ipdt=pidmap.begin(); ipdt != pidmap.end() ; ++ipdt) {
       s_length[i] = ipdt->second->size() ; 
       dest[i] = ipdt->first ;
       s_data[i] = &((*ipdt->second)[0]) ;
       i++ ; 
    }
   
    num_r = srcmap.size() ;  
    if (num_r > 0) MYCALLOC(r_length,int *,num_r,sizeof(int)) ;
    if (num_r > 0) MYCALLOC(src,int *,num_r,sizeof(int)) ;
    if (num_r > 0) MYCALLOC(r_data,Barycentric **,num_r,sizeof(Barycentric *)) ;
    
    // compute lengths of messages (no. of items) that will be sent     
	// new
    for (srcit=srcmap.begin(), i=0 ; srcit != srcmap.end() ; ++srcit, ++i) {
      src[i]      = srcit->first ; 
      r_length[i] = srcit->second ;
    }
	// new
    
    for(i=0 ; i < num_r ; i++) {
       MYCALLOC(r_data[i],Barycentric *,r_length[i],sizeof(Barycentric)) ;
    }
     
    // send/recv all the messages                       
    com_sr_datatype(comm,num_s,num_r,dest,src,s_length,r_length,
                    s_data,r_data,mpibaryctype,mypid) ;
    
    for(i=0 ; i < num_r ; i++) {
       for(j=0 ; j < r_length[i] ; j++) {
         locid = baryc2locvrtxmap[r_data[i][j]] ; 
         if ( newgid[locid] != -1) printf("%d> Error: remote global id %d\n",mypid,locid);  
         newgid[locid] = -(r_data[i][j].newgid + globoffsets[src[i]-1]) ;  
       } 
    }
    
    // Free allocated  memory 
    for(i=0 ; i < num_s ; i++) {
         free(s_data[i]) ;
    }
    for(i=0 ; i < num_r ; i++) {
         free(r_data[i]) ;
    }
    if (num_s > 0) {
      free(s_length) ;
      free(s_data) ;
      free(dest) ;
    } 
    if (num_r > 0) {
      free(r_length) ;
      free(r_data) ;
      free(src) ;
    } 
    free(globoffsetsm1) ; 
    MPI_Type_free(&mpibaryctype) ;
    
    /**  fcomment the following - it is for debugging purposes                
    printf("Num recvs: %d\n",num_r) ; 
    for(i=0 ; i < num_r ; i++) {
       printf("r_length[%d] %d\n",i,r_length[i]) ; 
       for(j=0 ; j < r_length[i] ; j++) {
       printf("========%3d(%d) %3d(%d) %3d(%d) [%d->%d]\n",
               r_data[i][j].gvrtx[0],
               r_data[i][j].coord[0],
               r_data[i][j].gvrtx[1],
               r_data[i][j].coord[1],
               r_data[i][j].gvrtx[2],
               r_data[i][j].coord[2],
               src[i],mypid) ; 
       }
    }
    **/
    
    /* comment the following - it is for debugging purposes */ 
    /*
    for(i=0 ; i < numprocs ; i++) {
       if (i == mypid) {
         for(locid = 1 ; locid <= numverts ; locid++) {
            printf("GlobID|%d> (%d,%d) \n",mypid,locid,newgid[locid]) ; 
         }
       }
       MPI_Barrier(comm) ;
    }
    */
    
 
}
