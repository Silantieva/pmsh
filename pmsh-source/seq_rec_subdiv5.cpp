#include "defns5.h"

static int asorted_key_compare(const void *key1, const void *key2) ;

void seq_rec_subdiv(
int          lp,
int          rp,
int          l,
int          r,
int         *indx,
pSortedKey   keys)
{
    int   mp       ;
    int   n,m      ;
    int   axes = 0 ;

    if ( (rp - lp) > 0 ) {
      n = r - l + 1 ;
      inertia_xform(n,&(keys[l])) ;
      qsort(&(keys[l]),n,sizeof(SortedKey),asorted_key_compare) ;
      m  = ceil( ((double) l + r)/2.0 )  ;
      mp = ceil(( (double) lp + rp)/2.0) ;
      seq_rec_subdiv(lp,mp-1,l,m-1,indx,keys) ;
      seq_rec_subdiv(mp,rp,m,r,indx,keys) ;
    }
    else {
      indx[lp] = l ;
    }
}

static int asorted_key_compare(
const void *key1,
const void *key2)
{
    pSortedKey k1, k2 ;

    k1 = (pSortedKey) key1 ;
    k2 = (pSortedKey) key2 ;

    if ( k1->xyz[0] < k2->xyz[0] ) 
       return(-1) ; 
    else if ( k1->xyz[0] > k2->xyz[0]) 
       return(1) ; 
    else if (k1->xyz[1] < k2->xyz[1] )
       return(-1) ; 
    else if (k1->xyz[1] > k2->xyz[1] ) 
       return(1) ; 
    else if ( k1->xyz[2] < k2->xyz[2] )
       return(-1) ; 
    else 
       return(1) ; 

    /* return( ((k1->xyz[0] < k2->xyz[0]) ? -1 : 1) ); */
}

