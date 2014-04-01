#include <iostream>
#include <fstream>

using namespace std;

#include "defns5.h" 

int numlevels = 0 ; 
int maxbarycoord = 0 ; 

int main (int argc, char ** argv)
{

   Ng_Mesh *mesh;             // Define pointer to a new Netgen Mesh 
   Ng_STL_Geometry *stl_geom; // Define pointer to STL Geometry
   Ng_Result ng_res;          // Result of Netgen Operations
   int numprocs    ;          // number of processors 
   int  mypid      ;          // mypid 
   char filename[100]  ;   // geometry file
   char *inpfile ;
   int    tp,ts,tv ; 

   
   MPI_Init(&argc,&argv) ;
   MPI_Comm_rank(MPI_COMM_WORLD,&mypid);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   // mypid = 103 ; 
   
   inpfile = argv[1] ;   
   sprintf(filename,"%s.stl",inpfile)  ;
   numprocs = atoi(argv[2]) ; 
   
   // initialize the Netgen Core library
   Ng_Init();

   // create the mesh structure
   mesh = Ng_NewMesh();

   // Read in the STL File
   stl_geom = Ng_STL_LoadGeometry(filename);
   if(!stl_geom) {
      cout << "MY: Error reading in STL File: " << filename << endl;
        return 1;
   }
   cout << "MY: Successfully loaded STL File: " << filename << endl;


   // Set the Meshing Parameters to be used
   Ng_Meshing_Parameters mp;
   mp.maxh = atof(argv[3]) ;
   mp.fineness = atof(argv[4]) ;
   //mp.maxh = 1e6;
   //mp.fineness = 1.0 ;
   mp.second_order = 0;
 
   cout << "MY: Initialise the STL Geometry structure...." << endl;
   // timer başlangıç
   ng_res = Ng_STL_InitSTLGeometry(stl_geom);
   if(ng_res != NG_OK) {
      cout << "Error Initialising the STL Geometry....Aborting!!" << endl;
         return 1;
   }
   
   cout << "MY: Start Edge Meshing...." << endl;
   ng_res = Ng_STL_MakeEdges(stl_geom, mesh, &mp);
   if(ng_res != NG_OK) {
      cout << "Error in Edge Meshing....Aborting!!" << endl;
         return 1;
   }

   cout << "MY: Start Surface Meshing...." << endl;
   ng_res = Ng_STL_GenerateSurfaceMesh(stl_geom, mesh, &mp);
   if(ng_res != NG_OK)  {
      cout << "Error in Surface Meshing....Aborting!!" << endl;
         return 1;
   }
   tp = Ng_GetNP(mesh) ; 
   ts = Ng_GetNSE(mesh);
   if (mypid == 0) {
     printf("MY: Global coarse surface mesh - No. of P/S: %d %d\n",tp,ts) ;
   }
    
   cout << "MY: Start Volume Meshing...." << endl;
   ng_res = Ng_GenerateVolumeMesh (mesh, &mp);
   if(ng_res != NG_OK)  {
      cout << "MY: Error in Volume Meshing....Aborting!!" << endl;
        return 1;
   }
   // Barrier 
   //timer serial finish 
   tp = Ng_GetNP(mesh) ; 
   ts = Ng_GetNSE(mesh);
   tv = Ng_GetNE(mesh) ; 
   if (mypid == 0) { 
     printf("MY: Global coarse volume mesh - No. of P/S/V: %d %d %d\n",tp,ts,tv) ; 
   }
   
   numlevels = atoi(argv[5]) ;
   maxbarycoord =  1  ; 
   maxbarycoord =  1  << numlevels ; 
   
   Ng_Mesh *submesh ;
   submesh = Ng_NewMesh();
   
   if (0) {  
     pSortedKey keys ;
     int numkeys ;
     int indx[numprocs+1] ;
     compute_keys(mesh,&numkeys,&keys) ;
     seq_rec_subdiv(0,numprocs-1,0,numkeys-1,indx,keys) ;
     indx[numprocs] = numkeys ; 
     printf("MY: Partitioning complete:  numprocs: %d\nMY:",numprocs) ; 
     for(int k=0 ; k <= numprocs ; k++) printf(" %d ",indx[k]) ; 
     printf("\n") ;
     cre_part_surf_mesh(mesh,numprocs,mypid,indx,keys,submesh) ;
   } 
   else {
      partwmetis(mesh,numprocs) ; 
      cre_part_surf_mesh(mesh,numprocs,mypid,submesh) ;
   }
  
   printf("MY: [%2d] Submeshes creation complete..\n",mypid) ;

   unirefine(submesh,numlevels) ; 
   computeadj(mypid) ;  
   //MPI_Finalize() ;
   //exit(0) ;

   int psc[3], tpsc[3] ;
   psc[0]=psc[1]=psc[2]=0 ;
   tpsc[0]=tpsc[1]=tpsc[2]=0 ;
   mp.maxh = atof(argv[3]) ;
   mp.fineness = atof(argv[4]) ;
   mp.second_order = 0; 
   Ng_DeleteMesh (mesh) ;  // NNN
   ng_res = Ng_GenerateVolumeMesh(submesh, &mp); 
   // Barrier
   //timer parallel finish 
   psc[0] = Ng_GetNP(submesh) ; 
   psc[1] = Ng_GetNSE(submesh);
   psc[2] = Ng_GetNE(submesh) ;  
   printf("MY: [%2d] No. of P/S/V: %d %d %d\n",mypid,psc[0],psc[1],psc[2]) ;
   MPI_Reduce(psc,tpsc,3,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD) ; 
   if (mypid == 0) {
      printf("MY: Total No. of P/S/V: %d %d %d\n",tpsc[0],tpsc[1],tpsc[2]) ;
   } 

   
   sprintf(filename,"%s.%d.vol",inpfile,mypid) ; 
   //printf("MY: [%2d] Saving Mesh in VOL Format\n",mypid);  
   //Ng_SaveMesh(submesh,filename);
 
   outgeomview(mypid,numprocs,submesh,inpfile) ;
   verifysubmesh(mesh,submesh,mypid) ;
   // outbaryc(submesh,mypid) ; 
   com_barycoords(submesh,MPI_COMM_WORLD,numprocs,mypid) ; 
   
   MPI_Finalize() ;
   exit(0) ;   
}
