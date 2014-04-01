#include <iostream>
#include <fstream>

using namespace std;

#include "defns5.h" 

int main (int argc, char ** argv)
{

   nglib::Ng_Mesh *mesh;             // Define pointer to a new Netgen Mesh 
   Ng_STL_Geometry *stl_geom; // Define pointer to STL Geometry
   Ng_Result ng_res;          // Result of Netgen Operations
   int numprocs    ;          // number of processors 
   int  mypid      ;          // mypid 
   char geofilename[100], filename[100], filename2[100]  ;   // geometry file
   char *inpfile ;
   int    tp,ts,tv ; 
   int surf_p, surf_s, surf_v;
   int vol_p, vol_s, vol_v;
   int numlevels = 0 ; 

   int *geoid;

   double t1, t2, t3, t4, tpar_, tpar;
   
   MPI_Init(&argc,&argv) ;
   MPI_Comm_rank(MPI_COMM_WORLD,&mypid);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   mypid = atoi(argv[6]); 
   
   inpfile = argv[1] ;   
   sprintf(geofilename,"%s.geo",inpfile)  ;
   numprocs = atoi(argv[2]) ; 
   
   // initialize the Netgen Core library
   Ng_Init();

   // create the mesh structure
   mesh = Ng_NewMesh();

   Ng_CSG_Geometry *geom = new Ng_CSG_Geometry(); //loaded geometry

   nglib::Ng_Meshing_Parameters mp;
   mp.maxh = atof(argv[3]) ;
   mp.fineness = atof(argv[4]) ;
   //mp.maxh = 1e6;
   //mp.fineness = 1.0 ;

   MPI_Barrier(MPI_COMM_WORLD);
   t1 = MPI_Wtime();
   ng_res = Ng_CSG_GenerateMesh (geofilename, mesh, geom, mp, nglib::MESH_SURFACE,
                                 &surf_p, &surf_v, &surf_s, geoid); //loaded geometry to geom, generated volume mesh

   ng_res = Ng_CSG_GenerateMesh (geofilename, mesh, geom, mp, nglib::MESH_VOLUME,
                                 &vol_p, &vol_v, &vol_s, geoid); //loaded geometry to geom, generated volume mesh

   if (mypid == 0) {
     printf("MY: Global coarse surface mesh - No. of P/S: %d %d\n",surf_p,surf_s) ;
   }
    
   cout << "MY: Start Volume Meshing...." << endl;

   ng_res = Ng_CSG_GenerateVolumeMesh (mesh, mp, &vol_p, &vol_v, &vol_s, geoid);

  if(ng_res != NG_OK)  {
      cout << "MY: Error in Volume Meshing....Aborting!!" << endl;
        return 1;
   }
   // Barrier 
   //timer serial finish 

   sprintf(filename, "mesh.init");
   //My_Ng_Export_OF (mesh, geofilename, filename);

   if (mypid == 0) { 
     printf("MY: Global coarse volume mesh - No. of P/S/V: %d %d %d\n",vol_p,vol_s,vol_v) ; 
   }
  
   MPI_Barrier(MPI_COMM_WORLD);
   t2 = MPI_Wtime();
   
   nglib::Ng_Mesh *submesh ;
   submesh = Ng_NewMesh();
   
   if (1) {  
     pSortedKey keys ;
     int numkeys ;
     int indx[numprocs+1] ;
     compute_keys(mesh,&numkeys,&keys) ;
     seq_rec_subdiv(0,numprocs-1,0,numkeys-1,indx,keys) ;
     indx[numprocs] = numkeys ; 
     //printf("MY: Partitioning complete:  numprocs: %d\nMY:",numprocs) ; 
     //for(int k=0 ; k <= numprocs ; k++) printf(" %d ",indx[k]) ; 
     //printf("\n") ;
     cre_part_surf_mesh(mesh,numprocs,mypid,indx,keys,submesh) ;
   } 
   else {
      partwmetis(mesh,numprocs) ; 
      cre_part_surf_mesh(mesh,numprocs,mypid,submesh) ;
   }
   printf("MY: [%2d] Submeshes creation complete..\n",mypid) ;

   printf("MY: [%2d] Submeshes creation complete..\n",mypid) ;
   numlevels = atoi(argv[5]) ;
   unirefine(submesh,numlevels) ; 

   sprintf(filename,"%s.%d.vol",inpfile,mypid) ;
   //printf("MY: [%2d] Saving Mesh in VOL Format\n",mypid);
   //Ng_SaveMesh(submesh,filename);
   //exit(0);
   int psc[3], tpsc[3] ;
   psc[0]=psc[1]=psc[2]=0 ;
   tpsc[0]=tpsc[1]=tpsc[2]=0 ;
   Ng_DeleteMesh (mesh) ;  // NNN

   //sprintf(filename,"%s.%d.vol",inpfile,mypid) ;
   //printf("MY: [%2d] Saving Mesh in VOL Format\n",mypid);
   //Ng_SaveMesh(submesh,filename);

   //outgeomview(mypid,numprocs,submesh,inpfile) ;
   // tetgen call

#ifdef TETGEN
   tetgenio in, out;
   sprintf(filename,"%s.%d",inpfile,mypid);
   sprintf(filename2,"%s.%d.off",inpfile,mypid);
   in.load_off(filename2);
   tetrahedralize("pqYNEFO", &in, NULL, filename);
   // read number from output.off file to psc[0], psc[1], psc[2]
   // ???
#endif

#ifdef NETGEN
   mp.fineness = atof(argv[4]) ;
   mp.maxh = atof(argv[3]) ;
   mp.second_order = 0; 

   ng_res = Ng_CSG_GenerateVolumeMesh (submesh, mp, &psc[0], &psc[2], &psc[1], geoid);

   sprintf(filename,"%s.%d.vol",inpfile,mypid) ;
   printf("MY: [%2d] Saving Mesh in VOL Format\n",mypid);
   Ng_SaveMesh(submesh,filename);

#endif

   // Barrier
   //timer parallel finish 
   MPI_Barrier(MPI_COMM_WORLD);
   t3 = MPI_Wtime();

   //outgeomview(mypid,numprocs,submesh,inpfile) ;
   sprintf(filename, "mesh.%d",mypid);
   //My_Ng_Export_OF (submesh, geofilename, filename);

   printf("MY: [%2d] No. of P/S/V: %d %d %d\n",mypid,psc[0],psc[1],psc[2]) ;
   MPI_Reduce(psc,tpsc,3,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD) ; 
   if (mypid == 0) {
      printf("MY: Total No. of P/S/V: %d %d %d\n",tpsc[0],tpsc[1],tpsc[2]) ;
   } 

   if (mypid == 0) {
      printf("RESULT: %d %d %.3f %.3f %d %d %d %d %d %d %d %d %d %.2f %.2f\n",
                numprocs, numprocs, mp.maxh, mp.fineness, numlevels,
                surf_p, surf_s, vol_p, vol_s, vol_v, // after surf & after volume mesh
                tpsc[0], tpsc[1], tpsc[2],
                t2-t1, t3-t2); // after refinement
   }

   MPI_Finalize() ;
   exit(0) ;   
}

