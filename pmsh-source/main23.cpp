#include <iostream>
#include <fstream>

using namespace std;

#include "defns5.h" 
#include <string.h> /* you need to include string to use strcmp*/

int numlevels = 0 ;
int maxbarycoord = 0 ;

Ng_STL_Geometry *stl_geom; // Define pointer to STL Geometry

char geofilename[100];
int filetype = 0; //1 = stl, 2 = geo, 0 = error
/*
char * getFileExt(int fileType) {
  char ret[3] = "";
  if (fileType == 1)
     strcpy(ret, "stl");
  if (fileType == 2)
    strcpy(ret, "geo");
}
*/

void printhelp() {
	cerr << "Invalid syntax. The parameters are: " << endl << 
		"	--stl: use stl type geometry as input. " << endl <<
        	"	--geo: use csg type geometry as input. (either geo or stl is required) " << endl <<
		"	--numprocs: number of processes (required). " << endl <<
		" 	--fineness: fineness level in netgen (required). " << endl <<
		"	--maxh: mamxh level in netgen (required). " << endl << 
		"	--numlevels: number of refinement levels in the submesh (required). " << endl << 
		"	--outvol: if set, netgen files will be printed for each submesh. " << endl << 
		"	--outof: if set, openfoam output will be printed. " << endl <<
		"	--outgeom: if set, geomview output will be printed. " << endl << 
		"	--help: this help will be printed. " << endl << endl <<
		"e.g. ./pmsh_netgen --stl onera-m6.stl --numcores 2048 --fineness 0.025 --maxh 0.1 --outgeom" << endl;
}

int main (int argc, char ** argv)
{

   nglib::Ng_Mesh *mesh;             // Define pointer to a new Netgen Mesh 
   Ng_Result ng_res;          // Result of Netgen Operations
   int  mypid      ;          // mypid 
   char filename[100], filename2[100]  ;   // geometry file
   char *inpfile ;
   int    tp,ts,tv ; 
   int surf_p, surf_s, surf_v;
   int vol_p, vol_s, vol_v;
  int numlevels = -1; // -1 = error
  double fineness = -1; // -1 = error
  double maxh = -1; // -1 = error
  bool outvol = false;
  bool outgeom = false;
  bool outof = false;
  int numprocs = 0;          // number of processors, 0 = error
  bool error = false;

  if (argc <= 1) {
    printhelp();
    return 1;
  }

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--stl") == 0) {
      filetype = 1;
      sprintf(geofilename,"%s",argv[++i])  ;
    } else  if (strcmp(argv[i], "--geo") == 0) {
      filetype = 2;
      sprintf(geofilename,"%s",argv[++i])  ;
    } else  if (strcmp(argv[i], "--numprocs") == 0) { 
      numprocs = atoi(argv[++i]);
    } else  if (strcmp(argv[i], "--fineness") == 0) { 
      fineness = atof(argv[++i]);
    } else  if (strcmp(argv[i], "--maxh") == 0) { 
      maxh = atof(argv[++i]);
    } else  if (strcmp(argv[i], "--numlevels") == 0) { 
      numlevels = atoi(argv[++i]);
    } else  if (strcmp(argv[i], "--outvol") == 0) { 
      outvol = true;
    } else  if (strcmp(argv[i], "--outof") == 0) { 
      outof = true;
    } else  if (strcmp(argv[i], "--outgeom") == 0) { 
      outgeom = true;
    } else  if (strcmp(argv[i], "--help") == 0) { 
      printhelp();
      return 0;
    }
  }
  
// CHECK PARAMETERS
  if (filetype == 0) {
    cerr << "ERROR: Invalid filetype." << endl;
    error = true;
  }
  
  if (numprocs <= 0) {
    cerr << "ERROR: Invalid number of processors." << endl;
    error = true;
  }
  
  if (fineness < 0) {
    cerr << "ERROR: Invalid fineness level." << endl;
    error = true;
 
  }
  
  if (maxh < 0) { 
    cerr << "ERROR: Invalid maxh level." << endl;
    error = true; 
  }
  
  if (numlevels < 0) {
    cerr << "ERROR: Invalid number of refinement levels." << endl;
    error = true;
  }
  
  if (error) {
    printhelp();
    return 1;
  } 


   int *geoid;

   double t1, t2, t3, t4, tpar_, tpar;
   
   MPI_Init(&argc,&argv) ;
   MPI_Comm_rank(MPI_COMM_WORLD,&mypid);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   
// PRINT PARAMETERS
  if (mypid == 0) {  
    cout << "The parameters are as follows: " << endl
        << "Filename: " << geofilename << endl
        << "Filetype: " << filetype << endl
        << "Number of procs: " << numprocs << endl
        << "Fineness: " << fineness << ", maxh: " << maxh << endl
        << "Number of refinement levels: " << numlevels << endl
        << "Output: ";
    outvol ? cout << "Netgen mesh, " : cout << "0";
    outof ? cout << "Openfoam, " : cout << "";
    outgeom ? cout << "Geomview" : cout << "";
    cout << endl;
  }
      
   // initialize the Netgen Core library
   Ng_Init();

   // create the mesh structure
   mesh = Ng_NewMesh();
   
   nglib::Ng_Meshing_Parameters mp;
   mp.maxh = maxh ;
   mp.fineness = fineness ;

   MPI_Barrier(MPI_COMM_WORLD);
   t1 = MPI_Wtime();
   
  if (filetype == 1) {
    stl_geom = Ng_STL_LoadGeometry(geofilename);
    if (!stl_geom) { cerr << "ERROR: Error reading STL file " << geofilename << endl; return 1; }
    
    double t_ = MPI_Wtime();
    ng_res = Ng_STL_InitSTLGeometry(stl_geom);
    if (ng_res != NG_OK) { cerr << "ERROR: Error initializing STL geometry. " << endl; return 1; }    
    
    ng_res = Ng_STL_MakeEdges(stl_geom, mesh, &mp);
    if (ng_res != NG_OK) { cerr << "ERROR: Error making edges from STL geometry. " << endl; return 1; } 
    
    ng_res = Ng_STL_GenerateSurfaceMesh(stl_geom, mesh, &mp);
    if (ng_res != NG_OK) { cerr << "ERROR: Error generating surface mesh from STL geometry. " << endl; return 1; }
  if (mypid == 0) {printf("Generating surface mesh took %.2f\n",MPI_Wtime()-t_);}    
        
    surf_p = Ng_GetNP(mesh);
    surf_s = Ng_GetNSE(mesh);
    t_ = MPI_Wtime();
    ng_res = Ng_GenerateVolumeMesh (mesh, &mp);
    if (ng_res != NG_OK) { cerr << "ERROR: Error generating volume mesh from STL geometry. " << endl; return 1; }
  if (mypid == 0) {printf("Generating volume mesh took 3  %.2f\n",MPI_Wtime()-t_);}    
    
    vol_p = Ng_GetNP(mesh);
    vol_s = Ng_GetNSE(mesh);
    vol_v = Ng_GetNE(mesh);
  }
  if (filetype == 2) {
    Ng_CSG_Geometry *geom = new Ng_CSG_Geometry(); //loaded geometry
    double t_ = MPI_Wtime();
    ng_res = Ng_CSG_GenerateMesh (geofilename, mesh, geom, mp, nglib::MESH_SURFACE,
            &surf_p, &surf_s, &surf_v); //loaded geometry to geom, generated volume mesh
    if (ng_res != NG_OK) { cerr << "ERROR: Error generating surface mesh from CSG geometry. " << endl; return 1; }
  if (mypid == 0) {printf("Generating surface mesh took %.2f\n",MPI_Wtime()-t_);}    

/*
    t_ = MPI_Wtime();
    ng_res = Ng_CSG_GenerateMesh (geofilename, mesh, geom, mp, nglib::MESH_VOLUME,
            &vol_p, &vol_s, &vol_v); //loaded geometry to geom, generated volume mesh
    if (ng_res != NG_OK) { cerr << "ERROR: Error generating volume mesh from CSG geometry. " << endl; return 1; }
  if (mypid == 0) {printf("Generating volume mesh took 2  %.2f\n",MPI_Wtime()-t_);}    
*/
    t_ = MPI_Wtime();

    ng_res = Ng_CSG_GenerateVolumeMesh (mesh, mp, &vol_p, &vol_s, &vol_v);
  if (mypid == 0) {printf("Generating volume mesh took 3  %.2f\n",MPI_Wtime()-t_);}    
    if (ng_res != NG_OK) { cerr << "ERROR: Error generating volume mesh from CSG geometry. " << endl; return 1; }
  }

  
    sprintf(filename,"%s.serial",geofilename);
   outgeomview(mypid,numprocs,mesh,filename) ;
   

#ifdef DEBUG
  printf("DEBUG: [%5d]Global coarse surface mesh - No. of P/S: %d %d\n",surf_p,surf_s) ;
#endif
    
  maxbarycoord =  1  ;
  maxbarycoord =  1  << numlevels ;
  

  // Barrier 
  //timer serial finish 
  MPI_Barrier(MPI_COMM_WORLD);
  t2 = MPI_Wtime();

  nglib::Ng_Mesh *submesh ;
  submesh = Ng_NewMesh();
   
  partwmetis(mesh,numprocs) ; 
  cre_part_surf_mesh(mesh,numprocs,mypid,submesh) ;

  // TEST
  if (outof) {
    sprintf(filename, "mesh.init");
    My_Ng_Export_OF (mesh, geofilename, filename);
  }

  unirefine(submesh,numlevels) ; 
  computeadj(mypid) ;
  

//  sprintf(filename, "after-ref.%s.%d.vol",geofilename,mypid);
//  Ng_SaveMesh(submesh,filename);

#ifdef DEBUG
  printf("DEBUG: [%5d]Generating volume mesh\n",mypid);
#endif
   
  unsigned long int tpsc[3], psc2[3];
  int psc[3];
  psc[0]=psc[1]=psc[2]=0 ;
  tpsc[0]=tpsc[1]=tpsc[2]=0 ;
  Ng_DeleteMesh (mesh) ;  // NNN

  mp.fineness = fineness ;
  mp.maxh = maxh ;
  mp.second_order = 0; 
  
  if (mypid == 0) {
    printf("RESULT: SUBMESH ELEMENT #: %d\n", Ng_GetNSE(submesh));
  }

  if (filetype == 1) {
    double t_ = MPI_Wtime();
    ng_res = Ng_GenerateVolumeMesh(submesh, &mp);  
  if (mypid == 0) {printf("Generating volume mesh parallel took %.2f\n",MPI_Wtime()-t_);}    
    if (ng_res != NG_OK) { cerr << "ERROR: Error generating volume submesh from STL geometry. " << endl; return 1; }        
    psc[0] = Ng_GetNP(submesh);
    psc[1] = Ng_GetNSE(submesh);
    psc[2] = Ng_GetNE(submesh);
  }
  if (filetype == 2) {
    double t_ = MPI_Wtime();
    ng_res = Ng_CSG_GenerateVolumeMesh (submesh, mp, &psc[0], &psc[1], &psc[2]);
  if (mypid == 0) {printf("Generating volume mesh parallel took %.2f\n",MPI_Wtime()-t_);}    
    if (ng_res != NG_OK) { cerr << "ERROR: Error generating volume submesh from CSG geometry. " << endl; return 1; }        
  }

   //Barrier
   //timer parallel finish 
   MPI_Barrier(MPI_COMM_WORLD);
   t3 = MPI_Wtime();

  if (outvol) {
    sprintf(filename,"%s.%d.vol",geofilename,mypid) ;
    Ng_SaveMesh(submesh,filename);
  }
  
  if (outgeom) {
    outgeomview(mypid,numprocs,submesh,geofilename) ;
  }
  
  if (outof) {
    sprintf(filename, "processor%d",mypid);
    My_Ng_Export_OF (submesh, geofilename, filename);
  }

   psc2[0] = psc[0];
   psc2[1] = psc[1];
   psc2[2] = psc[2];
   MPI_Reduce(psc2,tpsc,3,MPI_UNSIGNED_LONG,MPI_SUM,0,MPI_COMM_WORLD) ; 
#ifdef DEBUG
   if (mypid == 0) {
      printf("MY: Total No. of P/S/V: %lu %lu %lu\n",tpsc[0],tpsc[1],tpsc[2]) ;
   } 
   outbaryc(submesh,mypid);
#endif
  com_barycoords(submesh,MPI_COMM_WORLD,numprocs,mypid) ;

   if (mypid == 0) {
      printf("RESULT: %s %d %d %.3f %.3f %d %d %d %d %lu %lu %lu %.2f %.2f\n",
                geofilename, numprocs, numprocs, mp.maxh, mp.fineness, numlevels,
                surf_p, surf_s, vol_v, // after surf & after volume mesh
                tpsc[0], tpsc[1], tpsc[2],
                t2-t1, t3-t2); // after refinement
   }

   MPI_Finalize() ;
   exit(0) ;   
   
   

}

