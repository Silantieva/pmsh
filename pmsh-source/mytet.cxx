///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetcall.cxx                                                               //
//                                                                           //
// An example of how to call TetGen from another program by using the data   //
// type "tetgenio" and function "tetrahedralize()" of TetGen libaray.        //
//                                                                           //
// In order to run this example, you need the library of TetGen, you can get //
// the source code as well as the user's manul of TetGen from:               //
//                                                                           //
//            http://tetgen.berlios.de/index.html                            //
//                                                                           //
// Section 2 of the user's manual contains the information of how to compile //
// TetGen into a libaray.                                                    //
//                                                                           //
// The geometry used in this example (illustrated in Section 3.3 .1, Figure  //
// 12 of the user's manual) is a rectangluar bar consists of 8 points and 6  //
// facets (which are all rectangles). In additional, there are two boundary  //
// markers defined on its facets.                                            //
//                                                                           //
// This code illustrates the following basic steps:                          //
//   - at first create an input object "in", and set data of the geometry    //
//     into it.                                                              //
//   - then call function "tetrahedralize()" to create a quality mesh of the //
//     geometry with output in another object "out".                         //
// In addition, It outputs the geometry in the object "in" into two files    //
// (barin.node and barin.poly), and outputs the mesh in the object "out"     //
// into three files (barout.node, barout.ele, and barout.face).  These files //
// can be visualized by TetView.                                             //
//                                                                           //
// To compile this code into an executable program, do the following steps:  //
//   - compile TetGen into a library named "libtet.a" (see Section 2.1 of    //
//     the user's manula for compiling);                                     //
//   - Save this file into the same directory in which you have the files    //
//     "tetgen.h" and "libtet.a";                                            //
//   - compile it using the following command:                               //
//                                                                           //
//     g++ -o test tetcall.cxx -L./ -ltet                                    //
//                                                                           //
//     which will result an executable program named "test".                 //
//                                                                           //
// Please send your quesions, comments to Hang Si <si@wias-berlin.de>        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// main()  Create and refine a mesh using TetGen library.                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "tetgen.h" // Defined tetgenio, tetrahedralize().

int main(int argc, char *argv[])
{
  tetgenio in, out;
  tetgenbehavior b;
  b.geomview = 1;
  b.quality++;
  b.plc = 1;

  printf(" naber ");
  in.load_off(argv[1]);
  strcpy(b.outfilename, argv[2]);
  tetrahedralize("pqYNEFO", &in, NULL, "out");
  //out.save_nodes("naber.nodes");
  //out.save_faces("naber.faces");
  //tetrahedralize(&b, &in, NULL);
  

  
  return 0;
}
