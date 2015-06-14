# Parallel Mesh Generation for OpenFOAM with Netgen #

![https://pmsh.googlecode.com/svn/prace.png](https://pmsh.googlecode.com/svn/prace.png)

This work was supported by the PRACE-3IP project funded in part
by the EUs 7th Framework Programme (FP7/2012-2014) under grant agreement
no. RI-312763.

---


OpenFOAM is an open source computational  fuid dynamics (CFD) package with a large user
base from many areas of engineering and science. An enablement tool called PMSH  is
developed to generate multi-billion element unstructured tetrahedral meshes for OpenFOAM.
PMSH is developed as a wrapper code around the popular open source sequential Netgen mesh
generator. Parallelization of the mesh generation process is carried out in five main
stages: (i) generation of a coarse volume mesh (ii) partitioning of the coarse mesh to get
sub-meshes, each of which is processed by a processor (iii) extraction and refinement of
coarse surface mesh to produce fine surface sub-meshes (iv) re-meshing of each fine surface
sub-mesh to get a final fine mesh (v) matching of partition boundary vertices followed by
global vertex numbering.  Test results obtained  on an SGI Altix ICE X system with 8192 cores and 14 TB of total memory confirm that our  approach does indeed enable us to generate multi-billion element meshes in a  scalable way.


---


## pmsh-netgen installation tutorial ##

# Install tcl to $INSTALLDIR

# Install tk to $INSTALLDIR

# Install togl to $INSTALLDIR

# Configure netgen-5.1 with :
```
        ./configure --with-tk=$INSTALLDIR/lib --with-tcl=$INSTALLDIR/lib --with-tclinclude=$INSTALLDIR/include --with-tkinclude=$INSTALLDIR/include/ --with-togl=$INSTALLDIR/lib/ --prefix=$INSTALLDIR
```
# Edit nglib/Makefile :
> Change the following line
```
bin_PROGRAMS = ng_vol$(EXEEXT) ng_stl$(EXEEXT)
```
> to
```
bin_PROGRAMS =
```
# Edit ng/Makefile. Netgen-5.1 has a bug. You have to change the following line:
```
        -L$(TK_BIN_DIR)/Togl1.7 $(TOGLLIBDIR) -lTogl $(LIBGLU) $(TK_LIB_SPEC) $(TCL_LIB_SPEC) $(MPI_LIBS) $(FFMPEG_LIBS) $(JPEGLIB_LIBS) $(PKG_LIBS) $(MKL_LIBS)
```
to
```
        -L$(TK_BIN_DIR)/Togl1.7 $(TK_BIN_DIR)/Togl1.7/libTogl1.7.so $(LIBGLU)  $(TK_LIB_SPEC) $(TCL_LIB_SPEC) $(MPI_LIBS) $(FFMPEG_LIBS) $(JPEGLIB_LIBS) $(PKG_LIBS)
```
# Copy everyting from NETGEN\_PATCH to netgen-5.1 source directory

# Under netgen source directory
```
        make
        make install
```

# You can now install pmsh-netgen. First edit the following lines in pmsh-source/Makefile
```
        NETGENDIR = #insert your installdir here
        SOURCEDIR = #insert netgen sourcedir here/libsrc
```

# make

# In order to run pmsh-netgen, you may use the following commands:
```
        ./pmsh_netgen --help
        mpirun -np 4 ./pmsh_netgen --stl onera-m6.stl --numcores 4 --fineness 0.025 --maxh 0.1 --outgeom
        mpirun -np 4 ./pmsh_netgen --geo cube.geo --numcores 3 --fineness 0.025 --maxh 0.01 --outof
```