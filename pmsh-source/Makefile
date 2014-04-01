CC        = mpic++ -g 
NETGENDIR = /home/seren/INSTALL
SOURCEDIR = /home/seren/FINAL/netgen-curie-jan2014/libsrc
INCDIRS   = -I$(NETGENDIR)/include  -I$(SOURCEDIR)/interface
LIBDIRS   = -L$(NETGENDIR)/lib -L./
LIBS      = -lnglib -linterface -lmetis -lm # -lmpi++ -lmpi
SWITCHES = -Wall -DSELF_CHECK
PREDCXXFLAGS = -O
#CC = icpc
OBJS      = 				compute_keys5.o    \
					eig_jacobi5.o       \
					inertia_xform5.o    \
					seq_rec_subdiv5.o   \
					vis5.o             \
					verify5.o          \
					unirefine7.o        \
					cresurmesh5.o 	\
					partwmetis.o \
					barycentric2.o \
					output.o \
					com_sr_datatype.o \
					com_num_recvs.o \
					com_length_from.o

MAINOBJS = 				main23netgen.o

PROGRAMS =				pmsh_netgen pmsh_netgen_export

all:        mynetgen

install:		$(PROGRAMS)
			install -m 0755 $(PROGRAMS) JOBSEND/

mynetgen:  		$(OBJS) $(MAINOBJS)
			$(CC) $(OBJS) main23netgen.o $(INCDIRS)  $(LIBDIRS) $(LIBS) -o pmsh_netgen
			$(CC) $(OBJS) main23netgenexport.o $(INCDIRS) $(LIBDIRS) $(LIBS) -o pmsh_netgen_export
#			$(CC) $(OBJS) main23solo.o $(INCDIRS) $(LIBDIRS) $(LIBS) -o solo_netgen

#main23solo.o:       	main23solo.cpp defns5.h
#			$(CC) -c $(INCDIRS) -DNETGEN main23solo.cpp -o main23solo.o

main23netgen.o:		main23.cpp defns5.h
			$(CC) -c $(INCDIRS) -DNETGEN main23.cpp -o main23netgen.o $(LIBS)
			$(CC) -c $(INCDIRS) -DNETGEN -DEXPORT main23.cpp -o main23netgenexport.o $(LIBS)

compute_keys5.o:	compute_keys5.cpp defns5.h
			$(CC) -c $(INCDIRS) compute_keys5.cpp
			
eig_jacobi5.o:		eig_jacobi5.cpp defns5.h
			$(CC) -c $(INCDIRS) eig_jacobi5.cpp 
		
inertia_xform5.o:	inertia_xform5.cpp defns5.h
			$(CC) -c $(INCDIRS) inertia_xform5.cpp 
			
seq_rec_subdiv5.o:	seq_rec_subdiv5.cpp defns5.h
			$(CC) -c $(INCDIRS) seq_rec_subdiv5.cpp 
			
cresurmesh5.o:		cresurmesh5.cpp defns5.h
			$(CC) -c $(INCDIRS) cresurmesh5.cpp
			
vis5.o:			vis5.cpp defns5.h
			$(CC) -c $(INCDIRS) vis5.cpp			
			 
verify5.o:		verify5.cpp defns5.h
			$(CC) -c $(INCDIRS) verify5.cpp
		
unirefine7.o:		unirefine7.cpp defns5.h
			$(CC) -c $(INCDIRS) unirefine7.cpp

partwmetis.o:           partwmetis.cpp defns5.h
			$(CC) -c $(INCDIRS) partwmetis.cpp

barycentric2.o:		barycentric2.cpp defns5.h
			$(CC) -c $(INCDIRS) barycentric2.cpp

output.o:		output.cpp defns5.h
			$(CC) -c $(INCDIRS) output.cpp

com_sr_datatype.o:	com_sr_datatype.cpp defns5.h
			$(CC) -c $(INCDIRS) com_sr_datatype.cpp

com_num_recvs.o:	com_num_recvs.cpp defns5.h
			$(CC) -c $(INCDIRS) com_num_recvs.cpp

com_length_from.o:	com_length_from.cpp defns5.h
			$(CC) -c $(INCDIRS) com_length_from.cpp

