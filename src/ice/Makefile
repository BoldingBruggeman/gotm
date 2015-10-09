#
# Makefile to build the the ice interaction library - libice.a
#

include ../Rules.make

LIB	= ../../libice$(buildtype).a
LIB	:= $(LIBDIR)/libice$(buildtype).a

LIB1	= ../../libice_utils$(buildtype).a
LIB1	:= $(LIBDIR)/libice_utils$(buildtype).a

DOCSRC	=\
ice.F90			\

MOD	= \
${LIB1}(ice_uvic.o)	\
${LIB1}(ice_winton.o)	\
${LIB}(ice.o)

OBJ	= 

all: $(MOD) $(OBJ)
	$(MOVE_MODULES_COMMAND)

$(OBJ): $(MOD)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/ice.tex
	touch doc

clean:
	$(RM) -f ${LIB} ${LIB1} $(MODDIR)/ice.{m,mod} doc

realclean: clean
	$(RM) -f *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 1999 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
