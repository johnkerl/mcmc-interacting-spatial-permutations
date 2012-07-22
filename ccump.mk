# ================================================================
# Makefile for project ccump
# Automatically generated from "ccump.mki" at Tue Sep  8 09:22:09 2009

# yamm v1.0
# John Kerl
# 2002/05/04
# ================================================================


INCLUDE_DIRS =
LIB_DIRS = -L.
DEFINES =
MISC_CFLAGS =
MISC_LFLAGS = -lm
EXTRA_DEPS =
COMPILE_FLAGS = -c $(INCLUDE_DIRS) $(DEFINES) $(MISC_CFLAGS)
LINK_FLAGS =  $(LIB_DIRS) $(MISC_LFLAGS)

build: mk_obj_dir ./ccump

mk_obj_dir:
	mkdir -p ./ccump_objs

./ccump_objs/ccump.o:  ccump.c cump_params.h mtrand.h psdes.h rcmrand.h urandom.h util.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  ccump.c -o ./ccump_objs/ccump.o

./ccump_objs/cump_params.o:  cump_params.c cump_params.h mtrand.h psdes.h rcmrand.h urandom.h util.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  cump_params.c -o ./ccump_objs/cump_params.o

./ccump_objs/mtrand.o:  mtrand.c
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  mtrand.c -o ./ccump_objs/mtrand.o

./ccump_objs/rcmrand.o:  mtrand.h psdes.h rcmrand.c rcmrand.h urandom.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  rcmrand.c -o ./ccump_objs/rcmrand.o

./ccump_objs/urandom.o:  urandom.c urandom.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  urandom.c -o ./ccump_objs/urandom.o

./ccump_objs/util.o:  times.h util.c util.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  util.c -o ./ccump_objs/util.o

./ccump_objs/times.o:  times.c times.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  times.c -o ./ccump_objs/times.o

OBJS = \
	./ccump_objs/ccump.o \
	./ccump_objs/cump_params.o \
	./ccump_objs/mtrand.o \
	./ccump_objs/rcmrand.o \
	./ccump_objs/urandom.o \
	./ccump_objs/util.o \
	./ccump_objs/times.o

./ccump: $(OBJS) $(EXTRA_DEPS)
	gcc $(OPTLFLAGS) $(OBJS) -o ./ccump $(LINK_FLAGS)

clean:
	-@rm -f $(OBJS)
	-@rm -f ./ccump
