# ================================================================
# Makefile for project rvrcm
# Automatically generated from "rvrcm.mki" at Tue Sep  8 09:22:09 2009

# yamm v1.0
# John Kerl
# 2002/05/04
# ================================================================


INCLUDE_DIRS =
LIB_DIRS = -L.
DEFINES =
MISC_CFLAGS =
MISC_LFLAGS = -lrcm -lm
EXTRA_DEPS = ./librcm.a
COMPILE_FLAGS = -c $(INCLUDE_DIRS) $(DEFINES) $(MISC_CFLAGS)
LINK_FLAGS =  $(LIB_DIRS) $(MISC_LFLAGS)

build: mk_obj_dir ./rvrcm

mk_obj_dir:
	mkdir -p ./rcm_objs

./rcm_objs/rvrcm.o:  energy.h interactions.h mcmc_params.h mcmc_rvs.h metro_stats.h mtrand.h pmt.h points.h psdes.h rcmrand.h rho.h rvrcm.c seps.h stats.h times.h urandom.h util.h vector_misc.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  rvrcm.c -o ./rcm_objs/rvrcm.o

OBJS = \
	./rcm_objs/rvrcm.o

./rvrcm: $(OBJS) $(EXTRA_DEPS)
	gcc $(OPTLFLAGS) $(OBJS) -o ./rvrcm $(LINK_FLAGS)

clean:
	-@rm -f $(OBJS)
	-@rm -f ./rvrcm
