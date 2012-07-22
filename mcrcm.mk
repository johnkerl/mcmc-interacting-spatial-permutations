# ================================================================
# Makefile for project mcrcm
# Automatically generated from "mcrcm.mki" at Tue Sep  8 09:22:09 2009

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

build: mk_obj_dir ./mcrcm

mk_obj_dir:
	mkdir -p ./rcm_objs

./rcm_objs/mcrcm.o:  dotplot.h energy.h interactions.h mcmc_params.h mcmc_rvs.h mcrcm.c metro_stats.h metropolis.h mtrand.h pmt.h points.h psdes.h rcmrand.h seps.h stats.h thermalization.h times.h urandom.h util.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  mcrcm.c -o ./rcm_objs/mcrcm.o

OBJS = \
	./rcm_objs/mcrcm.o

./mcrcm: $(OBJS) $(EXTRA_DEPS)
	gcc $(OPTLFLAGS) $(OBJS) -o ./mcrcm $(LINK_FLAGS)

clean:
	-@rm -f $(OBJS)
	-@rm -f ./mcrcm
