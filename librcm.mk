# ================================================================
# Makefile for project librcm
# Automatically generated from "librcm.mki" at Tue Sep  8 09:22:09 2009

# yamm v1.0
# John Kerl
# 2002/05/04
# ================================================================


INCLUDE_DIRS =
LIB_DIRS =
DEFINES =
MISC_CFLAGS =
MISC_LFLAGS =
EXTRA_DEPS =
COMPILE_FLAGS = -c $(INCLUDE_DIRS) $(DEFINES) $(MISC_CFLAGS)
LINK_FLAGS =  $(LIB_DIRS) $(MISC_LFLAGS)

build: mk_obj_dir ./librcm.a

mk_obj_dir:
	mkdir -p ./librcm_objs

./librcm_objs/mcmc_params.o:  checks.h interactions.h mcmc_params.c mcmc_params.h metro_stats.h mtrand.h points.h psdes.h rcmrand.h urandom.h util.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  mcmc_params.c -o ./librcm_objs/mcmc_params.o

./librcm_objs/mcmc_rvs.o:  interactions.h mcmc_params.h mcmc_rvs.c mcmc_rvs.h metro_stats.h pmt.h points.h rho.h stats.h util.h vector_misc.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  mcmc_rvs.c -o ./librcm_objs/mcmc_rvs.o

./librcm_objs/util.o:  times.h util.c util.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  util.c -o ./librcm_objs/util.o

./librcm_objs/vector_misc.o:  vector_misc.c vector_misc.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  vector_misc.c -o ./librcm_objs/vector_misc.o

./librcm_objs/stats.o:  stats.c stats.h util.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  stats.c -o ./librcm_objs/stats.o

./librcm_objs/times.o:  times.c times.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  times.c -o ./librcm_objs/times.o

./librcm_objs/gasdev.o:  gasdev.c mtrand.h psdes.h rcmrand.h urandom.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  gasdev.c -o ./librcm_objs/gasdev.o

./librcm_objs/psdes.o:  psdes.c psdes.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  psdes.c -o ./librcm_objs/psdes.o

./librcm_objs/urandom.o:  urandom.c urandom.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  urandom.c -o ./librcm_objs/urandom.o

./librcm_objs/mtrand.o:  mtrand.c
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  mtrand.c -o ./librcm_objs/mtrand.o

./librcm_objs/rcmrand.o:  mtrand.h psdes.h rcmrand.c rcmrand.h urandom.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  rcmrand.c -o ./librcm_objs/rcmrand.o

./librcm_objs/points.o:  checks.h mtrand.h pmt.h points.c points.h psdes.h rcmrand.h urandom.h util.h wormidx.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  points.c -o ./librcm_objs/points.o

./librcm_objs/dotplot.o:  dotplot.c dotplot.h pmt.h points.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  dotplot.c -o ./librcm_objs/dotplot.o

./librcm_objs/energy.o:  checks.h energy.c energy.h interactions.h pmt.h points.h util.h worm_moves.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  energy.c -o ./librcm_objs/energy.o

./librcm_objs/interactions.o:  interactions.c interactions.h util.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  interactions.c -o ./librcm_objs/interactions.o

./librcm_objs/pmt.o:  mtrand.h pmt.c pmt.h points.h psdes.h rcmrand.h urandom.h util.h wormidx.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  pmt.c -o ./librcm_objs/pmt.o

./librcm_objs/rho.o:  pmt.h points.h rho.c rho.h util.h vector_misc.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  rho.c -o ./librcm_objs/rho.o

./librcm_objs/thermalization.o:  thermalization.c thermalization.h util.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  thermalization.c -o ./librcm_objs/thermalization.o

./librcm_objs/metropolis.o:  energy.h interactions.h mcmc_params.h metro_stats.h metropolis.c metropolis.h mtrand.h pmt.h points.h psdes.h rcmrand.h urandom.h util.h worm_moves.h wormidx.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  metropolis.c -o ./librcm_objs/metropolis.o

./librcm_objs/metro_stats.o:  metro_stats.c metro_stats.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  metro_stats.c -o ./librcm_objs/metro_stats.o

./librcm_objs/worm_moves.o:  worm_moves.c worm_moves.h
	gcc $(OPTCFLAGS) -Wall -Werror $(COMPILE_FLAGS)  worm_moves.c -o ./librcm_objs/worm_moves.o

OBJS = \
	./librcm_objs/mcmc_params.o \
	./librcm_objs/mcmc_rvs.o \
	./librcm_objs/util.o \
	./librcm_objs/vector_misc.o \
	./librcm_objs/stats.o \
	./librcm_objs/times.o \
	./librcm_objs/gasdev.o \
	./librcm_objs/psdes.o \
	./librcm_objs/urandom.o \
	./librcm_objs/mtrand.o \
	./librcm_objs/rcmrand.o \
	./librcm_objs/points.o \
	./librcm_objs/dotplot.o \
	./librcm_objs/energy.o \
	./librcm_objs/interactions.o \
	./librcm_objs/pmt.o \
	./librcm_objs/rho.o \
	./librcm_objs/thermalization.o \
	./librcm_objs/metropolis.o \
	./librcm_objs/metro_stats.o \
	./librcm_objs/worm_moves.o

./librcm.a: $(OBJS) $(EXTRA_DEPS)
	ar r librcm.a ./librcm.a $(OBJS)

clean:
	-@rm -f $(OBJS)
	-@rm -f ./librcm.a
