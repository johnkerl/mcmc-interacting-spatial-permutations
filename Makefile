# This file was automatically generated from Makefile.in by genmk.

opt:
	export OPTCFLAGS="-O3" OPTLFLAGS=""; make -ef librcm.mk
	export OPTCFLAGS="-O3" OPTLFLAGS=""; make -ef mcrcm.mk
	export OPTCFLAGS="-O3" OPTLFLAGS=""; make -ef rvrcm.mk
	export OPTCFLAGS="-O3" OPTLFLAGS=""; make -ef ccump.mk

profile:
	export OPTCFLAGS="-g -O3 -pg" OPTLFLAGS="-g -pg"; make -ef librcm.mk
	export OPTCFLAGS="-g -O3 -pg" OPTLFLAGS="-g -pg"; make -ef mcrcm.mk
	export OPTCFLAGS="-g -O3 -pg" OPTLFLAGS="-g -pg"; make -ef rvrcm.mk
	export OPTCFLAGS="-g -O3 -pg" OPTLFLAGS="-g -pg"; make -ef ccump.mk

gcov:
	export OPTCFLAGS="-g -fprofile-arcs -ftest-coverage" OPTLFLAGS="-g -fprofile-arcs -ftest-coverage"; make -ef librcm.mk
	export OPTCFLAGS="-g -fprofile-arcs -ftest-coverage" OPTLFLAGS="-g -fprofile-arcs -ftest-coverage"; make -ef mcrcm.mk
	export OPTCFLAGS="-g -fprofile-arcs -ftest-coverage" OPTLFLAGS="-g -fprofile-arcs -ftest-coverage"; make -ef rvrcm.mk
	export OPTCFLAGS="-g -fprofile-arcs -ftest-coverage" OPTLFLAGS="-g -fprofile-arcs -ftest-coverage"; make -ef ccump.mk

debug:
	export OPTCFLAGS="-g" OPTLFLAGS="-g"; make -ef librcm.mk
	export OPTCFLAGS="-g" OPTLFLAGS="-g"; make -ef mcrcm.mk
	export OPTCFLAGS="-g" OPTLFLAGS="-g"; make -ef rvrcm.mk
	export OPTCFLAGS="-g" OPTLFLAGS="-g"; make -ef ccump.mk

build:
	make -f librcm.mk
	make -f mcrcm.mk
	make -f rvrcm.mk
	make -f ccump.mk

mk:
	yamm librcm.mki
	yamm mcrcm.mki
	yamm rvrcm.mki
	yamm ccump.mki

install:
	make -f librcm.mk                      install
	make -f mcrcm.mk                       install
	make -f rvrcm.mk                       install
	make -f ccump.mk                       install

clean:
	make -f librcm.mk                      clean
	make -f mcrcm.mk                       clean
	make -f rvrcm.mk                       clean
	make -f ccump.mk                       clean

tags: .PHONY
	ctags *.[ch]

.PHONY:

over: clean mk build
