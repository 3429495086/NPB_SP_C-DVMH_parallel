SHELL=/bin/sh
BENCHMARK=lu
BENCHMARKU=LU
VEC=
ACTION=exec

include ../config/make.def

include ../sys/make.common

#------------------- Lists of sources and headers -----------------------------#

HEAD = 	npbparams.h ${COMMON}/print_results.h ${COMMON}/timers.h \
		${COMMON}/wtime.h ${COMMON}/type.h


# List of files which can be merged into a single AST.
MERGE_SRC = lu.c read_input.c \
			domain.c setcoeff.c setbv.c exact.c setiv.c \
			erhs.c ssor$(VEC).c rhs$(VEC).c l2norm.c \
			jacld.c blts$(VEC).c jacu.c buts$(VEC).c error.c \
			pintgr.c verify.c ${COMMON}/print_results.c \
			${COMMON}/wtime.c

# List of files which can not be imported by TSAR,
# for example due to symbol name redifinition.
SINGLE_SRC = ${COMMON}/c_timers.c

#------------------------------------------------------------------------------#

include ../sys/make.vars


# npbparams.h is included by applu.incl
# The following rule should do the trick but many make programs (not gmake)
# will do the wrong thing and rebuild the world every time (because the
# mod time on header.h is not changed. One solution would be to 
# touch header.h but this might cause confusion if someone has
# accidentally deleted it. Instead, make the dependency on npbparams.h
# explicit in all the lines below (even though dependence is indirect). 

# applu.incl: npbparams.h

${PROGRAM}: config
	@if [ x$(VERSION) = xvec ] ; then	\
		${MAKE} VEC=_vec ${ACTION};		\
	elif [ x$(VERSION) = xVEC ] ; then	\
		${MAKE} VEC=_vec ${ACTION};		\
	else					\
		${MAKE} ${ACTION};			\
	fi

include ../sys/make.tsar

lu.ll lu.ast lu.o:					lu.c applu.incl npbparams.h
blts$(VEC).ll blts$(VEC).ast blts$(VEC).o:	blts$(VEC).c
buts$(VEC).ll buts$(VEC).ast buts$(VEC).o:	buts$(VEC).c
erhs.ll erhs.ast erhs.o:				erhs.c applu.incl npbparams.h
error.ll error.ast error.o:				error.c applu.incl npbparams.h
exact.ll exact.ast exact.o:				exact.c applu.incl npbparams.h
jacld.ll jacld.ast jacld.o:				jacld.c applu.incl npbparams.h
jacu.ll jacu.ast jacu.o:				jacu.c applu.incl npbparams.h
l2norm.ll l2norm.ast l2norm.o:			l2norm.c
pintgr.ll pintgr.ast pintgr.o:			pintgr.c applu.incl npbparams.h
read_input.ll read_input.ast read_input.o:	read_input.c applu.incl npbparams.h
rhs$(VEC).ll rhs$(VEC).ast rhs$(VEC).o:		rhs$(VEC).c applu.incl npbparams.h
setbv.ll setbv.ast setbv.o:				setbv.c applu.incl npbparams.h
setiv.ll setiv.ast setiv.o:				setiv.c applu.incl npbparams.h
setcoef.ll setcoeff.ast setcoeff.o:		setcoeff.c applu.incl npbparams.h
ssor$(VEC).ll ssor$(VEC).ast ssor$(VEC).o:	ssor$(VEC).c applu.incl npbparams.h
domain.ll domain.ast domain.o:			domain.c applu.incl npbparams.h
verify.ll verify.ast verify.o:			verify.c applu.incl npbparams.h

clean:
	- rm -f npbparams.h *.ast *.ll
	- rm -f *.o *~
