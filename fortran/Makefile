F77  = ifort
FFLAGS  = -O -fpp -qopenmp
LIB_FLAGS = -L$(FFTW_DIR) -lfftw3f -lfftw3 -lfftw3l -lm
#INCLUDES = -I$(FFTW3_INCDIR)
#OPENMP_LIB_FLAGS = -Wl,-static -L$(FFTW3_LIBDIR) -lfftw3f_omp -lfftw3_omp -lfftw3l_omp -lfftw3f -lfftw3 -lfftw3l -lm -Wl,-dy
#LIB_FLAGS = -L/home2/czhao/lib/fftw-3.3.3/lib -lfftw3f -lfftw3 -lfftw3l -lm -Wl,-dy
INCLUDES = -I$(FFTW_INC)
#OPENMP_LIB_FLAGS = -L/home2/czhao/lib/fftw-3.3.3/lib -lfftw3f_omp -lfftw3_omp -lfftw3l_omp -lfftw3f -lfftw3 -lfftw3l -lm -Wl,-dy
OPENMP_LIB_FLAGS = -L${FFTW_DIR} -lfftw3f_omp -lfftw3_omp -lfftw3l_omp -lfftw3f -lfftw3 -lfftw3l -lm -Wl,-dy

EXECUTABLES = \
 pro_EZmock_pk_cf
# pro_EZmock_v0  pro_EZmock_simPDF pro_EZmock_pk_cf pro_fittingEZmock
#    pro_compute_pk_box pro_compute_pk_box_CICbinning  pro_compute_pk_divide_box_CICbinning \
#    pro_compute_pk_divide_box_CICbinning_zdist
#   pro_sdss_pk_fftw
# pro_CF_code_box   pro_mono_covar pro_compute_chi2  pro_ratio_wiggle_vs_nowiggle \

default:	${EXECUTABLES}

all:	${EXECUTABLES}


clean:	
	rm -f ${EXECUTABLES} *.[oadD] *.[mM][oO][dD]

#pro_EZmock_v0: pro_EZmock_v0.o  correlation_func_subs.o  numerical_recipe_subs.o
#	${F77} ${FFLAGS} ${LDFLAGS} -o $@ $^ ${OPENMP_LIB_FLAGS} ${INCLUDES}
#pro_EZmock_simPDF: pro_EZmock_simPDF.o  correlation_func_subs.o  numerical_recipe_subs.o
#	${F77} ${FFLAGS} ${LDFLAGS} -o $@ $^ ${OPENMP_LIB_FLAGS} ${INCLUDES}
pro_EZmock_pk_cf: pro_EZmock_pk_cf.o  correlation_func_subs.o  numerical_recipe_subs.o
	${F77} ${FFLAGS} ${LDFLAGS} ${INCLUDES} -o $@ $^ ${OPENMP_LIB_FLAGS} 
#pro_fittingEZmock: pro_fittingEZmock.o
#	${F77} ${FFLAGS} ${LDFLAGS} -o $@ $^ ${OPENMP_LIB_FLAGS} ${INCLUDES}

#pro_ratio_wiggle_vs_nowiggle: pro_ratio_wiggle_vs_nowiggle.o
#	${F77} ${FFLAGS} ${LDFLAGS} -o $@ $^ ${LIB_FLAGS} ${INCLUDES}

#pro_compute_chi2: pro_compute_chi2.o
#	${F77} ${FFLAGS} ${LDFLAGS} -o $@ $^ ${LIB_FLAGS} ${INCLUDES}

#pro_mono_covar: pro_mono_covar.o
#	${F77} ${FFLAGS} ${LDFLAGS} -o $@ $^ ${LIB_FLAGS} ${INCLUDES}

#pro_CF_code_box: pro_CF_code_box.o
#	${F77} ${FFLAGS} ${LDFLAGS} -o $@ $^ ${LIB_FLAGS} ${INCLUDES}


.f.o:
	${F77} ${FFLAGS} ${INCLUDES} -c $< ${LIB_FLAGS} 

