bin_PROGRAMS = statgen-mi.out

AM_CXXFLAGS = $(BOOST_CPPFLAGS) -ggdb -Wall -std=c++17
statgen_mi_out_SOURCES = statgen-mi/annotations.h statgen-mi/boost_random.cc statgen-mi/boost_random.h statgen-mi/cargs.cc statgen-mi/cargs.h statgen-mi/config_interpreter.cc statgen-mi/config_interpreter.h statgen-mi/custom_handlers.cc statgen-mi/custom_handlers.h statgen-mi/custom_handlers_emmaxkin.cc statgen-mi/custom_handlers_impute2.cc statgen-mi/custom_handlers_placeholder.cc statgen-mi/custom_handlers_plink.cc statgen-mi/custom_handlers_snptest.cc statgen-mi/custom_handlers_userdefined.cc statgen-mi/genotype_draw_methods.h statgen-mi/main.cc statgen-mi/mi_calculator.cc statgen-mi/mi_calculator.h statgen-mi/parameters.cc statgen-mi/parameters.h statgen-mi/prob_vector.cc statgen-mi/prob_vector.h statgen-mi/program_description.h statgen-mi/submission_formatter.cc statgen-mi/submission_formatter.h statgen-mi/transpose_matrix.cc statgen-mi/transpose_matrix.h statgen-mi/utilities.h
statgen-mi_out_LDADD = $(BOOST_LDFLAGS) -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lyaml-cpp
EXTRA_DIST = examples
statgen_midir = $(prefix)/@PACKAGE@/examples
statgen_mi_DATA = examples/example.sh examples/mi.imputeddataset.config examples/mi.program.config examples/mi.queue.config examples/test.gen examples/test.sample
dist_doc_DATA = README.md
ACLOCAL_AMFLAGS = -I m4
## TAP support
TEST_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(top_srcdir)/tap-driver.sh
TESTS = tests/placeholder.test
EXTRA_DIST = $(TESTS)
