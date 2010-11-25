# set root of cxxtest tree and add it to include directories
SET(CXXTEST_INCLUDE_DIR "${CONVERTERS_SOURCE_DIR}/../cxxtest")
INCLUDE_DIRECTORIES(${CXXTEST_INCLUDE_DIR})

# test generator
SET(CXXTEST_TESTGEN_EXECUTABLE "${CXXTEST_INCLUDE_DIR}/python/scripts/cxxtestgen")

#=============================================================
# CXXTEST_ADD_TEST (public macro)
#=============================================================
MACRO(CXXTEST_ADD_TEST _cxxtest_testname _cxxtest_outfname)
    SET(_cxxtest_real_outfname ${CMAKE_CURRENT_BINARY_DIR}/${_cxxtest_outfname})

    # generate test file
    ADD_CUSTOM_COMMAND(
        OUTPUT  ${_cxxtest_real_outfname}
        DEPENDS ${ARGN}
        COMMAND ${CXXTEST_TESTGEN_EXECUTABLE}
        --error-printer -o ${_cxxtest_real_outfname} ${ARGN}
    )

    SET_SOURCE_FILES_PROPERTIES(${_cxxtest_real_outfname} PROPERTIES GENERATED true)
    
    #compile test file
    ADD_EXECUTABLE(${_cxxtest_testname} ${_cxxtest_real_outfname})
    
    # run compiled test
    ADD_TEST(${_cxxtest_testname} ${CMAKE_CURRENT_BINARY_DIR}/${_cxxtest_testname})
    
ENDMACRO(CXXTEST_ADD_TEST)
