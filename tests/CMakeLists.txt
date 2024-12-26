
add_custom_target(mathutils_tests)

set(TESTS
        test_get_version
        test_unit
        test_angles
        test_interp1d
        test_lookup_table
        test_interp2d
        test_vector2d
        test_vector3d
        test_transform3d
        test_matrixMN
        test_integration1d
        test_fft
        test_spline
        test_vectorN
	    test_geometry
	    test_linear_system_solver
        test_integration2d
        test_polynomial
        test_regular_grid_interpolator
        test_regular_nearest
        )

IF(MSVC)
    ADD_DEFINITIONS("/MP")
    MESSAGE(STATUS "MSVC /MP option activated")
ENDIF(MSVC)

foreach(test ${TESTS})

    message(STATUS "...add ${test}")

    add_executable(${test} ${test}.cpp)
    if(${MATHUTILS_SKIP_LONG_TESTS})
        target_compile_definitions(${test} PUBLIC
                SKIP_LONG_TESTS)
    endif()
    target_link_libraries(${test} MathUtils::MathUtils gtest_main)

    add_dependencies(mathutils_tests ${test})

    gtest_discover_tests(${test})

endforeach()

# tests with Boost dependency
if (MATHUTILS_BUILD_BOOST_TESTS)
    include(Add_Boost)
    set(TESTS
            test_BoostFunctions
            test_ChebyshevSeries1d
            test_ChebyshevSeries2d
            test_ChebyshevSeries3d
            test_convolution
            test_PowerSeries1d
            test_PowerSeries2d
            test_PowerSeries3d
            )

    foreach(test ${TESTS})

        message(STATUS "...add ${test}")

        add_executable(${test} ${test}.cpp)
        if(${MATHUTILS_SKIP_LONG_TESTS})
            target_compile_definitions(${test} PUBLIC
                    SKIP_LONG_TESTS)
        endif()
        target_link_libraries(${test} MathUtils::MathUtils Boost::boost gtest_main)

        add_dependencies(mathutils_tests ${test})

        gtest_discover_tests(${test})

    endforeach()
endif()
