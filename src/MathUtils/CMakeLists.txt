
# HEADER_FILES
set(MATHUTILS_HEADERS
		${CMAKE_CURRENT_BINARY_DIR}/version.h
        Constants.h
        Angles.h
        Unit.h
        Interp.h
        Interp1d.h
        LookupTable1D.h
        Interp2d.h
        LookupTable2D.h
        Maths.h
        VectorGeneration.h
        Check.h
        Vector2d.h
        Vector3d.h
        Transform3d.h
        Matrix.h
        Integration1d.h
        FFT.h
        Complex.h
        StdVector.h
        Functions.h
        Convolution.h
        VectorN.h
	    Geometry.h
        Integration2d.h
        Integration2dTriangle.h
        QuadratureTables.h
        Integrand.h
        Polynomial.h
        BoostFunctions.h
		ChebyshevSeries1d.h
        ChebyshevSeries2d.h
        ChebyshevSeries3d.h
		PowerSeries1d.h
        PowerSeries2d.h
        PowerSeries3d.h
		RegularGridInterpolator.h
        )

# Building target
add_library(MathUtils INTERFACE)
add_dependencies(MathUtils check_git_${PROJECT_NAME}) # For git_watcher to fetch git informations before effective build
target_link_libraries(MathUtils INTERFACE Eigen3::Eigen)

add_library(MathUtils::MathUtils ALIAS MathUtils)

list(APPEND MATHUTILS_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/..")

target_include_directories(MathUtils INTERFACE
        $<BUILD_INTERFACE:${MATHUTILS_INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>)

target_include_directories(MathUtils INTERFACE
		$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/..>)

#include(GNUInstallDirs)
#install(TARGETS MathUtils EXPORT MathUtils-targets
#        DESTINATION ${CMAKE_INSTALL_LIBDIR})
#
#install(FILES ${MATHUTILS_HEADERS}
#        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MathUtils)
