
# Here we declare the different PATH, TAG and PATCH to get the FRyDoM dependencies

# Eigen
set(eigen_URL https://gitlab.com/libeigen/eigen.git  CACHE STRING "eigen repository URL")
set(eigen_TAG 3.3.7 CACHE STRING "eigen version")

# Units
set(units_URL https://github.com/nholthaus/units.git)
set(units_TAG v2.3.1 CACHE STRING "units version")

# Boost
set(boost_URL https://boostorg.jfrog.io/artifactory/main/release/1.75.0/source/boost_1_75_0.tar.gz)
set(boost_TAG 1.75 CACHE STRING "Boost version")
set(boost_FIND_TAG 1.66 CACHE STRING "Minimal version of Boost to find on the system")

# GoogleTest
set(googletest_URL https://github.com/google/googletest.git)
set(googletest_TAG release-1.10.0 CACHE STRING "googletest version")

# Boost
#set(boost_TAG 1.66)