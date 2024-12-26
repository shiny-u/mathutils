include(FetchContent)

FetchContent_Declare(units
  GIT_REPOSITORY ${units_URL}
  GIT_TAG ${units_TAG}
)
FetchContent_GetProperties(units)
if(NOT units_POPULATED)
  message(STATUS "******* FETCHING units dependency from ${PROJECT_NAME} (requested version: ${units_TAG}) *******")
  set(BUILD_TESTS OFF CACHE BOOL "")
  set(BUILD_DOCS OFF CACHE BOOL "")
  FetchContent_Populate(units)
  add_subdirectory(${units_SOURCE_DIR} ${units_BINARY_DIR})
endif()

