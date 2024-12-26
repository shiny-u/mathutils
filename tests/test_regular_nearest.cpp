//
// Created by frongere on 18/05/22.
//

#include "MathUtils/RegularGridNearest.h"
#include "MathUtils/VectorGeneration.h"

#include <gtest/gtest.h>

using namespace mathutils;

TEST(RegularGridNearest, test_regular_grid_nearest) {

  using Array = boost::multi_array<int, 2>;
  Array var(boost::extents[3][3]);

  std::vector<double> xcoord = {0, 1, 2};
  std::vector<double> ycoord = {0, 1, 2};

  var[0][0] = 0;
  var[0][1] = 1;
  var[0][2] = 2;
  var[1][0] = 1;
  var[1][1] = 2;
  var[1][2] = 3;
  var[2][0] = 2;
  var[2][1] = 3;
  var[2][2] = 4;

  RegularGridNearest<int, 2, double> nearest;

  nearest.AddCoord(xcoord);
  nearest.AddCoord(ycoord);

  nearest.AddVar(var);

  ASSERT_EQ(nearest.Nearest({0., 0.}).val(), 0);
  ASSERT_EQ(nearest.Nearest({0.25, 0.25}).val(), 0);
  ASSERT_EQ(nearest.Nearest({0.5, 0.25}).val(), 1);
  ASSERT_EQ(nearest.Nearest({0.75, 0.25}).val(), 1);
  ASSERT_EQ(nearest.Nearest({1.25, 0.25}).val(), 1);
  ASSERT_EQ(nearest.Nearest({1.75, 0.25}).val(), 2);

  ASSERT_EQ(nearest.Nearest({0.25, 0.75}).val(), 1);
  ASSERT_EQ(nearest.Nearest({0.75, 0.75}).val(), 2);
  ASSERT_EQ(nearest.Nearest({1.25, 0.75}).val(), 2);
  ASSERT_EQ(nearest.Nearest({1.75, 0.75}).val(), 3);

  ASSERT_EQ(nearest.Nearest({0.25, 1.25}).val(), 1);
  ASSERT_EQ(nearest.Nearest({0.75, 1.25}).val(), 2);
  ASSERT_EQ(nearest.Nearest({1.25, 1.25}).val(), 2);
  ASSERT_EQ(nearest.Nearest({1.75, 1.25}).val(), 3);

  ASSERT_EQ(nearest.Nearest({0.25, 1.75}).val(), 2);
  ASSERT_EQ(nearest.Nearest({0.75, 1.75}).val(), 3);
  ASSERT_EQ(nearest.Nearest({1.25, 1.75}).val(), 3);
  ASSERT_EQ(nearest.Nearest({1.75, 1.75}).val(), 4);

}

TEST(RegularGridNearest, test_get_surrounding_nodes) {
  using Array = boost::multi_array<int, 2>;
  Array var(boost::extents[3][3]);

  std::vector<double> xcoord = {0, 1, 2};
  std::vector<double> ycoord = {0, 1, 2};

  var[0][0] = 0;
  var[0][1] = 1;
  var[0][2] = 2;
  var[1][0] = 1;
  var[1][1] = 2;
  var[1][2] = 3;
  var[2][0] = 2;
  var[2][1] = 3;
  var[2][2] = 4;

  RegularGridNearest<int, 2> nearest;

  nearest.AddCoord(xcoord);
  nearest.AddCoord(ycoord);

  nearest.AddVar(var);

  auto nodes = nearest.GetSurroundingGridNodes({1.5, 1.5});

  ASSERT_EQ(nodes[0].val(), 2);
  ASSERT_EQ(nodes[1].val(), 3);
  ASSERT_EQ(nodes[2].val(), 3);
  ASSERT_EQ(nodes[3].val(), 4);

}