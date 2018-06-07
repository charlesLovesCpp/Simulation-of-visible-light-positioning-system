/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2016  <copyright holder> <email>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#ifndef COMMON_INCLUDE_H
#define COMMON_INCLUDE_H

// define the commonly included file to avoid a long include list
// for Eigen
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
using namespace Eigen;

// for Sophus
#include <sophus/se3.h>
#include <sophus/so3.h>
using Sophus::SO3;
using Sophus::SE3;

// for cv
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp> 
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/calib3d/calib3d.hpp>
using cv::Mat;

// std 
#include <vector>
#include <list>
#include <queue>
#include <memory>
#include <string>
#include <iostream>
#include <set>
#include <unordered_map>
#include <map>
#include <cmath>
#include <fstream>
#include <chrono>	 //计时
#include <streambuf>   //excel读取
using namespace std; 

// my
#include "rand.h"
using namespace g2orand;

//g2o 
#include <g2o/core/base_vertex.h>
#include <g2o/core/base_edge.h>
#include <g2o/core/base_binary_edge.h>
#include <g2o/core/base_unary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_dogleg.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>
#include <g2o/types/sba/types_six_dof_expmap.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/stuff/misc.h"
#include "g2o/stuff/macros.h"

// typedef
typedef Eigen::Matrix<double,6,1> Vector6d;
typedef Eigen::Matrix<double,6,6> Matrix6d;
#endif