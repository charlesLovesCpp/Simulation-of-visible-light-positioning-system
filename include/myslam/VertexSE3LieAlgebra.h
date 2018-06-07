#ifndef VLP_G2O_VERTEXSE3_H
#define VLP_G2O_VERTEXSE3_H

#include "myslam/common_include.h"


namespace vlp {
  namespace ba {
    class VertexSE3LieAlgebra: public g2o::BaseVertex<6,SE3> {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      bool read ( istream& is ) {}
      bool write ( ostream& os ) const {}
      
      virtual void setToOriginImpl() {
	_estimate = Sophus::SE3();
      }
      virtual void oplusImpl ( const double* update ) {//在edge里定义了update的含义
// 	Sophus::SE3 up ( Sophus::SO3 ( update[3], update[4], update[5] ),//这是错的，不能用
// 				    Eigen::Vector3d ( update[0], update[1], update[2] ) );
	Vector6d upVec;
	upVec <<update[3], update[4], update[5], update[0], update[1], update[2];//1-3是旋转 而SE3 是平移在前
	SE3 up = SE3::exp(upVec);
	_estimate =   up * _estimate  ;//左乘是因为这里取左扰动
      }
    };
  } // end namespace
} // end namespace

#endif