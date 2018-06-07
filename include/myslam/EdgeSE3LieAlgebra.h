#ifndef VLP_G2O_EDGESE3_H
#define VLP_G2O_EDGESE3_H

#include "myslam/common_include.h"
#include "myslam/VertexSE3LieAlgebra.h"

namespace vlp {
  namespace ba {
    // 给定误差求J_R^{-1}的近似
    Matrix6d JRInv( SE3 e )
    {
	Matrix6d J;
	J.block(0,0,3,3) = SO3::hat(e.so3().log()); //(11.10的矩阵I+1/2×[...])
	J.block(0,3,3,3) = SO3::hat(e.translation());
	J.block(3,0,3,3) = Eigen::Matrix3d::Zero(3,3);
	J.block(3,3,3,3) = SO3::hat(e.so3().log());
	J = J*0.5 + Matrix6d::Identity();
	return J;
    }

    class EdgeSE3LieAlgebra: public g2o::BaseBinaryEdge<6, SE3, VertexSE3LieAlgebra, VertexSE3LieAlgebra> {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      
      bool read ( istream& is) {}
      bool write ( ostream& os) const {}
      virtual void computeError() {
	Sophus::SE3 v1 = ( static_cast<VertexSE3LieAlgebra*>(_vertices[0]) )->estimate();
	Sophus::SE3 v2 = ( static_cast<VertexSE3LieAlgebra*>(_vertices[1]) )->estimate();
	_error = ( _measurement.inverse() * v1.inverse() * v2).log();
      }
      virtual void linearizeOplus() {
	Sophus::SE3 v1 = ( static_cast<VertexSE3LieAlgebra*>(_vertices[0]) )->estimate();
	Sophus::SE3 v2 = ( static_cast<VertexSE3LieAlgebra*>(_vertices[1]) )->estimate();
	Matrix6d J = JRInv( SE3::exp( _error ) );
	_jacobianOplusXi = -  J*v2.inverse().Adj();//inverse掉的原因是因为作者的Vertex是T_cw,而书本里的推导是T_wc
	_jacobianOplusXj =  J*v2.inverse().Adj();
      }
    };
    
 
  } // end namespace
} // end namespace

#endif