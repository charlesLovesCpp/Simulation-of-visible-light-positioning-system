#ifndef VLP_G2O_EDGESE3PROJECTXYZ2UVPOSEONLY_H
#define VLP_G2O_EDGESE3PROJECTXYZ2UVPOSEONLY_H

#include "myslam/common_include.h"
#include "myslam/VertexSE3LieAlgebra.h"
#include "camera.h"

namespace vlp {
  namespace ba {
       class EdgeSE3ProjectXYZ2UVPoseOnly: public g2o::BaseUnaryEdge<2, Eigen::Vector2d, vlp::ba::VertexSE3LieAlgebra >
    {
    public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	virtual void computeError() {
	  const VertexSE3LieAlgebra* pose = static_cast<const VertexSE3LieAlgebra*> ( _vertices[0] );
	  _error = _measurement - camera_->world2pixel ( point_, pose->estimate().inverse());//T_wc版本
// 	  _error = _measurement - camera_->world2pixel ( point_, pose->estimate());//T_cw版本
	}
	virtual void linearizeOplus() {
	  VertexSE3LieAlgebra* pose = static_cast<VertexSE3LieAlgebra*> ( _vertices[0] );
	  SE3 T_w_c ( pose->estimate() );
	  
// 	  T_wc版本
	  SE3 T_c_w = T_w_c.inverse();
	  
	  Vector3d xyz_trans = T_c_w * point_ ;
	  double x = xyz_trans[0];
	  double y = xyz_trans[1];
	  double z = xyz_trans[2];
	  double z_2 = z*z;
	  
	  //T_wc版本 
	  Matrix<double, 2,3 >de_p;
	  Matrix<double, 3,6 >dp_t;
	  
	  de_p( 0,0 ) = -1. / z * camera_->fx_;
	  de_p( 0,1 ) = 0;
	  de_p( 0,2 ) = x / z_2 * camera_->fx_ ;
	  de_p( 1,0 ) = 0;
	  de_p( 1,1 ) = -1. / z * camera_->fy_;
	  de_p( 1,2 ) = y  / z_2 * camera_->fy_;
	  

	  
	  //dp_t无J格式
	  Vector3d t = T_w_c.translation();
	  Matrix3d R_t = T_w_c.rotation_matrix().transpose();
	   
	  dp_t.fill(0); // [R_T * p_hat, -R_T ]
	  dp_t.block(0, 0, 3, 3) = R_t * SO3::hat(point_);          // P(i+1 : i+rows, j+1 : j+cols)
	  dp_t.block(0, 3, 3, 3) = -R_t;

// 	  //dp_t有J格式 -( (T_(-1)*P)_hat ) * J_l 结果是 -Tp_circle * J_l
// 	  Matrix6d J;//左Jacobi矩阵，公式在P234,7.85b
// 	  J.block(0,0,3,3) = SO3::hat(T_w_c.so3().log()); //由于J是对-hat求的
// 	  J.block(0,3,3,3) = SO3::hat(T_w_c.translation());
// 	  J.block(3,0,3,3) = Eigen::Matrix3d::Zero(3,3);
// 	  J.block(3,3,3,3) = SO3::hat(T_w_c.so3().log());
// 	  J =  Matrix6d::Identity() - J*0.5;
// 	  
// 	 Eigen::Matrix<double,4,6> Tp_circle;
// 	 Tp_circle.block(0,0,3,3) = Eigen::Matrix3d::Identity(3,3);
// 	 Tp_circle.block(0,3,3,3) = - SO3::hat(xyz_trans);
// 	 Tp_circle.block(3,0,1,3) = Eigen::Matrix3d::Zero(1,3);
// 	 Tp_circle.block(3,3,1,3) = Eigen::Matrix3d::Zero(1,3);
// 	  
// 	 dp_t =( - Tp_circle * J).block(0,0,3,6);
	  
	  
	  _jacobianOplusXi = de_p * dp_t;
	  
	  //T_cw版本
// 	  _jacobianOplusXi ( 0,0 ) =  x*y/z_2 *camera_->fx_;
// 	  _jacobianOplusXi ( 0,1 ) = - ( 1+ ( x*x/z_2 ) ) *camera_->fx_;
// 	  _jacobianOplusXi ( 0,2 ) = y/z * camera_->fx_;
// 	  _jacobianOplusXi ( 0,3 ) = -1./z * camera_->fx_;
// 	  _jacobianOplusXi ( 0,4 ) = 0;
// 	  _jacobianOplusXi ( 0,5 ) = x/z_2 * camera_->fx_;
// 
// 	  _jacobianOplusXi ( 1,0 ) = ( 1+y*y/z_2 ) *camera_->fy_;
// 	  _jacobianOplusXi ( 1,1 ) = -x*y/z_2 *camera_->fy_;
// 	  _jacobianOplusXi ( 1,2 ) = -x/z *camera_->fy_;
// 	  _jacobianOplusXi ( 1,3 ) = 0;
// 	  _jacobianOplusXi ( 1,4 ) = -1./z *camera_->fy_;
// 	  _jacobianOplusXi ( 1,5 ) = y/z_2 *camera_->fy_;
	}
	
	virtual bool read( std::istream& in ){}
	virtual bool write(std::ostream& os) const {};
	
	Vector3d point_;
	myslam::Camera* camera_;
    };
  }
}
#endif