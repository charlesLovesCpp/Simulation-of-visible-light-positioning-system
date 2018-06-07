#ifndef VLP_G2O_SIMULATOR_H
#define VLP_G2O_SIMULATOR_H

#include "myslam/common_include.h"
#include "myslam/camera.h"
#include <boost/concept_check.hpp>

namespace vlp {
  namespace simuator {
    class Simulator {
    public:
      
      enum  MotionType {
      MO_STRAIGHT, MO_LEFT, MO_RIGHT,
      MO_NUM_ELEMS
      };
      
      /**
	* \brief simulated LED landmarks
	*/	
      struct Landmark {
	int id;
	Eigen::Vector3d xyzLED;
      };
      typedef std::vector<Landmark, Eigen::aligned_allocator<Landmark>> LandmarksVector;
      
      /**
	* \brief simulated pose
	*/	    
      struct Pose {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	int id;
	SE3 truePose;
	SE3 simulatePose;
	
	LandmarksVector landmarks;
      };
      typedef std::vector<Pose, Eigen::aligned_allocator<Pose> >  PosesVector;
      
      /**
	* \brief simulated init pose for BA
	*/	
      struct InitPose {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	int id;
	SE3 initPose;
      };
      typedef std::vector<InitPose, Eigen::aligned_allocator<InitPose> >  InitPosesVector;
      
      /**
	* \brief simulated odometry
	*/	
      struct PoseEdge {
	int from;
	int to;
	SE3 trueTransf;
	SE3 simulateTransf;
	SE3 noiseTransf;
	Matrix6d information;
      };
      typedef std::vector<PoseEdge, Eigen::aligned_allocator<PoseEdge>> PoseEdgeVector;
      
      /**
	* \brief simulated observation
	*/	
      struct LandmarkEdge {
	int from;
	int to;
	Eigen::Vector2d trueMeas;
	Eigen::Vector2d simulateMeas;
	Eigen::Matrix2d information;
      };
      typedef std::vector<LandmarkEdge, Eigen::aligned_allocator<LandmarkEdge>> LandmarkEdgeVector;
      
    public:
      Simulator();
      ~Simulator();
      
      void simulate(int numPoses, int numLED);

      const PosesVector& poses() const { return _poses;}
      const LandmarksVector& landmarks() const { return _landmarks;}
      const PoseEdgeVector& odometry() const { return _odometry;}
      const LandmarkEdgeVector& landmarkObservations() const { return _measurement;}
      
      const InitPosesVector& init_poses() const { return _init_poses;}
    protected:
      PosesVector _poses;
      LandmarksVector _landmarks;
      PoseEdgeVector _odometry;
      LandmarkEdgeVector _measurement;
      
      InitPosesVector _init_poses;
      
      Pose generateNewPose( const Pose& prev, const SE3 trueMotion, const Eigen::Vector3d& transNoise, const Eigen::Vector3d rotNoise);
      SE3 getMotion(int motionDirection, double stepLen);
      SE3 sampleTransformation(const SE3& trueMotion_, const Eigen::Vector3d& transNoise, const Eigen::Vector3d rotNoise);
    };
    
    struct CameraMeas {
	int fromLed;
	Eigen::Vector2d uv;
	Eigen::Matrix2d information;
     };
    typedef std::vector<CameraMeas, Eigen::aligned_allocator<CameraMeas>> CameraMeasVector;
    
    struct State {
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      int id;
      SE3 pose;
      SE3 init_pose;
      SE3 odometry;
      Matrix6d odom_information;
      CameraMeasVector cam_meas;
      
      Simulator::LandmarkEdgeVector observations;;
    };
    typedef std::list<State, Eigen::aligned_allocator<State>> StateList;
    typedef std::vector<State, Eigen::aligned_allocator<State>> StateVector;
    
    struct Led {
      int id;
      Eigen::Vector3d xyz;
    };
    typedef std::vector<Led, Eigen::aligned_allocator<Led>> LedLayout;
    
    // 生成一组测量值；给定真实观测位置和Ledlayout； mode=1: 全可视 mode=0：FOV内可视
    int generateCameraMeas( State& state, LedLayout Ledlayout, Eigen::Vector2d measNoise, myslam::Camera cam, int mode);
    // 生成一组led；给定led总数量，以均匀分布生成led; 
    void generateLedLayout( int ledNum, Eigen::Vector3d& area , LedLayout& ledLayout);
    // 生成一组led；mode=1: 按间隔大小 mode=0：按led每轴个数
    void generateLedLayout( Eigen::Vector2d& vector, Eigen::Vector3d& area,  LedLayout& ledLayout, int mode );
    // 生成估计初始值；初始值为第一个状态对odometry积分
    void generateInitPose( StateVector& states );
    // 检测该led在给定状态是否可观测
    bool isObservable( State state, Eigen::Vector2d uv);
    
    void generateStates( StateVector& states,  Vector6d initVector, Vector3d transNoise, Vector3d rotNoise , int num_of_state, int stepLen,  int mode);
    State generateNewState( State preState, SE3 trueMotion, Eigen::Vector3d transNoise, Eigen::Vector3d rotNoise);
    SE3 getMotion(int direction, double stepLen);
    SE3 perturb( SE3 trueMotion, Eigen::Vector3d transNoise,  Eigen::Vector3d rotNoise);
    void saveSampledData( StateVector states, LedLayout ledLayout, Vector3d area );
    void saveResult( vector<Vector3d> result, char filename[] );
    void loadSampledData( StateVector& states, LedLayout& ledLayout , Vector3d& area );
    void showData( const StateVector& states, const vector<Vector3d> stateWindowVector, const LedLayout& ledLayout, Eigen::Vector3d area );
    void generateLedLayout( LedLayout& ledLayout );
    
  } // end namespace
} // end namespace

#endif

