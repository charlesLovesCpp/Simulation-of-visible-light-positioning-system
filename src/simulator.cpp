#include "myslam/simulator.h"
#include <boost/concept_check.hpp>

namespace vlp {
  namespace simuator {
    
    void generateStates( StateVector& states,  Vector6d initVector, Vector3d transNoise, Vector3d rotNoise , int num_of_state, int stepLen,  int mode) 
    {
      states.clear();
      cerr << "Simulator: sampling nodes ...";
      
      State state;
      state.id = 0;
      state.pose = SE3::exp(initVector);
      state.init_pose = state.pose;
      //初始odometry设为0
      state.odometry = SE3::exp(initVector.setZero());
      states.push_back( state );
      
      while ( (int)states.size() < num_of_state ) {
	State nextState = generateNewState(states.back(), getMotion(0, stepLen), transNoise, rotNoise);//默认走直线，由mode决定
	states.push_back(nextState);
      }
      //states[0].init_pose =  states[1].init_pose;
      cerr << "done." << endl;  
    }
    
    State generateNewState( State preState, SE3 trueMotion, Eigen::Vector3d transNoise, Eigen::Vector3d rotNoise) {
      
     Matrix6d covariance;
	  covariance.fill(0.);
	  covariance(0, 0) = transNoise[0]*transNoise[0];
	  covariance(1, 1) = transNoise[1]*transNoise[1];
	  covariance(2, 2) = transNoise[2]*transNoise[2];
	  covariance(3, 3) = rotNoise[0]*rotNoise[0];
	  covariance(4, 4) = rotNoise[1]*rotNoise[1];
	  covariance(5, 5) = rotNoise[2]*rotNoise[2];
      Matrix6d information = covariance.inverse();
      State nextState;
      //T_wc版本
      nextState.id = preState.id + 1;
      nextState.pose = preState.pose * trueMotion;
      nextState.odometry = perturb( trueMotion, transNoise, rotNoise );
      nextState.odom_information = information;
      nextState.init_pose = nextState.odometry * preState.init_pose;

      return nextState;
    }

    SE3 getMotion(int direction, double stepLen) {
      switch( direction ) {
	case 0:
	  return SE3( SO3(0, 0, 0), Eigen::Vector3d(stepLen, 0, 0));//这样写虽然是不对的，但是因为旋转均为0，所以结果相同
	  break;
	case 1:
	  return SE3( SO3(0, 0, 0), Eigen::Vector3d(0, stepLen, 0));
	  break;
	case 2:
	  return SE3( SO3(0, 0, 0), Eigen::Vector3d(0, -stepLen, 0));
	  break;
	default:
	  cerr << "Unknown motion direction" << endl;
      }
    }

    SE3 perturb( SE3 trueMotion, Eigen::Vector3d transNoise,  Eigen::Vector3d rotNoise) {
      Vector6d noiseVector;
      noiseVector << trueMotion.log()[0]+Rand::gauss_rand(0.0, transNoise[0]), 
				trueMotion.log()[1]+Rand::gauss_rand(0.0, transNoise[1]), 
				trueMotion.log()[2]+Rand::gauss_rand(0.0, transNoise[2]), 
				trueMotion.log()[3]+Rand::gauss_rand(0.0, rotNoise[0]), 
				trueMotion.log()[4]+Rand::gauss_rand(0.0, rotNoise[1]),
				trueMotion.log()[5]+Rand::gauss_rand(0.0, rotNoise[2]);//se(3)的平移在前，旋转在后
      return SE3::exp(noiseVector);
    }
    
    bool isObservable( State state, Eigen::Vector2d uv)
    {
      if( uv[0]>0 && uv[0]<650 && uv[1]>0 && uv[1]<500 ) {	
	//cout<<"测试pose "<<i<<", LED "<<j<<" : "<<uv.transpose()<<endl;
	return true;
      }
      return false;
    }
    
    int generateCameraMeas( State& state, LedLayout ledlayout, Eigen::Vector2d measNoise, myslam::Camera cam, int mode ) 
    { 
      //所有测量值的向量
      CameraMeasVector& camMeasVector = state.cam_meas;
      camMeasVector.clear();
      //可见led的数量
      int num_of_visible_led = 0;
      //测量噪声的信息矩阵
      Matrix2d covariance; covariance.fill(0.);
	    covariance(0, 0) = measNoise[0]*measNoise[0];
	    covariance(1, 1) = measNoise[1]*measNoise[1];
      Matrix2d information = covariance.inverse();
      
      for (size_t i=0; i<ledlayout.size(); ++i) {
	Vector2d uv = cam.world2pixel( ledlayout[i].xyz, state.pose.inverse());
	
	if ( mode == 1) {
	  num_of_visible_led++;
	  // 全可视
	  camMeasVector.push_back( CameraMeas() );
	  CameraMeas& camMeas = camMeasVector.back();
	  
	  camMeas.fromLed = ledlayout[i].id;
	  camMeas.uv = Vector2d(uv[0]+Rand::gauss_rand(0.0, measNoise[0]), uv[1]+Rand::gauss_rand(0.0, measNoise[1]));
	  camMeas.information = information;
	} else {
	  if ( isObservable( state, uv ) ) {
	    num_of_visible_led++;
	    //如果观测得到该led，进行测量并存入state的cam_meas中
	    camMeasVector.push_back( CameraMeas() );
	    CameraMeas& camMeas = camMeasVector.back();
	    
	    camMeas.fromLed = ledlayout[i].id;
	    camMeas.uv = Vector2d(uv[0]+Rand::gauss_rand(0.0, measNoise[0]), uv[1]+Rand::gauss_rand(0.0, measNoise[1]));
	    camMeas.information = information;
	  }
	}
	
      }
      return num_of_visible_led;
    }
    
    // 随机均匀分布生成led layout
    void generateLedLayout( int ledNum, Eigen::Vector3d& area , LedLayout& ledLayout)
    {
      ledLayout.clear();
      
      for ( size_t i=0; i<ledNum; ++i) {
	ledLayout.push_back( Led() );
	Led& led = ledLayout.back();	
	
	led.id = i;
	led.xyz = Eigen::Vector3d(Rand::uniform_rand( 0, area[0] ),
						    Rand::uniform_rand( 0, area[1] ),
						    area[2]);
      }   
    }
    
    void generateLedLayout( LedLayout& ledLayout )
    {
      ledLayout.clear();
      
      ledLayout.push_back( Led() );
      Led& led0 = ledLayout.back();	
      led0.id = 0;
      led0.xyz = Eigen::Vector3d(20, 80, 50);
      
       ledLayout.push_back( Led() );
      Led& led1 = ledLayout.back();	
      led1.id = 1;
      led1.xyz = Eigen::Vector3d(80, 20, 50);
      
      ledLayout.push_back( Led() );
      Led& led2 = ledLayout.back();	
      led2.id = 2;
      led2.xyz = Eigen::Vector3d(20, 20, 50);
      
      ledLayout.push_back( Led() );
      Led& led3 = ledLayout.back();	
      led3.id = 3;
      led3.xyz = Eigen::Vector3d(80, 80, 50);
//       
//       ledLayout.push_back( Led() );
//       Led& led4 = ledLayout.back();	
//       led4.id = 4;
//       led4.xyz = Eigen::Vector3d(50, 50, 50);
//       
//       ledLayout.push_back( Led() );
//       Led& led5 = ledLayout.back();	
//       led5.id = 5;
//       led5.xyz = Eigen::Vector3d(65, 65, 50);
//       
//       ledLayout.push_back( Led() );
//       Led& led6 = ledLayout.back();	
//       led6.id = 6;
//       led6.xyz = Eigen::Vector3d(35, 35, 50);
//       
//       ledLayout.push_back( Led() );
//       Led& led7 = ledLayout.back();	
//       led7.id = 7;
//       led7.xyz = Eigen::Vector3d(10, 50, 50);
//       
//       ledLayout.push_back( Led() );
//       Led& led8 = ledLayout.back();	
//       led8.id = 8;
//       led8.xyz = Eigen::Vector3d(80, 50, 50);
    }
    
    void generateLedLayout( Eigen::Vector2d& vector, Eigen::Vector3d& area,  LedLayout& ledLayout, int mode ) 
    {
      ledLayout.clear();
      
      Vector2d gab;
      Vector2d num;
      switch( mode ) {
	case 0:
	  gab = vector;
	  num = Vector2d( (int)area[0]/gab[0], (int)area[1]/gab[1] );
	  break;
	case 1:
	  num = vector;
	  gab = Vector2d( area[0]/num[0], area[1]/num[1] );
	  break;
	default:
	  break;
      }

      int id = 0;
      for(size_t i=0; i<num[0]; ++i) {
	for(size_t j=0; j<num[1]; ++j) {
	  ledLayout.push_back( Led() );
	  Led& led = ledLayout.back();	
	
	  led.id = id;
	  id++;
	  led.xyz = Vector3d(5+i*gab[0],5+ j*gab[1], area[2]);
	}
      }  
    }
    
    void generateInitPose( StateVector& states ) 
    {
      states[0].init_pose = states[0].pose;
      for ( size_t i=1; i<states.size(); ++i ) {
	states[i].init_pose = states[i].odometry *  states[i-1].pose;
      }
    }

    void saveSampledData( StateVector states, LedLayout ledLayout, Vector3d area) 
    { 
      ofstream fout;
      //真实姿态
      char filenametemp1[] ="truePose.txt";
      fout.open( filenametemp1 );  
      //真实坐标
      for (size_t i=0; i<states.size(); ++i) {
	SE3 truePose =states[i].pose;
	fout<<truePose.log().transpose()<<endl;
      }
      fout.flush();
      fout.close();
      //初始姿态
      char filenametemp2[] ="initialPose.txt";
      fout.open( filenametemp2 );  
      for (size_t i=0; i<states.size(); ++i) {
	SE3 initPose = states[i].init_pose;
	fout<<initPose.log().transpose()<<endl;
      }
      fout.flush();
      fout.close();
      //里程计测量值
      char filenametemp3[] ="odometry.txt";
      fout.open( filenametemp3 );  
      for (size_t i=0; i<states.size(); ++i) {
	SE3 odometry = states[i].odometry;
	fout<<odometry.log().transpose()<<endl;
      }
      
      fout.flush();
      fout.close();
      //led布局
      char filenametemp4[] ="ledLayout.txt";
      fout.open( filenametemp4 );  
      fout<<area[0]<<" "<<area[1]<<" "<<area[2]<<endl;
      for (size_t i=0; i<ledLayout.size(); ++i) {
	Vector3d xyz = ledLayout[i].xyz;
	fout<<xyz.transpose()<<endl;
      }
      
      fout.flush();
      fout.close();
      
      //对应led观测值
      char filenametemp5[] ="cameraMeas.txt";
      fout.open( filenametemp5 );  
      for (size_t i=0; i<states.size(); ++i) {
	fout<<states[i].cam_meas.size()<<" ";
	for ( size_t j=0; j<states[i].cam_meas.size(); ++j ) {
	  fout<<states[i].cam_meas[j].fromLed<<" "<<states[i].cam_meas[j].uv.transpose()<<" ";
	}
	fout<<endl;
      }
      
      fout.flush();
      fout.close();
      
    }
    
    void saveResult( vector<Vector3d> result, char filename[] )
    {
      ofstream fout;

      strcat(filename,".txt");
      fout.open( filename );  

      for (size_t i=0; i<result.size(); ++i) {
	fout<<result[i].transpose()<<endl;
      }
      fout.flush();
      fout.close();
    }
    void loadSampledData( StateVector& states, LedLayout& ledLayout, Vector3d& area)
    {
      ifstream fin_truePose, fin_initialPose, fin_odometry, fin_cameraMeas, fin_ledLayout;  
    
      fin_truePose.open("truePose.txt",ios::in);
      fin_initialPose.open("initialPose.txt",ios::in);
      fin_odometry.open("odometry.txt",ios::in);
      fin_cameraMeas.open("cameraMeas.txt",ios::in);
      fin_ledLayout.open("ledLayout.txt",ios::in);
      
      if( fin_truePose.fail() || fin_initialPose.fail() || fin_odometry.fail() || fin_cameraMeas.fail() || fin_ledLayout.fail() ) {  
	cout << "文件不存在." << endl;  
	fin_truePose.close(); 
	fin_initialPose.close(); 
	fin_odometry.close(); 
	fin_cameraMeas.close(); 
	fin_ledLayout.close(); 
	return;
      } 
      
      Vector3d transNoise(0.2, 0.2, 0.2);
      Vector3d rotNoise(DEG2RAD(2.), DEG2RAD(2.), DEG2RAD(2.));
      Vector2d measNoise(2, 2);
      
      Matrix6d covariance1;
	  covariance1.fill(0.);
	  covariance1(0, 0) = transNoise[0]*transNoise[0];
	  covariance1(1, 1) = transNoise[1]*transNoise[1];
	  covariance1(2, 2) = transNoise[2]*transNoise[2];
	  covariance1(3, 3) = rotNoise[0]*rotNoise[0];
	  covariance1(4, 4) = rotNoise[1]*rotNoise[1];
	  covariance1(5, 5) = rotNoise[2]*rotNoise[2];
      Matrix6d information1 = covariance1.inverse();
    
      Matrix2d covariance2; covariance2.fill(0.);
	    covariance2(0, 0) = measNoise[0]*measNoise[0];
	    covariance2(1, 1) = measNoise[1]*measNoise[1];
      Matrix2d information2 = covariance2.inverse();	
      
      int index_state = 0;
      while ( states.size() < 100 ) {
	double trans1, trans2, trans3, rot1, rot2, rot3;
	Vector6d true_se3, init_se3, odom_se3;
	
	fin_truePose>>trans1>>trans2>>trans3>>rot1>>rot2>>rot3;
	true_se3<< trans1, trans2, trans3, rot1, rot2, rot3;
	
	fin_initialPose>>trans1>>trans2>>trans3>>rot1>>rot2>>rot3;
	init_se3<< trans1, trans2, trans3, rot1, rot2, rot3;
	
	fin_odometry>>trans1>>trans2>>trans3>>rot1>>rot2>>rot3;
	odom_se3<< trans1, trans2, trans3, rot1, rot2, rot3;

	CameraMeasVector camMeasVector;
	int num_of_visible_led;
	fin_cameraMeas>>num_of_visible_led;
	for ( size_t i=0; i<num_of_visible_led; ++i ) {
	  int fromLed;
	  double u, v;
	  fin_cameraMeas>>fromLed>>u>>v;
	  CameraMeas camMeas;
	  camMeas.fromLed = fromLed;
	  camMeas.uv = Vector2d( u, v );
	  camMeas.information = information2;
	  
	  camMeasVector.push_back(camMeas);
	}
	
	State state;
	state.id = index_state;
	state.pose = SE3::exp(true_se3);
	state.init_pose = SE3::exp(init_se3);
	state.odometry = SE3::exp(odom_se3);
	state.odom_information = information1;
	state.cam_meas = camMeasVector;
	index_state++;
	
// 	cout<<"load到的真实姿态： "<<state.pose.translation().transpose()<<endl;
// 	cout<<"load到的odometry： "<<state.odometry.translation().transpose()<<endl;
	
	states.push_back(state);
      }
      
      int index_led = 0;
      fin_ledLayout>>area[0]>>area[1]>>area[2];
      while( !fin_ledLayout.eof() ) {
	double x, y, z;
	Vector3d xyz;
	
	fin_ledLayout>>x>>y>>z;
	
	Led led;
	led.id = index_led;
	led.xyz = Vector3d( x, y, z );
	ledLayout.push_back( led );
	index_led++;
      }
      cout<<"读取到有LED： "<<index_led<<endl;
      
      fin_truePose.close(); 
      fin_initialPose.close(); 
      fin_odometry.close(); 
      fin_cameraMeas.close(); 
      fin_ledLayout.close(); 
    }
    
    void showData( const StateVector& states, const vector<Vector3d> estimate, const LedLayout& ledLayout, Eigen::Vector3d area ) 
    {
      //背景框
      Mat plot_image( area[1]+10, area[0]+10, CV_8UC3 );  
      cv::namedWindow("Plot"); 
      for (size_t y=0; y<plot_image.rows; y++)        //遍历每一行每一列并设置其像素值  
      {  
	for (size_t x=0; x<plot_image.cols; x++)  
	{
	  unsigned char* row_ptr = plot_image.ptr<unsigned char> ( y );  // row_ptr是第y行的头指针
	  unsigned char* data_ptr = &row_ptr[ x*plot_image.channels() ]; // data_ptr 指向待访问的像素数据，这里的&是指c语言的取地址
	  for ( int c = 0; c != plot_image.channels(); c++ )
	  {
	      data_ptr[c] = 255; // data为I(x,y)第c个通道的值
	  }
	}  
      } 
      cv::rectangle( plot_image,  cv::Rect( 5, 5, area[0], area[1] ), cv::Scalar( 0, 0, 0 ), 1, 8 );  
      
      // 所有LED画下来
      for ( size_t i=0; i<ledLayout.size(); ++i) {
	cv::Point p(5+ledLayout[i].xyz[0], 5+ledLayout[i].xyz[1]);
	cv::circle(plot_image, p, 2, cv::Scalar(0, 0, 255));
      }
      
      // 可视的LED变成实心
      for ( size_t i=0; i<states.size(); ++i ) {
	for ( size_t j=0; j<states[i].cam_meas.size(); ++j ) {
	  //cout<<"观察到LED了："<<states[i].cam_meas.size()<<endl;
	  cv::Point p(5+ledLayout[states[i].cam_meas[j].fromLed].xyz[0], 5+ledLayout[states[i].cam_meas[j].fromLed].xyz[1]);
	  cv::circle(plot_image, p, 3, cv::Scalar(0, 0, 255), -1);
	}
      }
      
      //True
      cv::Point pre_point, cur_point;
      pre_point = cv::Point( 5+states[0].pose.translation()[0], 5+states[0].pose.translation()[1] );
      for (size_t i=1; i<states.size(); ++i) {
	cur_point =  cv::Point( 5+states[i].pose.translation()[0], 5+states[i].pose.translation()[1] );
	cv::line(plot_image, pre_point, cur_point, cv::Scalar(255, 0, 0));  
	pre_point = cur_point;
      }
      //Init
      pre_point = cv::Point( 5+states[0].init_pose.translation()[0], 5+states[0].init_pose.translation()[1] );
      for (size_t i=1; i<states.size(); ++i) {
	cur_point =  cv::Point( 5+states[i].init_pose.translation()[0], 5+states[i].init_pose.translation()[1] );
	cv::line(plot_image, pre_point, cur_point, cv::Scalar(0, 0, 255));  
	pre_point = cur_point;
      }
      //Estimate
      pre_point = cv::Point(5+estimate[0][0], 5+estimate[0][1]);
      for (size_t i=1; i<estimate.size(); ++i) {
	if( estimate[i][0] == 9 );
	cur_point = cv::Point(5+estimate[i][0], 5+estimate[i][1]);
	cv::line(plot_image, pre_point, cur_point, cv::Scalar(0, 255, 0));  
	pre_point = cur_point;
      }

      cv::imshow( "Plot", plot_image );  
      cv::waitKey( 0 );  
    }
    
    
    
    
    
    
    
    Simulator::Simulator()
	{
	  time_t seed = time(0);//出现随机数相同的原因是：两次计算之间相隔时间太短了，对应的seed一样
	  //srand((unsigned)seed);
	  Rand::seed_rand(static_cast<unsigned int>(seed));
	}

    Simulator::~Simulator() {}

    void Simulator::simulate(int numNodes, int numLED) {
      
      myslam::Camera cam(520.9, 521.0, 325.1, 249.7, 0); //相机的FOV为X:(-62.5,62.5), Y:(-48,48)
      Matrix<double, 3, 3> K; 
      K << 520.9, 0, 325.1, 0, 521.0, 249.7, 0,  0, 1;
	      
      //Eigen::Vector3d transNoise(0.2, 0.2, 0.1);
      Vector3d transNoise(0.2, 0.2, 0.2);
      Vector3d rotNoise(DEG2RAD(2.), DEG2RAD(2.), DEG2RAD(2.));
      Vector2d landmarkNoise(2, 2);
      
      int step = 0;
      double stepLen = 1.0;
      //以空间大小为标准设置LED
      Vector3d boundArea(1000, 500, 100);//空间总大小，原点在左下角
//       Vector2d tidyLED(20, 20);//x,y轴LED个数
//       Vector2d tidyLedDistance(boundArea[0]/tidyLED[0], boundArea[1]/tidyLED[1]);//x,y轴led间隔
      
      // 也可以选择以间隔大小来布局LED
      Vector2d tidyLedDistance(60, 40);//x,y轴led间隔,62.5和48分别是极限值
      Vector2d tidyLED(boundArea[0]/tidyLedDistance[0], boundArea[1]/tidyLedDistance[1]);//x,y轴LED个数
      
      Matrix6d covariance;
	  covariance.fill(0.);
	  covariance(0, 0) = transNoise[0]*transNoise[0];
	  covariance(1, 1) = transNoise[1]*transNoise[1];
	  covariance(2, 2) = transNoise[2]*transNoise[2];
	  covariance(3, 3) = rotNoise[0]*rotNoise[0];
	  covariance(4, 4) = rotNoise[1]*rotNoise[1];
	  covariance(5, 5) = rotNoise[2]*rotNoise[2];
      Matrix6d information = covariance.inverse();

      Simulator::PosesVector& poses = _poses;
      poses.clear();
      Simulator::LandmarksVector& landmarks = _landmarks;
      landmarks.clear();
      //2.生成第一个pose作为初始pose
      Simulator::Pose firstPose;//初始化第一个pose，为I
      firstPose.id = 0;
      Vector6d initVector;
      //出发点设置××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××
      initVector << 100, 250, 0, 0, 0, 0;//se(3)的平移在前，旋转在后
      
      firstPose.truePose = SE3::exp(initVector);
      firstPose.simulatePose = SE3::exp(initVector);
      poses.push_back(firstPose);
      cerr << "Simulator: sampling nodes ...";
      //3. 确定方向，生成全部的pose //右乘 
      SE3 trueMotion;
      while ((int)poses.size() < numNodes) {
// 	double ranDir = Rand::uniform_rand(0,1) ;
// 	if (ranDir > 0.5) {
// 	  trueMotion = getMotion(MO_STRAIGHT, stepLen);
// 	} else if (ranDir <=  0.5 && ranDir > 0.25) {
// 	  trueMotion = getMotion(MO_LEFT, stepLen);
// 	} else {
// 	  trueMotion = getMotion(MO_RIGHT, stepLen);
// 	}
// 	Simulator::Pose nextPose = generateNewPose(poses.back(), trueMotion, transNoise, rotNoise);//默认走直线
	
	Simulator::Pose nextPose = generateNewPose(poses.back(), getMotion(MO_STRAIGHT, stepLen), transNoise, rotNoise);//默认走直线
	poses.push_back(nextPose);
	step++;
      }
      cerr << "done." << endl;
      
      cerr << "Simulator: creating landmark LEDs ...";
      
      Simulator::Landmark landmark1;
      landmark1.id = 0;
      landmark1.xyzLED = Vector3d(100, 250, 100);//现象：坐标越大，精度越高
      landmarks.push_back(landmark1);
      
      Simulator::Landmark landmark2;
      landmark2.id = 1;
      landmark2.xyzLED = Vector3d(300, 250, 100);
      landmarks.push_back(landmark2);
      
      Simulator::Landmark landmark3;
      landmark3.id = 2;
      landmark3.xyzLED = Vector3d(200, 150, 100);
      landmarks.push_back(landmark3);
       
      Simulator::Landmark landmark4;
      landmark4.id = 3;
      landmark4.xyzLED = Vector3d(200, 350, 100);
      landmarks.push_back(landmark4);
     
//       Simulator::Landmark landmark5;
//       landmark5.id = 4;
//       landmark5.xyzLED = Vector3d(50, 180, 100);
//       landmarks.push_back(landmark5);
      
      //随机布局led
//       int randomLayoutLED = numLED;//400个为FOV里平均4LED的配置
//        int ledID = 0;
//       for (size_t i=0; i<randomLayoutLED; ++i) {
// 	Simulator::Landmark createLandmark;
// 	createLandmark.id = ledID;
// 	Vector3d pos(Rand::uniform_rand(0,1000), Rand::uniform_rand(0,1000), boundArea[2]);
// 	createLandmark.xyzLED = pos;//(x,y,z)
// 	landmarks.push_back(createLandmark);
// 	ledID++;
//       }
      
//       //空间为主的LED布局
//       int ledID = 0;
//       for(size_t i=0; i<tidyLED[0]; ++i) {
// 	for(size_t j=0; j<tidyLED[1]; ++j) {
// 	  Simulator::Landmark createLandmark;
// 	  createLandmark.id = ledID;
// 	  Vector3d xyzLEDmid(i*tidyLedDistance[0], j*tidyLedDistance[1], boundArea[2]);
// 	  createLandmark.xyzLED = xyzLEDmid;//(x,y,z)
// 	  //cout<<"LED id 号： "<<i*tidyLED[0]+j<<"  对应坐标: "<<xyzLEDmid.transpose()<<endl;
// 	  landmarks.push_back(createLandmark);
// 	  ledID++;
// 	}
//       }     
     
     cerr << "done." << endl;
      
      // 遍历所有LED看是否在可视范围，如果在，添加进Pose的landmarks里
      for ( size_t i = 0; i<poses.size(); ++i) {
	for (size_t j=0; j<landmarks.size(); ++j) {  
	  Vector2d uv = cam.world2pixel( landmarks[j].xyzLED, poses[i].truePose.inverse());//T_wc版本
// 	  Vector2d uv = cam.world2pixel( landmarks[j].xyzLED, poses[i].truePose.inverse());//T_cw版本
	  
	  //if( uv[0]>0 && uv[0]<650 && uv[1]>0 && uv[1]<500 ) {	    
	    cout<<"测试pose "<<i<<", LED "<<j<<" : "<<uv.transpose()<<endl;
	    poses[i].landmarks.push_back(landmarks[j]);
	  //}
	}
      }
      
      cerr << "Simulator: Adding odometry measurements ... ";
      for (size_t i = 1; i < poses.size(); ++i) {
	const Pose& prev = poses[i-1];
	const Pose& p = poses[i];

	_odometry.push_back(PoseEdge());//生成一个GridEdge,放进去，然后再取出来
	PoseEdge& edge = _odometry.back();

	//T_wc 版本
	edge.from = prev.id;
	edge.to = p.id;
	edge.trueTransf =  prev.truePose.inverse() *p.truePose;//Twc(-1)*Twc' = Tcc'，这是从相当于前面pose - 后面pose,是个所谓的负值
	//带有噪声的odometry测量值noiseTransf
	edge.noiseTransf = sampleTransformation(edge.trueTransf, transNoise,  rotNoise);
	edge.simulateTransf = prev.simulatePose.inverse() * p.simulatePose  ;//里程计的测量值不应该是这个值，而应该是真实里程值加上高斯误差，simulatePose其实只是初始值而已
	edge.information = information;
	
// 	//T_cw版本
// 	edge.from = prev.id;
// 	edge.to = p.id;
// 	edge.trueTransf =  p.truePose * prev.truePose.inverse();//Twc(-1)*Twc' = Tcc'，这是从相当于前面pose - 后面pose,是个所谓的负值
// 	//带有噪声的odometry测量值noiseTransf
// 	edge.noiseTransf = sampleTransformation(edge.trueTransf, transNoise,  rotNoise);
// 	edge.simulateTransf = p.simulatePose * prev.simulatePose.inverse() ;//里程计的测量值不应该是这个值，而应该是真实里程值加上高斯误差，simulatePose其实只是初始值而已
// 	edge.information = information;
      }
      
      cerr << "done." << endl;
      
      //创建完pose和odometry后，建立初始pose
      InitPose init_pose;
      //初始值也需要优化的情况
      Vector3d transNoise2(0, 0, 0);
       Vector6d initVector2;
       
      //初始误差固定的情况。这里是16的原因是想把误差固定起来作比较
      initVector2 << 16, 50, 0, 0, 0, 0;//se(3)的平移在前，旋转在后
      //init_pose.initPose = SE3::exp(initVector2);
      
      //初始值是带高斯误差的真实值
      SE3 noiseInit = poses[1].simulatePose;
      init_pose.initPose = noiseInit;

      //初始值是正确值的情况
      //init_pose.initPose = poses[0].truePose;
      
      init_pose.id = 0;
      _init_poses.push_back(init_pose);
      for (size_t i = 0; i < _odometry.size(); ++i) {
	InitPose init_pose;
	const PoseEdge& odom = _odometry[i];
	//init_pose.initPose =  _init_poses.back().initPose *  odom.noiseTransf ;//T_wc版本
	init_pose.initPose =  odom.noiseTransf * _init_poses.back().initPose  ;//T_cw版本 左乘造成的误差最大，因为推导得 位移为R't*t',而右乘为t'+t,受了两次误差的影响，所以误差更大
	init_pose.id = i+1;
	_init_poses.push_back(init_pose);
      }
      
      cerr << "Simulator: Adding LEDs observation  measurements ... ";
      _measurement.clear();
      
      Matrix2d covarianceLM; covarianceLM.fill(0.);
	    covarianceLM(0, 0) = landmarkNoise[0]*landmarkNoise[0];
	    covarianceLM(1, 1) = landmarkNoise[1]*landmarkNoise[1];
      Matrix2d informationLM = covarianceLM.inverse();
	    
      for (size_t i = 0; i < poses.size(); ++i) {
	const Pose& p = poses[i];
	for (size_t j=0; j<p.landmarks.size(); j++) {
	  const Landmark& landmark = p.landmarks[j];
	  
	  _measurement.push_back(LandmarkEdge());
	  LandmarkEdge& le = _measurement.back();

	  Vector2d uv2 = cam.world2pixel(landmark.xyzLED, p.truePose.inverse()); //T_wc版本
// 	  Vector2d uv2 = cam.world2pixel(landmark.xyzLED, p.truePose); //T_cw版本
	  le.from = landmark.id;						
	  le.to = p.id;
	  le.trueMeas = uv2;
	  le.simulateMeas = Vector2d(uv2[0]+Rand::gauss_rand(0.0, landmarkNoise[0]), uv2[1]+Rand::gauss_rand(0.0, landmarkNoise[1]));
	  le.information = informationLM;
	}
      }
      cerr << "done." << endl;
    }

    Simulator::Pose Simulator::generateNewPose( const Simulator::Pose& prev, const SE3 trueMotion, const Eigen::Vector3d& transNoise, const Eigen::Vector3d rotNoise) {
      Simulator::Pose next;
      next.id = prev.id + 1;
      //T_wc版本
      next.truePose =  prev.truePose * trueMotion ;
      SE3 noiseMotion = sampleTransformation(trueMotion, transNoise,  rotNoise);
      next.simulatePose =  prev.truePose * noiseMotion;
      
      //T_cw版本
//       next.truePose =  trueMotion  * prev.truePose ;
//       SE3 noiseMotion = sampleTransformation(trueMotion, transNoise,  rotNoise);
//       next.simulatePose =  noiseMotion * prev.truePose;
      
      return next;
    }

    SE3 Simulator::getMotion(int motionDirection, double stepLen) {
      switch( motionDirection ) {
	case MO_STRAIGHT:
	  return SE3( SO3(0, 0, 0), Eigen::Vector3d(stepLen, 0, 0));//这样写虽然是不对的，但是因为旋转均为0，所以结果相同
	  break;
	case MO_LEFT:
	  return SE3( SO3(0, 0, 0), Eigen::Vector3d(0, stepLen, 0));
	  break;
	case MO_RIGHT:
	  return SE3( SO3(0, 0, 0), Eigen::Vector3d(0, -stepLen, 0));
	  break;
	default:
	  cerr << "Unknown motion direction" << endl;
      }
    }

    SE3 Simulator::sampleTransformation(const SE3& trueMotion_, const Eigen::Vector3d& transNoise, const Eigen::Vector3d rotNoise) {
      Vector6d noiseVector;
      noiseVector << trueMotion_.log()[0]+Rand::gauss_rand(0.0, transNoise[0]), 
				trueMotion_.log()[1]+Rand::gauss_rand(0.0, transNoise[1]), 
				trueMotion_.log()[2]+Rand::gauss_rand(0.0, transNoise[2]), 
				trueMotion_.log()[3]+Rand::gauss_rand(0.0, rotNoise[0]), 
				trueMotion_.log()[4]+Rand::gauss_rand(0.0, rotNoise[1]),
				trueMotion_.log()[5]+Rand::gauss_rand(0.0, rotNoise[2]);//se(3)的平移在前，旋转在后
      return SE3::exp(noiseVector);
    }
  }
}

    //       //以LED为标准设置空间
//       Vector2d ledArea(10, 5);
//       Vector2d ledDistence(50, 50);//50x50 对应2~3led定位， initVector << 2*ledDistence[0], 2*ledDistence[1], 0, 0, 0, 0;
//       //Vector2d ledDistence(62, 40);//62x40 对应4led定位
      
      //LED为主的空间划分
//       int ledID = 0;
//       for(int i=0; i<ledArea[0]; ++i) {
// 	for(int j=0; j<ledArea[1]; ++j) {
// 	  Simulator::Landmark createLandmark;
// 	  createLandmark.id = ledID;
// 	  Vector3d xyzLEDmid(i*ledDistence[0], j*ledDistence[1], 100);
// 	  createLandmark.xyzLED = xyzLEDmid;//(x,y,z)
// 	  //cout<<"LED id 号： "<<i*ledArea[0]+j<<"  对应坐标: "<<xyzLEDmid.transpose()<<endl;
// 	  landmarks.push_back(createLandmark);
// 	  ledID++;
// 	}
//       }     


      //证明实验：LED空间越近，误差越大，尤其是Z轴
//       Simulator::Landmark landmark1;
//       landmark1.id = 0;
//       landmark1.xyzLED = Vector3d(100, 120, 100);//现象：坐标越大，精度越高
//       landmarks.push_back(landmark1);
//       
//       Simulator::Landmark landmark2;
//       landmark2.id = 1;
//       landmark2.xyzLED = Vector3d(110, 120, 100);
//       landmarks.push_back(landmark2);
//       
//       Simulator::Landmark landmark3;
//       landmark3.id = 2;
//       landmark3.xyzLED = Vector3d(210, 220, 100);
//       landmarks.push_back(landmark3);
//       
//       Simulator::Landmark landmark4;
//       landmark4.id = 3;
//       landmark4.xyzLED = Vector3d(230, 180, 100);
//       landmarks.push_back(landmark4);
/*      
       Simulator::Landmark landmark1;
      landmark1.id = 0;
      landmark1.xyzLED = Vector3d(100, 250, 100);//现象：坐标越大，精度越高
      landmarks.push_back(landmark1);
      
      Simulator::Landmark landmark2;
      landmark2.id = 1;
      landmark2.xyzLED = Vector3d(300, 250, 100);
      landmarks.push_back(landmark2);
      
      Simulator::Landmark landmark3;
      landmark3.id = 2;
      landmark3.xyzLED = Vector3d(200, 150, 100);
      landmarks.push_back(landmark3);
      
      Simulator::Landmark landmark4;
      landmark4.id = 3;
      landmark4.xyzLED = Vector3d(200, 350, 100);
      landmarks.push_back(landmark4);

      Simulator::Landmark landmark5;
      landmark5.id = 4;
      landmark5.xyzLED = Vector3d(50, 200, 100);
      landmarks.push_back(landmark5);
      
     Simulator::Landmark landmark6;
      landmark6.id = 5;
      landmark6.xyzLED = Vector3d(350, 300, 100);
      landmarks.push_back(landmark6);
      
      Simulator::Landmark landmark7;
      landmark7.id = 6;
      landmark7.xyzLED = Vector3d(300, 120, 100);
      landmarks.push_back(landmark7);
      
      Simulator::Landmark landmark8;
      landmark8.id = 7;
      landmark8.xyzLED = Vector3d(230, 280, 100);
      landmarks.push_back(landmark8);
      
      Simulator::Landmark landmark9;
      landmark7.id = 8;
      landmark7.xyzLED = Vector3d(270, 260, 100);
      landmarks.push_back(landmark7);
      */
//       Simulator::Landmark landmark10;
//       landmark8.id = 9;
//       landmark8.xyzLED = Vector3d(80, 110, 100);
//       landmarks.push_back(landmark8);