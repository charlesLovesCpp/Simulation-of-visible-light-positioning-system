#include <iostream>

#include "myslam/common_include.h"
#include "myslam/config.h"
#include "myslam/simulator.h"
#include "myslam/camera.h"
#include "myslam/g2o_types.h"
#include "myslam/EdgeSE3LieAlgebra.h"
#include "myslam/VertexSE3LieAlgebra.h"
#include "myslam/EdgeSE3ProjectXYZ2UVPoseOnly.h"
#include "myslam/config.h"

using namespace vlp::simuator;
using namespace vlp::ba;

//统计误差
void RMS(vector<double>&, double&, double&);
void calError_xyz( StateVector states, vector<Vector3d> stateWindowVector,  double& err);
void calError_trans( StateVector states, vector<Vector3d> stateVector,  vector<double>& err_trans );
void calError_xy( StateVector states, vector<Vector3d> stateVector,  vector<double>& err_xy, vector<double>& err_z );
void optimizeBA(StateList&, LedLayout, int); 
void SWM(StateVector states, LedLayout ledLayout, int num_of_iter, int win_size, vector<Vector3d>& trans, vector<Vector3d>& rot, double& time);
void FPM(StateVector states, LedLayout ledLayout, vector<Vector3d>& trans, vector<Vector3d>& rot, double& time);
//window基础功能实现
void exp_window_implement();
void exp_window_size_vs_accuracy();
void exp_window_size_vs_robustness();
void exp_route_length_vs_accuracy();
void exp_window_vs_pnp_layout();
void exp_window_vs_pnp_random_layout();
void exp_window_vs_pnp_reguLayout_limFOV();
//***********************************
void exp_reguLayout_unlimFOV();
void exp_irreguLayout_unlimFOV();
void exp_reguLayout_limFOV();

void exp_winsize_accuracy();
void exp_winsize_robustness();
//***********************************
//下面是要被更新的
void batchOptimize(Simulator& simulator, int numIter, vector<Vector3d>& trans, double& time);
void pnpMethod(Simulator& simulator, vector<Vector3d>& statePnPVector, double& time);
void exp1_var();
void exp4_var();

void computeError( StateVector states, vector<Vector3d> trans,  vector<Vector3d> rot, vector<double>& err_trans, vector<double>& err_rot );
void computeRMSE( vector<vector<double>> err_SWM_trans, vector<double>& RMSE_SWM_trans);

//全局变量
myslam::Camera cam(520.9, 521.0, 325.1, 249.7, 0);
Mat K = ( cv::Mat_<double> ( 3,3 ) << 520.9, 0, 325.1, 0, 521.0, 249.7, 0, 0, 1 );
      
  
int main ( int argc, char** argv ) {
  cout<<"the program starts"<<endl;
  chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
  
  time_t seed = time(0);//出现随机数相同的原因是：两次计算之间相隔时间太短了，对应的seed一样
  Rand::seed_rand(static_cast<unsigned int>(seed));

  //**************************************************************
  
  exp_window_implement();
//   exp_reguLayout_unlimFOV();
//   exp_irreguLayout_unlimFOV();
//   exp_reguLayout_limFOV();	
//   exp_winsize_accuracy();
//   exp_winsize_robustness();
  
  //exp_route_length_vs_accuracy();
  //exp_window_size_vs_robustness();

  //exp_window_vs_pnp_layout();
  //exp_window_vs_pnp_random_layout();
//   exp_window_vs_pnp_reguLayout_limFOV();
  //**************************************************************
  
  chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
  chrono::duration<double> time_used = chrono::duration_cast<chrono::duration<double>>( t2-t1 );
  double time = time_used.count();
  cout<<"运行时间："<<time<<endl;
  
  cout<<"the program ends"<<endl;
  return 0;
}

void exp_window_implement()
{
  StateVector states;
  int num_of_state = 80;
  Vector6d initVector;
  initVector << 10, 50, 0, 0, 0, 0;
  Vector3d transNoise(0.2, 0.2, 0.2);
  Vector3d rotNoise(DEG2RAD(2.), DEG2RAD(2.), DEG2RAD(2.));
  Vector2d measNoise(2, 2);
  double stepLen = 1.0;

  Vector3d area(100, 100, 50);
  int num_of_led = 4;
  LedLayout ledLayout;
  
  //生成状态
  generateStates( states, initVector, transNoise, rotNoise , num_of_state, stepLen, 0);
  //生成led布局
  //generateLedLayout( num_of_led,  area , ledLayout);
  // 生成一组led；mode=0: 按间隔大小 mode=1：按led每轴个数
  //Vector2d vec(60, 40);
  //generateLedLayout( vec, area, ledLayout, 0 );
  //生成固定LED布局
  generateLedLayout( ledLayout );
  //生成led观测  mode=1: 全可视 mode=0：FOV内可视
  for ( size_t i=0; i<states.size(); ++i ) {
    int num_of_visible_led = generateCameraMeas( states[i],  ledLayout, measNoise, cam, 1);
    cout<<"状态： "<<states[i].pose.translation().transpose();
    cout<<" 可见到 "<<num_of_visible_led<<"个LED";
    cout<<" 初始状态： "<<states[i].init_pose.translation().transpose();
    cout<<" 里程计： "<<states[i].odometry.translation().transpose()<<endl;
  }
  //保存数据
  saveSampledData( states, ledLayout, area);
  
  //Window Optimize
  vector<Vector3d> est_SWM_trans, est_SWM_rot;
  double time_window;
  SWM(states, ledLayout, 10,  10, est_SWM_trans, est_SWM_rot, time_window);
  
   //PnP Method
  vector<Vector3d> est_FPM_trans, est_FPM_rot;
  double timePnP;
  FPM(states, ledLayout, est_FPM_trans, est_FPM_rot, timePnP);
  
  //计算每个frame上的translation和rotation误差
  vector<double> err_SWM_trans, err_SWM_rot;
  computeError( states, est_SWM_trans, est_SWM_rot, err_SWM_trans, err_SWM_rot );
  
  vector<double> err_FPM_trans, err_FPM_rot;
  computeError( states, est_FPM_trans, est_FPM_rot, err_FPM_trans, err_FPM_rot );

  double sum_SWM_trans=accumulate(err_SWM_trans.begin(),err_SWM_trans.end(),0.0);
  double sum_SWM_rot=accumulate(err_SWM_rot.begin(),err_SWM_rot.end(),0.0);
  cout<<"sum_SWM_trans: "<<sum_SWM_trans<<" sum_SWM_rot: "<<sum_SWM_rot<<endl;
  
  double sum_FPM_trans=accumulate(err_FPM_trans.begin(),err_FPM_trans.end(),0.0);
  double sum_FPM_rot=accumulate(err_FPM_rot.begin(),err_FPM_rot.end(),0.0);
  cout<<"sum_FPM_trans: "<<sum_FPM_trans<<" sum_FPM_rot: "<<sum_FPM_rot<<endl;
  
  char filename_SWM_trans[20] ="SWM_trans";
  char filename_SWM_rot[20] ="SWM_rot";
  char filename_FPM_trans[20] ="FPM_trans";
  char filename_FPM_rot[20] ="FPM_rot";
  
  saveResult( est_SWM_trans, filename_SWM_trans );
  saveResult( est_SWM_rot, filename_SWM_rot );
  
  saveResult( est_FPM_trans, filename_FPM_trans );
  saveResult( est_FPM_rot, filename_FPM_rot );
  
  showData( states, est_SWM_trans, ledLayout, area );
  showData( states, est_FPM_trans, ledLayout, area );
  //作图
  
}

void computeError( StateVector states, vector<Vector3d> trans,  vector<Vector3d> rot, vector<double>& err_trans, vector<double>& err_rot ) 
{
  err_trans.clear(); err_rot.clear();
  for (size_t i=0;  i<states.size(); ++i) {
    Vector3d temp_trans( states[i].pose.translation()[0] - trans[i][0],
				      states[i].pose.translation()[1] - trans[i][1],  
				      states[i].pose.translation()[2] - trans[i][2]);
    Vector3d temp_rot( RAD2DEG(states[i].pose.log()[3]) - rot[i][0],
				  RAD2DEG(states[i].pose.log()[4]) - rot[i][1],  
				  RAD2DEG(states[i].pose.log()[5]) - rot[i][2]);
    
    err_trans.push_back( sqrt( (temp_trans[0]*temp_trans[0] + temp_trans[1]*temp_trans[1] + temp_trans[2]*temp_trans[2]) ) ); 
    err_rot.push_back( sqrt( (temp_rot[0]*temp_rot[0] + temp_rot[1]*temp_rot[1] + temp_rot[2]*temp_rot[2]) ) ); 
  }  
}

void exp_reguLayout_unlimFOV()
{
  ofstream fout;
  char filenametemp[]="exp_reguLayout_unlimFOV.txt";
  fout.open( filenametemp );  

  int num_of_state = 80;
  Vector6d initVector;
  initVector << 10, 50, 0, 0, 0, 0;
  Vector3d transNoise(0.2, 0.2, 0.2);
  Vector3d rotNoise(DEG2RAD(2.), DEG2RAD(2.), DEG2RAD(2.));
  Vector2d measNoise(2, 2);
  double stepLen = 1.0;
  Vector3d area(100, 100, 50);
  
  LedLayout ledLayout;
  generateLedLayout( ledLayout );//固定layout
  
  size_t win_min = 2;
  size_t win_max = num_of_state;
  size_t monte = 100;
  
  vector< vector<double> > err_SWM_trans;//存放各个frame对应的总路程误差
  vector< vector<double> > err_SWM_rot;
  vector< vector<double> > err_FPM_trans;//存放各个frame对应的总路程误差
  vector< vector<double> > err_FPM_rot;
  err_SWM_trans.resize( monte );
  err_SWM_rot.resize( monte );
  err_FPM_trans.resize( monte );
  err_FPM_rot.resize( monte );
  for( size_t iter=0; iter<monte; ++iter ) {
    //生成一组新的数据
    StateVector states;
    generateStates( states, initVector, transNoise, rotNoise , num_of_state, stepLen, 0);
    for ( size_t i=0; i<states.size(); ++i ) {
      int num_of_visible_led = generateCameraMeas( states[i],  ledLayout, measNoise, cam, 1);//mode=1: 全可视 mode=0：FOV内可视
    }
    //生成结束
    //Window Optimize
    vector<Vector3d> est_SWM_trans, est_SWM_rot;
    double time_window;
    SWM(states, ledLayout, 10,  10, est_SWM_trans, est_SWM_rot, time_window);
    
    //PnP Method
    vector<Vector3d> est_FPM_trans, est_FPM_rot;
    double timePnP;
    FPM(states, ledLayout, est_FPM_trans, est_FPM_rot, timePnP);

    //计算每个frame上的translation和rotation误差
    computeError( states, est_SWM_trans, est_SWM_rot, err_SWM_trans[iter], err_SWM_rot[iter] );
    computeError( states, est_FPM_trans, est_FPM_rot, err_FPM_trans[iter], err_FPM_rot[iter] );
  }//end monte
  
  vector<double> RMSE_SWM_trans, RMSE_SWM_rot;
  computeRMSE(err_SWM_trans, RMSE_SWM_trans);
  computeRMSE(err_SWM_rot, RMSE_SWM_rot);
  
  vector<double> RMSE_FPM_trans, RMSE_FPM_rot;
  computeRMSE(err_FPM_trans, RMSE_FPM_trans);
  computeRMSE(err_FPM_rot, RMSE_FPM_rot);
  
  for( size_t i=0; i<RMSE_SWM_trans.size(); ++i) {
    fout<<i<<" "<<RMSE_SWM_trans[i]<<" "<<RMSE_SWM_rot[i]<<" "<<RMSE_FPM_trans[i]<<" "<<RMSE_FPM_rot[i]<<endl;
//     fout<<i<<" "<<RMSE_SWM_trans[i]<<" "<<RMSE_SWM_rot[i]<<endl;
  }
  
  fout.flush();
  fout.close();
}

void computeRMSE( vector<vector<double>> err_SWM_trans, vector<double>& RMSE_SWM_trans) 
{
  RMSE_SWM_trans.clear();
  for( size_t i=0; i<err_SWM_trans[0].size(); ++i) {
    double temp = 0;
    for( size_t j=0; j<err_SWM_trans.size(); ++j) {
      temp += pow( err_SWM_trans[j][i], 2 );
    }
    RMSE_SWM_trans.push_back( sqrt( temp/err_SWM_trans.size() ) );
  }
}

void exp_irreguLayout_unlimFOV()
{
  ofstream fout;
  char filenametemp[]="exp_irreguLayout_unlimFOV.txt";
  fout.open( filenametemp );  

  int num_of_state = 80;
  Vector6d initVector;
  initVector << 10, 50, 0, 0, 0, 0;
  Vector3d transNoise(0.2, 0.2, 0.2);
  Vector3d rotNoise(DEG2RAD(2.), DEG2RAD(2.), DEG2RAD(2.));
  Vector2d measNoise(2, 2);
  double stepLen = 1.0;
  Vector3d area(100, 100, 50);
  
  size_t win_min = 2;
  size_t win_max = num_of_state;
  size_t monte = 100;
  
  vector< vector<double> > err_SWM_trans;//存放各个frame对应的总路程误差
  vector< vector<double> > err_SWM_rot;
  vector< vector<double> > err_FPM_trans;//存放各个frame对应的总路程误差
  vector< vector<double> > err_FPM_rot;
  err_SWM_trans.resize( monte );
  err_SWM_rot.resize( monte );
  err_FPM_trans.resize( monte );
  err_FPM_rot.resize( monte );
  for( size_t iter=0; iter<monte; ++iter ) {
    //生成一组新的数据
    StateVector states;
    LedLayout ledLayout;
    int num_of_led = 6;
    generateStates( states, initVector, transNoise, rotNoise , num_of_state, stepLen, 0);
    generateLedLayout( num_of_led,  area , ledLayout);
    for ( size_t i=0; i<states.size(); ++i ) {
      int num_of_visible_led = generateCameraMeas( states[i],  ledLayout, measNoise, cam, 1);//mode=1: 全可视 mode=0：FOV内可视
    }
    //生成结束
    //Window Optimize
    vector<Vector3d> est_SWM_trans, est_SWM_rot;
    double time_window;
    SWM(states, ledLayout, 10,  10, est_SWM_trans, est_SWM_rot, time_window);
    
    //PnP Method
    vector<Vector3d> est_FPM_trans, est_FPM_rot;
    double timePnP;
    FPM(states, ledLayout, est_FPM_trans, est_FPM_rot, timePnP);

    //计算每个frame上的translation和rotation误差
    computeError( states, est_SWM_trans, est_SWM_rot, err_SWM_trans[iter], err_SWM_rot[iter] );
    computeError( states, est_FPM_trans, est_FPM_rot, err_FPM_trans[iter], err_FPM_rot[iter] );
  }//end monte
  
  vector<double> RMSE_SWM_trans, RMSE_SWM_rot;
  computeRMSE(err_SWM_trans, RMSE_SWM_trans);
  computeRMSE(err_SWM_rot, RMSE_SWM_rot);
  
  vector<double> RMSE_FPM_trans, RMSE_FPM_rot;
  computeRMSE(err_FPM_trans, RMSE_FPM_trans);
  computeRMSE(err_FPM_rot, RMSE_FPM_rot);
  
  for( size_t i=0; i<RMSE_SWM_trans.size(); ++i) {
    fout<<i<<" "<<RMSE_SWM_trans[i]<<" "<<RMSE_SWM_rot[i]<<" "<<RMSE_FPM_trans[i]<<" "<<RMSE_FPM_rot[i]<<endl;
//     fout<<i<<" "<<RMSE_SWM_trans[i]<<" "<<RMSE_SWM_rot[i]<<endl;
  }
  
  fout.flush();
  fout.close();
}

void exp_reguLayout_limFOV() {
  ofstream fout;
  char filenametemp[]="exp_reguLayout_limFOV.txt";
  fout.open( filenametemp );  

  int num_of_state = 80;
  Vector6d initVector;
  initVector << 10, 50, 0, 0, 0, 0;
  Vector3d transNoise(0.2, 0.2, 0.2);
  Vector3d rotNoise(DEG2RAD(2.), DEG2RAD(2.), DEG2RAD(2.));
  Vector2d measNoise(2, 2);
  double stepLen = 1.0;
  Vector3d area(100, 100, 50);
  Vector2d gab(30, 30);
  
  LedLayout ledLayout;
  generateLedLayout( gab,  area,  ledLayout, 0);
  
  size_t win_min = 2;
  size_t win_max = num_of_state;
  size_t monte = 1;
  
  vector< vector<double> > err_SWM_trans;//存放各个frame对应的总路程误差
  vector< vector<double> > err_SWM_rot;
  vector< vector<double> > err_FPM_trans;//存放各个frame对应的总路程误差
  vector< vector<double> > err_FPM_rot;
  err_SWM_trans.resize( monte );
  err_SWM_rot.resize( monte );
  err_FPM_trans.resize( monte );
  err_FPM_rot.resize( monte );
  
  int flag = 0;
  
  for( size_t iter=0; iter<monte; ++iter ) {
    //生成一组新的数据
    StateVector states;
    generateStates( states, initVector, transNoise, rotNoise , num_of_state, stepLen, 0);
    for ( size_t i=0; i<states.size(); ++i ) {
      int num_of_visible_led = generateCameraMeas( states[i],  ledLayout, measNoise, cam, 0);//mode=1: 全可视 mode=0：FOV内可视
      cout<<"观测到："<<num_of_visible_led<<endl;
    }
    //生成结束
    //Window Optimize
    vector<Vector3d> est_SWM_trans, est_SWM_rot;
    double time_window;
    SWM(states, ledLayout, 10,  10, est_SWM_trans, est_SWM_rot, time_window);
    
    //PnP Method
    vector<Vector3d> est_FPM_trans, est_FPM_rot;
    double timePnP;
    FPM(states, ledLayout, est_FPM_trans, est_FPM_rot, timePnP);

    //计算每个frame上的translation和rotation误差
    computeError( states, est_SWM_trans, est_SWM_rot, err_SWM_trans[iter], err_SWM_rot[iter] );
    computeError( states, est_FPM_trans, est_FPM_rot, err_FPM_trans[iter], err_FPM_rot[iter] );
    
    //只记录一次数据，目的只是为了ledLayout
    if(flag == 0)  {
      saveSampledData( states, ledLayout, area);
      flag = 1;
    }  
  }//end monte
  vector<double> RMSE_SWM_trans, RMSE_SWM_rot;
  computeRMSE(err_SWM_trans, RMSE_SWM_trans);
  computeRMSE(err_SWM_rot, RMSE_SWM_rot);
  
  vector<double> RMSE_FPM_trans, RMSE_FPM_rot;
  computeRMSE(err_FPM_trans, RMSE_FPM_trans);
  computeRMSE(err_FPM_rot, RMSE_FPM_rot);
  
  for( size_t i=0; i<RMSE_SWM_trans.size(); ++i) {
    fout<<i<<" "<<RMSE_SWM_trans[i]<<" "<<RMSE_SWM_rot[i]<<" "<<RMSE_FPM_trans[i]<<" "<<RMSE_FPM_rot[i]<<endl;
//     fout<<i<<" "<<RMSE_SWM_trans[i]<<" "<<RMSE_SWM_rot[i]<<endl;
  }  
  
  
  
  fout.flush();
  fout.close();
}

void exp_winsize_accuracy()
{
  ofstream fout;
  char filenametemp[]="exp_winsize_accuracy.txt";
  fout.open( filenametemp );  

 
  int num_of_state = 50;
  Vector6d initVector;
  initVector << 25, 50, 0, 0, 0, 0;
  Vector3d transNoise(0.2, 0.2, 0.2);
  Vector3d rotNoise(DEG2RAD(2.), DEG2RAD(2.), DEG2RAD(2.));
  Vector2d measNoise(2, 2);
  double stepLen = 1.0;
  Vector3d area(100, 100, 50);
  
  LedLayout ledLayout;
  generateLedLayout( ledLayout );//固定layout
  
  size_t win_min = 2;
  size_t win_max = num_of_state;
  size_t monte = 100;
  
  vector< vector<double> > err_winsize;//各个window size对应的平均总路程RMSE
  err_winsize.resize( monte );
  
  for( size_t iter=0; iter<monte; ++iter ) {
    //生成一组新的数据
    StateVector states;
    generateStates( states, initVector, transNoise, rotNoise , num_of_state, stepLen, 0);
    for ( size_t i=0; i<states.size(); ++i ) {
      int num_of_visible_led = generateCameraMeas( states[i],  ledLayout, measNoise, cam, 1);//mode=1: 全可视 mode=0：FOV内可视
    }
    //生成结束
    err_winsize[iter].resize( win_max );

    for (size_t win_iter=win_min; win_iter < win_max; ++win_iter) {
      //Window Optimize
      vector<Vector3d> est_SWM_trans, est_SWM_rot;
      double time_window;
      SWM(states, ledLayout, 10,  win_iter, est_SWM_trans, est_SWM_rot, time_window);
     
      vector<double> temp_err_trans, temp_err_rot;
      computeError( states, est_SWM_trans, est_SWM_rot, temp_err_trans, temp_err_rot );// 计算一次monte试验中，各个路程上的误差大小，接下来将他们求和
      
      double sum_SWM_trans=accumulate(temp_err_trans.begin(),temp_err_trans.end(),0.0);
      double sum_SWM_rot=accumulate(temp_err_rot.begin(),temp_err_rot.end(),0.0);
      err_winsize[iter][win_iter] = sum_SWM_trans;
      
//       cout<<" sum_SWM_trans: "<<sum_SWM_trans<<" sum_SWM_rot: "<<sum_SWM_rot<<endl;
    }//end win_iter
  }//end monte

  vector<double>  rmse;
  rmse.resize(win_max);
  for( size_t i=win_min; i<win_max; ++i ) {
    double temp = 0;
    for ( size_t j=0; j<monte; ++j ) {
      temp += err_winsize[j][i];
    }
//     cout<<"temp是"<<temp<<endl;
    rmse[i] = temp / monte;
    fout<<i<<" "<<rmse[i]<<endl;
  }	
  
  fout.flush();
  fout.close();
}

void exp_winsize_robustness()
{
  ofstream fout;
  char filenametemp[]="exp_winsize_robustness.txt";
  fout.open( filenametemp );  

  int num_of_state = 50;
  Vector6d initVector;
  initVector << 25, 50, 0, 0, 0, 0;
  Vector3d transNoise(0.2, 0.2, 0.2);
  Vector3d rotNoise(DEG2RAD(2.), DEG2RAD(2.), DEG2RAD(2.));
  Vector2d measNoise(2, 2);
  double stepLen = 1.0;
  Vector3d area(100, 100, 50);
  
  size_t win_min = 2;
  size_t win_max = num_of_state;
  size_t monte = 50;
  
  vector<int> successCountor;//存放各个frame对应的总路程误差
  vector<double> successVector;
  successCountor.resize( win_max );
  successVector.resize( win_max );
  for( size_t iter=0; iter<monte; ++iter ) {
    //生成一组新的数据
    StateVector states;
    LedLayout ledLayout;
    int num_of_led = 8;
    generateStates( states, initVector, transNoise, rotNoise , num_of_state, stepLen, 0);
    generateLedLayout( num_of_led,  area , ledLayout);
    for ( size_t i=0; i<states.size(); ++i ) {
      int num_of_visible_led = generateCameraMeas( states[i],  ledLayout, measNoise, cam, 0);//mode=1: 全可视 mode=0：FOV内可视
    }
    //生成结束
    //每个数据都用不同的window测量一次
    for (size_t win_iter=win_min; win_iter < win_max; ++win_iter) {
      //Window Optimize
      vector<Vector3d> est_SWM_trans, est_SWM_rot;
      double time_window;
      SWM(states, ledLayout, 10,  win_iter, est_SWM_trans, est_SWM_rot, time_window);
     
      vector<double> temp_err_trans, temp_err_rot;
      computeError( states, est_SWM_trans, est_SWM_rot, temp_err_trans, temp_err_rot );// 计算一次monte试验中，各个路程上的误差大小，接下来将他们求和
      
      double ARMSE_trans=accumulate(temp_err_trans.begin(),temp_err_trans.end(),0.0) / temp_err_trans.size();
      double ARMSE_rot=accumulate(temp_err_rot.begin(),temp_err_rot.end(),0.0) / temp_err_trans.size();

      cout<<"当前monte："<<iter<<" 当前win_size："<<win_iter<<" 总路程平均误差："<<ARMSE_trans<<endl;
      if( ARMSE_trans < 10 ) {
	successCountor[win_iter]++;
      }
    }//end win_iter
  }//end monte
  
  for( size_t i=win_min; i<win_max; ++i ) {
    successVector[i] = (double)successCountor[i] / monte;
    fout<<i<<" "<<successVector[i]<<endl;
  }
  
  fout.flush();
  fout.close();
}


/*

void exp_window_vs_pnp_reguLayout_limFOV() {
  ofstream fout;
  char filenametemp[]="exp_window_vs_pnp_reguLayout_limFOV.txt";
  fout.open( filenametemp );  

  int num_of_state = 80;
  Vector6d initVector;
  initVector << 10, 50, 0, 0, 0, 0;
  Vector3d transNoise(0.2, 0.2, 0.2);
  Vector3d rotNoise(DEG2RAD(2.), DEG2RAD(2.), DEG2RAD(2.));
  Vector2d measNoise(2, 2);
  double stepLen = 1.0;
  Vector3d area(100, 100, 50);
  Vector2d gab(30, 30);
  
  LedLayout ledLayout;
//   generateLedLayout( ledLayout );//固定layout
  generateLedLayout( gab,  area,  ledLayout, 0);
  
  size_t win_min = 2;
  size_t win_max = num_of_state;
  size_t monte = 100;
  
  vector< vector<double> > err_win_trans;//存放各个frame对应的总路程误差
  vector< vector<double> > err_pnp_trans;//存放各个frame对应的总路程误差
  err_win_trans.resize( monte );
  err_pnp_trans.resize( monte );
  
  int flag = 0;
  
  for( size_t iter=0; iter<monte; ++iter ) {
    //生成一组新的数据
    StateVector states;
    generateStates( states, initVector, transNoise, rotNoise , num_of_state, stepLen, 0);
    for ( size_t i=0; i<states.size(); ++i ) {
      int num_of_visible_led = generateCameraMeas( states[i],  ledLayout, measNoise, cam, 0);//mode=1: 全可视 mode=0：FOV内可视
      cout<<"观测到："<<num_of_visible_led<<endl;
    }
    //生成结束
    vector<Vector3d> stateWindowVector;
    double time_window;
    windowOptimize(states, ledLayout, 10,  10, stateWindowVector, time_window);
    
    //PnP Method
    vector<Vector3d> statePnPVector;
    double timePnP;
    pnpMethod(states, ledLayout, statePnPVector, timePnP);

    err_win_trans[iter].resize( num_of_state );
    err_pnp_trans[iter].resize( num_of_state );
    calError_trans( states, stateWindowVector, err_win_trans[iter]);//返回各个step的xy,z误差
    calError_trans( states, statePnPVector, err_pnp_trans[iter]);//返回各个step的xy,z误差  
    
//     showData( states, stateWindowVector, ledLayout, area );
//     showData( states, statePnPVector, ledLayout, area );
    
    //只记录一次数据，目的只是为了ledLayout
    if(flag == 0)  {
      saveSampledData( states, ledLayout, area);
      flag = 1;
    }
    
  }//end monte
  
  vector<double>  rmse_win_trans;
  vector<double>  rmse_pnp_trans;
  rmse_win_trans.resize(num_of_state);
  rmse_pnp_trans.resize(num_of_state);
  for( size_t i=0; i<num_of_state; ++i ) {
    double temp_win_trans = 0;
    double temp_pnp_trans = 0;
    for ( size_t j=0; j<monte; ++j ) {
      temp_win_trans += err_win_trans[j][i];
      temp_pnp_trans += err_pnp_trans[j][i];
    }
    rmse_win_trans[i] = temp_win_trans / monte;
    rmse_pnp_trans[i] = temp_pnp_trans / monte;
    fout<<i<<" "<<rmse_win_trans[i]<<" "<<rmse_pnp_trans[i]<<endl;
  }	
  
  fout.flush();
  fout.close();
}

void exp_window_vs_pnp_random_layout()
{
  ofstream fout;
  char filenametemp[]="exp_window_vs_pnp_random_layout.txt";
  fout.open( filenametemp );  

  int num_of_state = 80;
  int num_of_led = 8;
  Vector6d initVector;
  initVector << 10, 50, 0, 0, 0, 0;
  Vector3d transNoise(0.2, 0.2, 0.2);
  Vector3d rotNoise(DEG2RAD(2.), DEG2RAD(2.), DEG2RAD(2.));
  Vector2d measNoise(2, 2);
  double stepLen = 1.0;
  Vector3d area(100, 100, 50);
  
  size_t win_min = 2;
  size_t win_max = num_of_state;
  size_t monte = 1;
  
  vector<int> num_of_visible_led;
  num_of_visible_led.resize(num_of_state);
  
  vector< vector<double> > err_win_xy;//存放各个frame对应的总路程误差
  vector< vector<double> > err_win_z;
  vector< vector<double> > err_pnp_xy;//存放各个frame对应的总路程误差
  vector< vector<double> > err_pnp_z;
  err_win_xy.resize( monte );
  err_win_z.resize( monte );
  err_pnp_xy.resize( monte );
  err_pnp_z.resize( monte );
  
  for( size_t iter=0; iter<monte; ++iter ) {
    //生成一组新的数据
    StateVector states;
    LedLayout ledLayout;
    generateStates( states, initVector, transNoise, rotNoise , num_of_state, stepLen, 0);
    generateLedLayout( num_of_led,  area , ledLayout);
    for ( size_t i=0; i<states.size(); ++i ) {
      num_of_visible_led[i] = generateCameraMeas( states[i],  ledLayout, measNoise, cam, 0);//mode=1: 全可视 mode=0：FOV内可视
      cout<<"state "<<i<<" 有 "<<num_of_visible_led[i]<<" 个led可视"<<endl;
    }
    //保存数据
    saveSampledData( states, ledLayout, area);
    //生成结束 
    vector<Vector3d> stateWindowVector;
    double time_window;
    windowOptimize(states, ledLayout, 10,  10, stateWindowVector, time_window);
    
    //PnP Method
    vector<Vector3d> statePnPVector;
    double timePnP;
    pnpMethod(states, ledLayout, statePnPVector, timePnP);//这里的PnP程序需要一点修改，如果可见点<4个或是误差大于5，都置0
    
    err_win_xy[iter].resize( num_of_state );
    err_win_z[iter].resize( num_of_state );
    err_pnp_xy[iter].resize( num_of_state );
    err_pnp_z[iter].resize( num_of_state );
    calError_xy( states, stateWindowVector, err_win_xy[iter], err_win_z[iter]);//返回各个step的xy,z误差
    calError_xy( states, statePnPVector, err_pnp_xy[iter], err_pnp_z[iter]);//返回各个step的xy,z误差  

    cout<<"到这里了吗"<<endl;
    vector<double>  rmse_win_xy;
    vector<double>  rmse_win_z;
    vector<double>  rmse_pnp_xy;
    vector<double>  rmse_pnp_z;
    rmse_win_xy.resize(num_of_state);
    rmse_win_z.resize(num_of_state);
    rmse_pnp_xy.resize(num_of_state);
    rmse_pnp_z.resize(num_of_state);
    for( size_t i=0; i<num_of_state; ++i ) {
      double temp_win_xy = 0;
      double temp_win_z = 0;
      double temp_pnp_xy = 0;
      double temp_pnp_z = 0;
      for ( size_t j=0; j<monte; ++j ) {
	temp_win_xy += err_win_xy[j][i];
	temp_win_z += err_win_z[j][i];
	temp_pnp_xy += err_pnp_xy[j][i];
	temp_pnp_z += err_pnp_z[j][i];
      }
      rmse_win_xy[i] = temp_win_xy / monte;
      rmse_win_z[i] = temp_win_z / monte;
      rmse_pnp_xy[i] = temp_pnp_xy / monte;
      rmse_pnp_z[i] = temp_pnp_z / monte;
      fout<<i<<" "<<num_of_visible_led[i]<<" "
			  <<stateWindowVector[i].transpose()<<" "
			  <<statePnPVector[i].transpose();
      fout<<endl;
	    
      cout<<"当前 "<<i<<" 的window_xy误差为："<<rmse_win_xy[i]<<" 以及z误差为 "<<rmse_win_z[i]<<endl;
    }	
    showData( states, stateWindowVector, ledLayout, area );
  }//end monte

  
  fout.flush();
  fout.close();
}

//窗口大小 vs 精度
void exp_window_size_vs_accuracy()
{
  ofstream fout;
  char filenametemp[]="exp_window_size_vs_accuracy.txt";
  fout.open( filenametemp );  

 
  int num_of_state = 50;
  Vector6d initVector;
  initVector << 25, 50, 0, 0, 0, 0;
  Vector3d transNoise(0.2, 0.2, 0.2);
  Vector3d rotNoise(DEG2RAD(2.), DEG2RAD(2.), DEG2RAD(2.));
  Vector2d measNoise(2, 2);
  double stepLen = 1.0;
  Vector3d area(100, 100, 50);
  
  LedLayout ledLayout;
  generateLedLayout( ledLayout );//固定layout
  
  size_t win_min = 2;
  size_t win_max = num_of_state;
  size_t monte = 10;
  
  vector< vector<double> > err;//存放各个frame对应的总路程误差
  err.resize( monte );
  for( size_t iter=0; iter<monte; ++iter ) {
    //生成一组新的数据
    StateVector states;
    generateStates( states, initVector, transNoise, rotNoise , num_of_state, stepLen, 0);
    for ( size_t i=0; i<states.size(); ++i ) {
      int num_of_visible_led = generateCameraMeas( states[i],  ledLayout, measNoise, cam, 1);//mode=1: 全可视 mode=0：FOV内可视
    }
    //生成结束
    //每个数据都用不同的window测量一次
    err[iter].resize( win_max );
    for (size_t win_iter=win_min; win_iter < win_max; ++win_iter) {
      vector<Vector3d> stateWindowVector;
      double time_window;
      windowOptimize(states, ledLayout, 10,  win_iter, stateWindowVector, time_window);
     
      //计算各个window size的总路程平均误差, 求sum( sqrt(x^2+y^2+z^2) )/N
      calError_xyz(states, stateWindowVector, err[iter][win_iter]);
      cout<<"当前monte："<<iter<<" 当前win_size："<<win_iter<<" calError_xyz输出："<<err[iter][win_iter]<<endl;
    }//end win_iter
  }//end monte

  cout<<"到这里了吗"<<endl;
  vector<double>  rmse;
  rmse.resize(win_max);
  for( size_t i=win_min; i<win_max; ++i ) {
    double temp = 0;
    for ( size_t j=0; j<monte; ++j ) {
      temp += err[j][i];
    }
    cout<<"temp是"<<temp<<endl;
    rmse[i] = temp / monte;
    fout<<i<<" "<<rmse[i]<<endl;
  }	
  
  fout.flush();
  fout.close();
}


//窗口大小 vs 鲁棒性（相同VO估计）
void exp_window_size_vs_robustness()
{
  ofstream fout;
  char filenametemp[]="exp_window_size_vs_robustness.txt";
  fout.open( filenametemp );  

  int num_of_state = 50;
  int num_of_led = 1;
  Vector6d initVector;
  initVector << 25, 50, 0, 0, 0, 0;
  Vector3d transNoise(0.2, 0.2, 0.2);
  Vector3d rotNoise(DEG2RAD(2.), DEG2RAD(2.), DEG2RAD(2.));
  Vector2d measNoise(2, 2);
  double stepLen = 1.0;
  Vector3d area(100, 100, 50);
  
  size_t win_min = 2;
  size_t win_max = num_of_state;
  size_t monte = 50;
  
  vector<int> successCountor;//存放各个frame对应的总路程误差
  vector<double> successVector;
  successCountor.resize( win_max );
  successVector.resize( win_max );
  for( size_t iter=0; iter<monte; ++iter ) {
    //生成一组新的数据
    StateVector states;
    LedLayout ledLayout;
    generateStates( states, initVector, transNoise, rotNoise , num_of_state, stepLen, 0);
    generateLedLayout( num_of_led,  area , ledLayout);
    for ( size_t i=0; i<states.size(); ++i ) {
      int num_of_visible_led = generateCameraMeas( states[i],  ledLayout, measNoise, cam, 0);//mode=1: 全可视 mode=0：FOV内可视
    }
    //生成结束
    //每个数据都用不同的window测量一次
    for (size_t win_iter=win_min; win_iter < win_max; ++win_iter) {
      vector<Vector3d> stateWindowVector;
      double time_window;
      windowOptimize(states, ledLayout, 10,  win_iter, stateWindowVector, time_window);
     
      //判断当前是否算是成功估计：小于5算成功，成功的话对应window_size +1，
      double err = 0;
      calError_xyz( states, stateWindowVector, err);
      cout<<"当前monte："<<iter<<" 当前win_size："<<win_iter<<" calError_xyz输出："<<err<<endl;
      if( err < 5 ) {
	successCountor[win_iter]++;
      }
    }//end win_iter
  }//end monte
  
  for( size_t i=win_min; i<win_max; ++i ) {
    successVector[i] = (double)successCountor[i] / monte;
    fout<<i<<" "<<successVector[i]<<endl;
  }
  
  fout.flush();
  fout.close();
}

//路程大小 vs 精度（相同VO估计）
void exp_route_length_vs_accuracy()
{
  //读取相同VO估计
  StateVector states;
  LedLayout ledLayout;
  Vector3d area;
  loadSampledData( states, ledLayout, area );
  
  ofstream fout;
  char filenametemp[]="exp_route_length_vs_accuracy.txt";
  fout.open( filenametemp );  
  
  size_t win_min = 2;
  size_t win_max = states.size();
    
  for (size_t win_iter=win_min; win_iter < win_max; ++win_iter) {
    size_t monte = 1;
    for ( size_t iter=0; iter<monte; ++iter) {
      vector<Vector3d> stateWindowVector;
      double time_window;
      windowOptimize(states, ledLayout, 10,  win_iter, stateWindowVector, time_window);
      
      //统计RMSE
      vector< vector<double> > err_init(3), err_window(3);
      vector<double>  mean_window(3), var_window(3);
      for (size_t i=0; i<3; ++i) {
	err_window[i].resize( states.size() );
	err_init[i].resize( states.size() );
	for (size_t j=0; j<states.size(); ++j) {
	  //为error赋值
	  err_window[i][j] = states[j].pose.translation()[i] - stateWindowVector[j][i];
	}
      } 
      fout<<win_iter<<" ";
      //fout<<iter<<" ";
      for (size_t i=0; i<3; ++i) {	
	RMS(err_window[i], mean_window[i], var_window[i]);
	fout<<var_window[i]<<" ";
      }
      fout<<time_window<<" ";
      fout<<endl;
    }//完成一组测量
  }//完成全部测量
  
  fout.flush();
  fout.close();
}

*/

void calError_xyz( StateVector states, vector<Vector3d> stateVector,  double& err ) {
  for( size_t i=0; i<states.size(); ++i ) {
    Vector3d err_xyz( states[i].pose.translation()[0] - stateVector[i][0], 
				 states[i].pose.translation()[1] - stateVector[i][1],
				 states[i].pose.translation()[2] - stateVector[i][2] );
    err = sqrt( (err_xyz[0]*err_xyz[0] + err_xyz[1]*err_xyz[1] + err_xyz[2]*err_xyz[2]) );
  }
}

void calError_trans( StateVector states, vector<Vector3d> stateVector,  vector<double>& err_trans ) {
  for( size_t i=0; i<states.size(); ++i ) {
    Vector3d err_xyz( states[i].pose.translation()[0] - stateVector[i][0], 
				 states[i].pose.translation()[1] - stateVector[i][1],
				 states[i].pose.translation()[2] - stateVector[i][2] );
    err_trans[i] = sqrt( (err_xyz[0]*err_xyz[0] + err_xyz[1]*err_xyz[1] + err_xyz[2]*err_xyz[2]) );
  }
}

void calError_xy( StateVector states, vector<Vector3d> stateVector,  vector<double>& err_xy,  vector<double>& err_z )
{
  for( size_t i=0; i<states.size(); ++i ) {
    Vector2d xy( states[i].pose.translation()[0] - stateVector[i][0], 
			  states[i].pose.translation()[1] - stateVector[i][1]);
    err_xy[i] = sqrt( (xy[0]*xy[0] + xy[1]*xy[1]) );
    err_z[i] = abs( states[i].pose.translation()[2] - stateVector[i][2] );
  }	
}

void FPM(StateVector states, LedLayout ledLayout, vector<Vector3d>& trans, vector<Vector3d>& rot, double& time) {
    chrono::steady_clock::time_point t1_PnP = chrono::steady_clock::now();
    
    //初始值默认是有的, 崩溃原因：初始值没算出来，后面却还想要赋值
//     Vector3d init_t(states[0].init_pose.translation()[0], states[0].init_pose.translation()[1], states[0].init_pose.translation()[2]);
//     statePnPVector.push_back(init_t);
    for (size_t i=0; i<states.size(); ++i) {
      vector<cv::Point3f> pts_3d;
      vector<cv::Point2f> pts_2d;
      for( size_t j=0; j<states[i].cam_meas.size(); ++j) {
	cv::Point3f xyz(ledLayout[states[i].cam_meas[j].fromLed].xyz[0], 
					ledLayout[states[i].cam_meas[j].fromLed].xyz[1], 
					ledLayout[states[i].cam_meas[j].fromLed].xyz[2]);
	cv::Point2f uv(states[i].cam_meas[j].uv[0], states[i].cam_meas[j].uv[1]);
	pts_3d.push_back( xyz );
	pts_2d.push_back( uv );
      }
      if (pts_3d.size() <= 3) {
	//置0
	Vector3d t0(0, 0, 0);
	Vector3d r0(0, 0, 0);
	trans.push_back(t0);
	trans.push_back(r0);
	//statePnPVector.push_back(statePnPVector.back());//默认为上一次的估计值
      } else {
	Mat r0, t0;
	cv::solvePnP ( pts_3d, pts_2d, K, Mat(), r0, t0, false, CV_EPNP);  //CV_EPNP适合用于大于四个点
	Vector3d t(-t0.at<double>(0,0), -t0.at<double>(0,1), -t0.at<double>(0,2));
	Vector3d r(RAD2DEG(-r0.at<double>(0,0)), RAD2DEG(-r0.at<double>(0,1)), RAD2DEG(-r0.at<double>(0,2)));
	double err = sqrt( pow(states[i].pose.translation()[0] - t[0], 2) + pow(states[i].pose.translation()[1] - t[1], 2) + pow(states[i].pose.translation()[2] - t[2], 2) );
	if( err < 5 ) {
	  trans.push_back(t);
	  rot.push_back(r);
	} else {
	  trans.push_back(t);
	  rot.push_back(r);
	  
// 	  测试success rate时才用
// 	  Vector3d t0(0, 0, 0);
// 	  statePnPVector.push_back(t0);
	}
// 	cout<<"PNP estimated Pose: "<<statePnPVector.back().transpose()<<endl;
      }
      chrono::steady_clock::time_point t2_PnP = chrono::steady_clock::now();
      chrono::duration<double> time_used_PnP = chrono::duration_cast<chrono::duration<double>>( t2_PnP-t1_PnP );
      time = time_used_PnP.count();
    }
    
    chrono::steady_clock::time_point t2_PnP = chrono::steady_clock::now();
    chrono::duration<double> time_used_PnP = chrono::duration_cast<chrono::duration<double>>( t2_PnP-t1_PnP );
    time = time_used_PnP.count();
}

void SWM(StateVector states, LedLayout ledLayout, int num_of_iter, int win_size, vector<Vector3d>& trans, vector<Vector3d>& rot, double& time) {
  chrono::steady_clock::time_point t1_window = chrono::steady_clock::now();
  
  StateVector stateVector;
  StateList stateList;//装window中待优化的状态

  for ( size_t i=0; i<win_size; ++i ) {
    stateList.push_back(states[i]);
  }
  //第一次优化
  optimizeBA( stateList, ledLayout, num_of_iter);//优化后，把优化的值重复赋值到winFrame里去，为下次优化做准备

  for (size_t state_pointer=win_size; state_pointer<states.size(); ++state_pointer) {
    //用上一次最后一个pose的值乘以里程计测量值得到更好的初始值
    states[state_pointer].init_pose = states[state_pointer].odometry * stateList.back().init_pose;
    stateList.push_back(states[state_pointer]);
    stateVector.push_back(stateList.front());
    stateList.pop_front();
    optimizeBA( stateList, ledLayout, num_of_iter);//优化后，把优化的值重复赋值到winFrame里去，为下次优化做准备
  }
  //最后一次优化好后，再把剩下的29个状态装进stateVector里
  list<State>::iterator itList;
  for( itList = stateList.begin(); itList != stateList.end(); ++itList ) {
    stateVector.push_back(*itList);
  }
  
  for (size_t i=0; i<states.size(); ++i) {
    Vector3d rotVector( RAD2DEG(stateVector[i].init_pose.log()[3]), RAD2DEG(stateVector[i].init_pose.log()[4]),  RAD2DEG(stateVector[i].init_pose.log()[5]) );
    trans.push_back(stateVector[i].init_pose.translation());
    rot.push_back(rotVector);
  }
    
  //window优化结束
  chrono::steady_clock::time_point t2_window = chrono::steady_clock::now();
  chrono::duration<double> time_used_window = chrono::duration_cast<chrono::duration<double>>( t2_window-t1_window );
  time = time_used_window.count();
}

void optimizeBA ( StateList& stateList, LedLayout ledLayout, int num_of_iter){
 
  typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1,-1> > Block; 
  std::unique_ptr<Block::LinearSolverType> linearSolver ( new g2o::LinearSolverCSparse<Block::PoseMatrixType>());
  std::unique_ptr<Block> solver_ptr ( new Block ( std::move(linearSolver)));   
  g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg ( std::move(solver_ptr));
  g2o::SparseOptimizer optimizer;
  optimizer.setAlgorithm ( solver );
    
  list<State>::iterator itList;
  
  //Adding vertexs
  vector<VertexSE3LieAlgebra*> vectices;//T_wc情况
  cerr << "Optimization: Adding robot poses ... ";

  for( itList = stateList.begin(); itList != stateList.end(); ++itList ) {
    VertexSE3LieAlgebra* v = new VertexSE3LieAlgebra;
    v->setId((*itList).id);
    v->setEstimate((*itList).init_pose);
    vectices.push_back(v);
    optimizer.addVertex(v);
  }
  cerr << "done." << endl;
  
  //Adding odometry edges
  cerr << " Adding odometry measurements ...";
  itList = stateList.begin();
  itList++;
  for( itList; itList != stateList.end(); ++itList ) {
    EdgeSE3LieAlgebra* odometry = new EdgeSE3LieAlgebra;
    odometry->vertices()[0] = optimizer.vertex((*itList).id-1);
    odometry->vertices()[1] = optimizer.vertex((*itList).id);
    odometry->setMeasurement((*itList).odometry);
    odometry->setInformation((*itList).odom_information);
    optimizer.addEdge(odometry);
  }
  cerr << "done." << endl;

  int index = 0;
  for ( itList = stateList.begin(); itList != stateList.end(); ++itList ) {
    for ( size_t j=0; j<(*itList).cam_meas.size(); ++j ) {
      EdgeSE3ProjectXYZ2UVPoseOnly* edge  = new  EdgeSE3ProjectXYZ2UVPoseOnly();
      edge->setId(index);
      edge->vertices()[0] = optimizer.vertex( (*itList).id );
      edge->camera_ = &cam;
      edge->point_ = ledLayout[(*itList).cam_meas[j].fromLed].xyz;
      edge->setMeasurement( (*itList).cam_meas[j].uv );
      edge->setInformation( (*itList).cam_meas[j].information );
      index++;
      optimizer.addEdge( edge );
    }
  }
   
  //Start optimizing**********************************************************************复制部分
  
  optimizer.setVerbose(false);
  optimizer.initializeOptimization();
  cout<<"calling optimizing ..."<<endl;
  optimizer.optimize(num_of_iter);
  cout<<"end !"<<endl;
  
  //优化出来的结果还要更新到stateList上去
  int i = 0;
  for(itList = stateList.begin(); itList != stateList.end(); ++itList ) {
    (*itList).init_pose = vectices[i]->estimate();
    //cout<<"estimated Pose: "<<i<<" : "<<vectices[i]->estimate().translation().transpose()<<endl;
    i++;
  }
}
      
/*
// 实验内容： window大小对定位鲁棒性影响
// (1)随机稀疏LED布局400个；
// (2)window大小2-100；
// (3)测试50组实验定位成功率，定位成功的定义为x,y,z 的RMSE各小于20
// (4)迭代次数设为30；
// (5)需要相同的初始位置估计（从VO中来）；
void exp_window_size_vs_robustness() {
ofstream fout;
  char filenametemp[]="exp_window_size_vs_robustness_output.txt";
  fout.open( filenametemp );  
  
  size_t win_min = 2;
  size_t win_max = 20;
  numOfTotalFrame  = 500; //做参数比较实验，窗口都用100的
  numOfLED = 300;
  double max_RMSE = 20;
  
  // 参数比较试验中，初始值都设置成相同，以排除不同初始值对结果的影响
  Simulator simulator;
  simulator.simulate(numOfTotalFrame, numOfLED);
  
  //记录下改组数据
  
  for (int win_iter=win_min; win_iter < win_max; ++win_iter) {
    size_t monte = 30;//每个数量led 迭代2次
    size_t success_time = 0;
    double success_rate = 0;
    for ( size_t iter=0; iter<monte; ++iter) {
      //Window Optimize
      vector<Vector3d> stateWindowVector;
      double timeWindow;
      windowOptimize(simulator, 5,  win_iter, stateWindowVector, timeWindow);//由于win_iter小造成卡顿
      
      //统计RMSE
      vector< vector<double> > err_init(3), err_batch(3), err_window(3), err_PnP(3);
      vector<double>  mean_window(3), var_window(3);
      for (size_t i=0; i<3; ++i) {
	err_window[i].resize(numOfTotalFrame);
	for (size_t j=0; j<numOfTotalFrame; ++j) {
	  //为error赋值
	  err_window[i][j] = (double)simulator.poses()[j].truePose.translation()[i] - (double)stateWindowVector[j].transpose()[i];
	}
      }
     
      fout<<"窗口大小： "<<win_iter<<" ";
      fout<<"仿真次数： "<<iter<<" ";
      for (size_t i=0; i<3; ++i) {	
	RMS(err_window[i], mean_window[i], var_window[i]);
	
	//次序为window|batch|pnp|vo
	//记录格式为  迭代次数 0  var_window(x) var_batch(x) var_PnP(x) var_init(x) 1 var_window(y) var_batch(y) var_PnP(y) var_init(y) 2 var_window(z) var_batch(z) var_PnP(z) var_init(z) 
	fout<<"&"<<" "
	      <<var_window[i]<<" ";
      }
      
      if ( var_window[0] < max_RMSE && var_window[1] < max_RMSE  && var_window[2] < max_RMSE ) {
	//估计成功，记录1
	success_time ++;
      }
      fout<<"消耗时间： "<<timeWindow<<" ";
      fout<<endl;
    }
    //完成一组测
    success_rate = (double)success_time / monte;
    fout<<"&&&&&&&&&&&&&&&& 该window大小的成功率： "<<success_rate<<"&&&&&&&&&&&&&&&& "<<endl;
    //cout<<"当前window大小： "<<win_iter<<endl;
  }
  fout.flush();
  fout.close();
  printf("文件输出完毕\r\n");
}

void exp_windoe_size_vs_accuracy() {
  ofstream fout;
  char filenametemp[]="exp_windoe_size_vs_accuracy_1led_.txt";
  fout.open( filenametemp );  
  
  numOfTotalFrame  = 50; //做参数比较实验，窗口都用50的
  size_t win_min = 2;
  size_t win_max = numOfTotalFrame;
  
  Simulator simulator;
  simulator.simulate(numOfTotalFrame, numOfLED);
  
  //读取原始数据；真实pose; 初始pose; vo测量值
  // read the raw data: true pose, initial pose, vo measurements.
  StateVector stateTrue;
  StateVector stateInit;
  StateVector stateVO;
  ifstream fin;  
 
  fin.open("true_pose.txt",ios::in);//ios::in 表示以只读的方式读取文件  
  if(fin.fail()) {  
    cout << "文件不存在." << endl;  
    fin.close(); 
    return;
  } 
  while (!fin.eof()) {
    double trans1, trans2, trans3, rot1, rot2, rot3;
    fin>>trans1>>trans2>>trans3>>rot1>>rot2>>rot3;
    Sophus::Vector6d se3;
    se3<< trans1, trans2, trans3, rot1, rot2, rot3;
    State state;
    state.pose = SE3::exp(se3);
    stateTrue.push_back(state);
    //cout<<"true poseeee"<<stateTrue.back().pose.translation().transpose()<<endl;
  }
  fin.close();  
  fin.open("initial_pose.txt",ios::in);//ios::in 表示以只读的方式读取文件  
  if(fin.fail()) {  
    cout << "文件不存在." << endl;  
    fin.close(); 
    return;
  } 
  while (!fin.eof()) {
    double trans1, trans2, trans3, rot1, rot2, rot3;
    fin>>trans1>>trans2>>trans3>>rot1>>rot2>>rot3;
    Sophus::Vector6d se3;
    se3<< trans1, trans2, trans3, rot1, rot2, rot3;
    State state;
    state.pose = SE3::exp(se3);
    stateInit.push_back(state);
    //cout<<"initial poseeee"<<stateInit.back().pose.translation().transpose()<<endl;
  }
  fin.close();  
  fin.open("vo_measurement.txt",ios::in);//ios::in 表示以只读的方式读取文件  
  if(fin.fail()) {  
    cout << "文件不存在." << endl;  
    fin.close(); 
    return;
  } 
  while (!fin.eof()) {
    double trans1, trans2, trans3, rot1, rot2, rot3;
    fin>>trans1>>trans2>>trans3>>rot1>>rot2>>rot3;
    Sophus::Vector6d se3;
    se3<< trans1, trans2, trans3, rot1, rot2, rot3;
    State state;
    state.odometry = SE3::exp(se3);
    stateVO.push_back(state);
    //cout<<"odometry"<<stateVO.back().pose.translation().transpose()<<endl;
  }
  fin.close();  
  
  // generate the led layout.
  
  // generate the led observation with respect to the true pose.
  

  for (int win_iter=win_min; win_iter < win_max; ++win_iter) {
    size_t monte = 1;
    for ( size_t iter=0; iter<monte; ++iter) {
      //不同的window跑出来的精度
      vector<Vector3d> stateWindowVector;
      double timeWindow;
      windowOptimize(simulator, 5,  win_iter, stateWindowVector, timeWindow);//由于win_iter小造成卡顿
      
      //统计RMSE
      vector< vector<double> > err_init(3), err_batch(3), err_window(3), err_PnP(3);
      vector<double>  mean_window(3), var_window(3);
      for (size_t i=0; i<3; ++i) {
	err_window[i].resize(numOfTotalFrame);
	for (size_t j=0; j<numOfTotalFrame; ++j) {
	  //为error赋值
	  err_window[i][j] = (double)simulator.poses()[j].truePose.translation()[i] - (double)stateWindowVector[j].transpose()[i];
	}
      }
     
      fout<<win_iter<<" ";
      //fout<<iter<<" ";
      for (size_t i=0; i<3; ++i) {	
	RMS(err_window[i], mean_window[i], var_window[i]);
	
	//次序为window|batch|pnp|vo
	//记录格式为  迭代次数 0  var_window(x) var_batch(x) var_PnP(x) var_init(x) 1 var_window(y) var_batch(y) var_PnP(y) var_init(y) 2 var_window(z) var_batch(z) var_PnP(z) var_init(z) 
	fout<<var_window[i]<<" ";
      }
      
      fout<<timeWindow<<" ";
      fout<<endl;
    }
    //完成一组测量
  }
  fout.flush();
  fout.close();
  printf("文件输出完毕\r\n");
}

void exp4_var() {
  ofstream fout;
  char filenametemp[]="exp4_output.txt";
  fout.open( filenametemp );  
  
  size_t led_min = 570;
  size_t led_max = 600;
  for (size_t led_iter=led_min; led_iter < led_max; ++led_iter) {
    int monte = 30;//每个数量led 迭代2次
    for ( size_t iter=0; iter<monte; ++iter) {
      //仿真数据
      Simulator simulator;
      simulator.simulate(numOfTotalFrame, led_iter);//生成不同数量的led
      //分别求出三种方法的结果
      
      //Window Optimize
      vector<Vector3d> stateWindowVector;
      double timeWindow;
      windowOptimize(simulator, 10,  10, stateWindowVector, timeWindow);
      
      //Batch Optimize
      vector<Vector3d> stateBatchVector;
      double timeBatch;
      batchOptimize(simulator, 30, stateBatchVector, timeBatch);
      
      //PnP Method
      vector<Vector3d> statePnPVector;
      double timePnP;
      pnpMethod(simulator, statePnPVector, timePnP);
      
      //统计RMSE
      vector< vector<double> > err_init(3), err_batch(3), err_window(3), err_PnP(3);
      vector<double> mean_PnP(3), var_PnP(3);
      vector<double> mean_batch(3), var_batch(3);
      vector<double>  mean_window(3), var_window(3);
      vector<double> mean_init(3), var_init(3);
      for (size_t i=0; i<3; ++i) {
	err_batch[i].resize(numOfTotalFrame);
	err_window[i].resize(numOfTotalFrame);
	err_PnP[i].resize(numOfTotalFrame);
	err_init[i].resize(numOfTotalFrame);
	for (size_t j=0; j<numOfTotalFrame; ++j) {
	  //为error赋值
	  err_batch[i][j] = (double)simulator.poses()[j].truePose.translation()[i] - (double)stateBatchVector[j].transpose()[i];
	  err_window[i][j] = (double)simulator.poses()[j].truePose.translation()[i] - (double)stateWindowVector[j].transpose()[i];
	  err_PnP[i][j] = (double)simulator.poses()[j].truePose.translation()[i] - (double)statePnPVector[j].transpose()[i];
	  err_init[i][j] = (double)simulator.poses()[j].truePose.translation()[i] -(double)simulator.init_poses()[j].initPose.translation()[i];
	}
      }
      fout<<led_iter<<" ";
      fout<<iter<<" ";
      for (size_t i=0; i<3; ++i) {
	RMS(err_PnP[i], mean_PnP[i], var_PnP[i]);
	RMS(err_batch[i], mean_batch[i], var_batch[i]);
	RMS(err_window[i], mean_window[i], var_window[i]);
	RMS(err_init[i], mean_init[i], var_init[i]);
	
	//次序为window|batch|pnp|vo
	//记录格式为  迭代次数 0  var_window(x) var_batch(x) var_PnP(x) var_init(x) 1 var_window(y) var_batch(y) var_PnP(y) var_init(y) 2 var_window(z) var_batch(z) var_PnP(z) var_init(z) 
	fout<<i<<" "
	      <<var_window[i]<<" "
	      <<var_batch[i]<<" "
	      <<var_PnP[i]<<" "
	      <<var_init[i]<<" ";
      }
      fout<<endl;
    }
  }
  fout.flush();
  fout.close();
  printf("文件输出完毕\r\n");
}
//在led规整分布下，不同方法行走直线的估计
void exp1_var() {
  ofstream fout;
  char filenametemp[]="exp1_var.txt";
  fout.open( filenametemp );  
  
  numOfTotalFrame = 200;
  int monte = 1;
  for ( size_t iter=0; iter<monte; ++iter) {
    //仿真数据
    Simulator simulator;
    simulator.simulate(numOfTotalFrame, numOfLED);
    //分别求出三种方法的结果
    
     //The initial data should be the same.
    //记录下改组数据
    ofstream fin;
    
    char filenametemp1[] ="true_pose.txt";
    fin.open( filenametemp1 );  
    //真实坐标
    for (int i=0; i<simulator.poses().size(); ++i) {
      SE3 truePose =simulator.poses()[i].truePose;
      fin<<truePose.log().transpose()<<endl;
    }
    fin.flush();
    fin.close();
    
    char filenametemp2[] ="initial_pose.txt";
    fin.open( filenametemp2 );  
    for (int i=0; i<simulator.poses().size(); ++i) {
      SE3 initPose = simulator.init_poses()[i].initPose;
      fin<<initPose.log().transpose()<<endl;
    }
    fin.flush();
    fin.close();
    
    char filenametemp3[] ="vo_measurement.txt";
    fin.open( filenametemp3 );  
    for (int i=0; i<simulator.odometry().size(); ++i) {
      SE3 odometry = simulator.odometry()[i].noiseTransf;
      fin<<odometry.log().transpose()<<endl;
    }
    
    fin.flush();
    fin.close();
    
    //Window Optimize
    vector<Vector3d> stateWindowVector;
    double timeWindow;
    windowOptimize(simulator, 10,  10, stateWindowVector, timeWindow);
    
    //Batch Optimize
    vector<Vector3d> stateBatchVector;
    double timeBatch;
    batchOptimize(simulator, 30, stateBatchVector, timeBatch);
    
    //PnP Method
    vector<Vector3d> statePnPVector;
    double timePnP;
    pnpMethod(simulator, statePnPVector, timePnP);
    
    //统计RMSE
    vector< vector<double> > err_init(3), err_batch(3), err_window(3), err_PnP(3);
    vector<double> mean_PnP(3), var_PnP(3);
    vector<double> mean_batch(3), var_batch(3);
    vector<double>  mean_window(3), var_window(3);
    vector<double> mean_init(3), var_init(3);
    for (size_t i=0; i<3; ++i) {
      err_batch[i].resize(numOfTotalFrame);
      err_window[i].resize(numOfTotalFrame);
      err_PnP[i].resize(numOfTotalFrame);
      err_init[i].resize(numOfTotalFrame);
      for (size_t j=0; j<numOfTotalFrame; ++j) {
	//为error赋值
	err_batch[i][j] = (double)simulator.poses()[j].truePose.translation()[i] - (double)stateBatchVector[j].transpose()[i];
	err_window[i][j] = (double)simulator.poses()[j].truePose.translation()[i] - (double)stateWindowVector[j].transpose()[i];
	err_PnP[i][j] = (double)simulator.poses()[j].truePose.translation()[i] - (double)statePnPVector[j].transpose()[i];
	err_init[i][j] = (double)simulator.poses()[j].truePose.translation()[i] -(double)simulator.init_poses()[j].initPose.translation()[i];
      }
    }
    fout<<iter<<" ";
    for (size_t i=0; i<3; ++i) {
      RMS(err_PnP[i], mean_PnP[i], var_PnP[i]);
      RMS(err_batch[i], mean_batch[i], var_batch[i]);
      RMS(err_window[i], mean_window[i], var_window[i]);
      RMS(err_init[i], mean_init[i], var_init[i]);
      
      //次序为window|batch|pnp|vo
      //记录格式为  迭代次数 0  var_window(x) var_batch(x) var_PnP(x) var_init(x) 1 var_window(y) var_batch(y) var_PnP(y) var_init(y) 2 var_window(z) var_batch(z) var_PnP(z) var_init(z) 
      fout<<i<<" "
	     <<var_window[i]<<" "
	     <<var_batch[i]<<" "
	     <<var_PnP[i]<<" "
	     <<var_init[i]<<" ";
    }
    fout<<endl;
    
    //显示数据
    cout<<"PnP estimated Pose: "<<statePnPVector.back().transpose()<<endl;
    cout<<"Batch estimated Pose: "<<stateBatchVector.back().transpose()<<endl;
    cout<<"Window estimated Pose: "<<stateWindowVector.back().transpose()<<endl;

    cout<<"PnP  method's solve time cost = "<<timePnP<<" seconds. "<<endl;
    cout<<"Batch method's solve time cost = "<<timeBatch<<" seconds. "<<endl;
    cout<<"Window  method's solve time cost = "<<timeWindow<<" seconds. "<<endl;

    cout<<"\nwindow error X: "<<var_window[0]<<"  vs  batch error X: "<<var_batch[0]<<endl;
    cout<<"window error Y: "<<var_window[1]<<"  vs  batch error Y: "<<var_batch[1]<<endl;
    cout<<"window error Z: "<<var_window[2]<<"  vs  batch error Z: "<<var_batch[2]<<"\n"<<endl;
    
    cout<<"  vs  PnP error X: "<<var_PnP[0]<<"  vs  VO error X: "<<var_init[0]<<endl;
    cout<<"  vs  PnP error Y: "<<var_PnP[1]<<"  vs  VO error Y: "<<var_init[1]<<endl;
    cout<<"  vs  PnP error Z: "<<var_PnP[2]<<"  vs  VO error Z: "<<var_init[2]<<endl;
  
    //作图
    Mat plot_image( 504, 1004, CV_8UC3 );  
    cv::namedWindow("Test"); 
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
      
    cv::rectangle( plot_image,  cv::Rect(2,2,1000,500), cv::Scalar( 0, 0, 0 ), 1, 8 );  
    
    for ( size_t i; i<simulator.landmarks().size(); ++i) {
      cv::Point p(2+simulator.landmarks()[i].xyzLED[0], 2+simulator.landmarks()[i].xyzLED[1]);
      cv::circle(plot_image, p, 2, cv::Scalar(0, 0, 255));
    }
    
    for ( size_t i; i<simulator.landmarkObservations().size(); ++i ) {
      cv::Point p(2+simulator.landmarks()[ simulator.landmarkObservations()[i].from ].xyzLED[0], 2+simulator.landmarks()[ simulator.landmarkObservations()[i].from ].xyzLED[1]);
      cv::circle(plot_image, p, 3, cv::Scalar(0, 0, 255), -1);
    }
    
    //True
    cv::Point pre_point, cur_point;
    pre_point = cv::Point(simulator.poses()[0].truePose.translation()[0], simulator.poses()[0].truePose.translation()[1]);
    for (size_t i=1; i<simulator.poses().size(); ++i) {
      cur_point = cv::Point(simulator.poses()[i].truePose.translation()[0], simulator.poses()[i].truePose.translation()[1]);
      cv::line(plot_image, pre_point, cur_point, cv::Scalar(0, 0, 0));  
      pre_point = cur_point;
    }
    //VO
    pre_point = cv::Point(simulator.init_poses()[0].initPose.translation()[0], simulator.init_poses()[0].initPose.translation()[1]);
    for (size_t i=1; i<simulator.poses().size(); ++i) {
      cur_point = cv::Point(simulator.init_poses()[i].initPose.translation()[0], simulator.init_poses()[i].initPose.translation()[1]);
      cv::line(plot_image, pre_point, cur_point, cv::Scalar(0, 255, 0));  
      pre_point = cur_point;
    }
    
    //Window
    pre_point = cv::Point(stateWindowVector[0].transpose()[0], stateWindowVector[0].transpose()[1]);
    for (size_t i=1; i<simulator.poses().size(); ++i) {
      cur_point = cv::Point(stateWindowVector[i].transpose()[0], stateWindowVector[i].transpose()[1]);
      cv::line(plot_image, pre_point, cur_point, cv::Scalar(0, 0, 255));  
      pre_point = cur_point;
    }
    //Batch
    pre_point = cv::Point(stateBatchVector[0].transpose()[0], stateBatchVector[0].transpose()[1]);
    for (size_t i=1; i<simulator.poses().size(); ++i) {
      cur_point = cv::Point(stateBatchVector[i].transpose()[0], stateBatchVector[i].transpose()[1]);
      cv::line(plot_image, pre_point, cur_point, cv::Scalar(255, 0, 255));  
      pre_point = cur_point;
    }
    //PnP
    pre_point = cv::Point(statePnPVector[0].transpose()[0], statePnPVector[0].transpose()[1]);
    for (size_t i=1; i<simulator.poses().size(); ++i) {
      cur_point = cv::Point(statePnPVector[i].transpose()[0], statePnPVector[i].transpose()[1]);
      cv::line(plot_image, pre_point, cur_point, cv::Scalar(227,207,87));  
      pre_point = cur_point;
    }
    
    cv::imshow( "Test", plot_image );  
    cv::waitKey( 0 );  
    
  }
  
  fout.flush();
  fout.close();
  printf("文件输出完毕\r\n");
}

void pnpMethod(Simulator& simulator, vector<Vector3d>& statePnPVector, double& time) {
    chrono::steady_clock::time_point t1_PnP = chrono::steady_clock::now();
    
    const Simulator::LandmarkEdgeVector& landmakrObsv = simulator.landmarkObservations();
    //初始值默认是有的, 崩溃原因：初始值没算出来，后面却还想要赋值
    Vector3d init_t(simulator.init_poses()[0].initPose.translation()[0], simulator.init_poses()[0].initPose.translation()[1], simulator.init_poses()[0].initPose.translation()[2]);
    statePnPVector.push_back(init_t);
    for (size_t i=1; i<simulator.poses().size(); ++i) {
      //需要在landmakrObsv中找属于状态i的观测值
      vector<cv::Point3f> pts_3d;
      vector<cv::Point2f> pts_2d;
      for (size_t j=0; j<landmakrObsv.size(); ++j) {
	if (landmakrObsv[j].to == i) { 
	    cv::Point3f xyzLED(simulator.landmarks()[landmakrObsv[j].from].xyzLED[0], 
					    simulator.landmarks()[landmakrObsv[j].from].xyzLED[1], 
					    simulator.landmarks()[landmakrObsv[j].from].xyzLED[2]);
	    cv::Point2f xyPro(landmakrObsv[j].simulateMeas[0], landmakrObsv[j].simulateMeas[1]);
	    pts_3d.push_back(xyzLED);
	    pts_2d.push_back(xyPro);
	}
      }    
      if (pts_3d.size() <= 3) {
	//在图中标一个红叉
	statePnPVector.push_back(statePnPVector.back());//默认为上一次的估计值
      } else {
	Mat r0, t0;
	cv::solvePnP ( pts_3d, pts_2d, K, Mat(), r0, t0, false, CV_EPNP);  //CV_EPNP适合用于大于四个点
	Vector3d t(-t0.at<double>(0,0), -t0.at<double>(0,1), -t0.at<double>(0,2));
	statePnPVector.push_back(t);
	//cout<<"PNP estimated Pose: "<<statePnPVector.back().transpose()<<endl;
      }
    }
    
    chrono::steady_clock::time_point t2_PnP = chrono::steady_clock::now();
    chrono::duration<double> time_used_PnP = chrono::duration_cast<chrono::duration<double>>( t2_PnP-t1_PnP );
    time = time_used_PnP.count();
}*/

void RMS(vector<double>& vec, double& mean, double& variance) {
  
  double sum = std::accumulate(std::begin(vec), std::end(vec), 0.0);  
  mean =  sum / vec.size(); //均值  
  
  double accum  = 0.0;  
  std::for_each (std::begin(vec), std::end(vec), [&](const double d) {  
    accum  += (d-mean)*(d-mean);  
  });  
  
  variance = accum/vec.size();//方差
  double stdev = sqrt(accum/(vec.size()-1)); //标准差
}

// void batchOptimize(Simulator& simulator, int numIter, vector<Vector3d>& trans, double& time) {
//   chrono::steady_clock::time_point t1_batch = chrono::steady_clock::now();
//   
//   typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1,-1> > Block;  // pose 维度为 6, landmark 维度为 3, 误差的维度是2
//   Block::LinearSolverType* linearSolver = new g2o::LinearSolverCSparse<Block::PoseMatrixType>(); // 线性方程求解器
//   Block* solver_ptr = new Block ( linearSolver );     // 矩阵块求解器
//   g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg ( solver_ptr );//对初始值要求高，但是只有它能处理多LED情况
//   g2o::SparseOptimizer optimizer;
//   optimizer.setAlgorithm ( solver );
//  
//   //Adding vertexs
//   vector<VertexSE3LieAlgebra*> vectices;//T_wc情况
//   cerr << "Optimization: Adding robot poses ... ";
//   for (int i=0; i<simulator.init_poses().size(); ++i) {
//     const Simulator::InitPose& p =simulator.init_poses()[i];//初始pose
//     const Simulator::Pose& ptrue =simulator.poses()[i];//pose
//     VertexSE3LieAlgebra* v = new VertexSE3LieAlgebra;
//     v->setId(p.id);
//     v->setEstimate(p.initPose );
//     vectices.push_back(v);
//     optimizer.addVertex(v);
//   }
//   cerr << "done." << endl;
//   
//   //Adding odometry edges
//   cerr << " Adding odometry measurements ...";
//   for (int i=0; i<simulator.odometry().size(); ++i) {
//     const Simulator::PoseEdge& poseEdge =simulator.odometry()[i];
//     
//     EdgeSE3LieAlgebra* odometry = new EdgeSE3LieAlgebra;
//     odometry->vertices()[0] = optimizer.vertex(poseEdge.from);
//     odometry->vertices()[1] = optimizer.vertex(poseEdge.to);
//     odometry->setMeasurement(poseEdge.noiseTransf);//用noise就对了！
//     odometry->setInformation(poseEdge.information);
//     optimizer.addEdge(odometry);
//   }
//   cerr << "done." << endl;
// 
//   int index = 0;
//   for (int i=0; i<simulator.landmarkObservations().size(); ++i) 
//   {
//         EdgeSE3ProjectXYZ2UVPoseOnly* edge  = new  EdgeSE3ProjectXYZ2UVPoseOnly();
//         edge->setId(index);
// 	edge->vertices()[0] = optimizer.vertex(simulator.landmarkObservations()[i].to);
//         edge->camera_ = &cam;
//         edge->point_ = simulator.landmarks()[ simulator.landmarkObservations()[i].from].xyzLED;
//         edge->setMeasurement( simulator.landmarkObservations()[i].simulateMeas );
//         edge->setInformation( simulator.landmarkObservations()[i].information );
//         optimizer.addEdge( edge );
// 	index++;
//     }
//    
//   //Start optimizing**********************************************************************复制部分
//     optimizer.setVerbose(false);
//     optimizer.initializeOptimization();
//     cout<<"calling optimizing ..."<<endl;
//     optimizer.optimize(numIter);
//     cout<<"end !"<<endl;
//     chrono::steady_clock::time_point t2_batch = chrono::steady_clock::now();
//     chrono::duration<double> time_used_batch = chrono::duration_cast<chrono::duration<double>>( t2_batch-t1_batch );
//     time = time_used_batch.count();
//     
//     for (size_t i=0; i<simulator.poses().size(); ++i) {
//       trans.push_back(vectices[i]->estimate().translation());
//     }
// }

//   //真实坐标
//   for (size_t i=0; i<simulator.poses().size(); ++i) {
//     SE3 truePose =simulator.poses()[i].truePose;
//     cout<<"true  pose : "<<i<<" : "<<truePose.translation().transpose()<<endl;
//   }

//10月18日周三下午1：40记录
   // 输出全部优化值
//   for (int i=0; i<simulator.poses().size(); ++i) {
//     cout<<"批量estimated Pose: "<<i<<" : "<<vectices[i]->estimate().translation().transpose()<<endl;
//   }
//   for (int i=0; i<stateVector.size(); ++i) {
//     cout<<"窗口estimated Pose: "<<i<<" : "<<stateVector[i].pose.translation().transpose()<<endl;
//   }

//10月16日周一晚上十点半记录
//     //从前向后显示stateList中的数据   
//     cout<<"listOne.begin() --- listOne.end():"<<endl;
//     list<State>::iterator itList;
//     for( itList = stateList.begin(); itList != stateList.end(); ++itList ) cout<<(*itList).id<<" ";
//     cout<<endl;
/* 以及
 //T_wc Sliding Window版本**********************************************************************复制部分
  // 初始化g2o
  typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1,-1> > Block;  // pose 维度为 6, landmark 维度为 3, 误差的维度是2
  Block::LinearSolverType* linearSolver = new g2o::LinearSolverCSparse<Block::PoseMatrixType>(); // 线性方程求解器
  Block* solver_ptr = new Block ( linearSolver );     // 矩阵块求解器
  g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg ( solver_ptr );//对初始值要求高，但是只有它能处理多LED情况
  g2o::SparseOptimizer optimizer;
  optimizer.setAlgorithm ( solver );
 
  //Adding vertexs
  vector<VertexSE3LieAlgebra*> vectices;//T_wc情况
  cerr << "Optimization: Adding robot poses ... ";
  for (int i=0; i<simulator.init_poses().size(); ++i) {
    const Simulator::InitPose& p =simulator.init_poses()[i];//初始pose
    const Simulator::Pose& ptrue =simulator.poses()[i];//pose
    VertexSE3LieAlgebra* v = new VertexSE3LieAlgebra;
    v->setId(p.id);
    v->setEstimate(p.initPose );
    vectices.push_back(v);
    optimizer.addVertex(v);
  }
  cerr << "done." << endl;
  
  //Adding odometry edges
  cerr << " Adding odometry measurements ...";
  for (int i=0; i<simulator.odometry().size(); ++i) {
    const Simulator::PoseEdge& poseEdge =simulator.odometry()[i];
    
    EdgeSE3LieAlgebra* odometry = new EdgeSE3LieAlgebra;
    odometry->vertices()[0] = optimizer.vertex(poseEdge.from);
    odometry->vertices()[1] = optimizer.vertex(poseEdge.to);
    odometry->setMeasurement(poseEdge.noiseTransf);//用noise就对了！
    odometry->setInformation(poseEdge.information);
    optimizer.addEdge(odometry);
  }
  cerr << "done." << endl;

  int index = 0;
  for (int i=0; i<simulator.landmarkObservations().size(); ++i) 
  {
        EdgeSE3ProjectXYZ2UVPoseOnly* edge  = new  EdgeSE3ProjectXYZ2UVPoseOnly();
        edge->setId(index);
	edge->vertices()[0] = optimizer.vertex(simulator.landmarkObservations()[i].to);
        edge->camera_ = &cam;
        edge->point_ = simulator.landmarks()[ simulator.landmarkObservations()[i].from].xyzLED;
        edge->setMeasurement( simulator.landmarkObservations()[i].simulateMeas );
        edge->setInformation( simulator.landmarkObservations()[i].information );
        optimizer.addEdge( edge );
	index++;
    }
   
  //Start optimizing**********************************************************************复制部分
  
  optimizer.setVerbose(true);
  optimizer.initializeOptimization();
  cout<<"calling optimizing ..."<<endl;
  optimizer.optimize(10);
  cout<<"end !"<<endl;
  
  //真实坐标
  for (int i=0; i<simulator.poses().size(); ++i) {
    SE3 truePose =simulator.poses()[i].truePose;
    cout<<"true  pose : "<<i<<" : "<<truePose.translation().transpose()<<endl;
  }
  
  //初始值
  for (int i=0; i<simulator.init_poses().size(); ++i) {
    SE3 initPose =simulator.init_poses()[i].initPose;
    cout<<"init  pose : "<<i<<" : "<<initPose.translation().transpose()<<endl;
  }
  
   // 输出优化值
  for (int i=0; i<simulator.poses().size(); ++i) {
    cout<<"estimated Pose: "<<i<<" : "<<vectices[i]->estimate().translation().transpose()<<endl;
  }

  double init_error = 0;
  for (size_t i = 0; i<numOfTotalFrame; ++i) {
    double e = (double)simulator.poses()[i].truePose.translation()[0] - simulator.init_poses()[i].initPose.translation()[0];
    init_error += e*e;
  }
  cout<<"init_error: "<<init_error<<endl;
  
  double error = 0;
  for (size_t i = 0; i<numOfTotalFrame; ++i) {
    double e = (double)simulator.poses()[i].truePose.translation()[0] - (double)vectices[i]->estimate().translation()[0];
    error += e*e;
  }
   cout<<"error: "<<error<<endl;
   */
//T_wc 版本Start optimizing**********************************************************************复制部分
//   // 初始化g2o
//   typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1,-1> > Block;  // pose 维度为 6, landmark 维度为 3, 误差的维度是2
//   Block::LinearSolverType* linearSolver = new g2o::LinearSolverCSparse<Block::PoseMatrixType>(); // 线性方程求解器
//   Block* solver_ptr = new Block ( linearSolver );     // 矩阵块求解器
//   g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg ( solver_ptr );//对初始值要求高，但是只有它能处理多LED情况
//   //g2o::OptimizationAlgorithmGaussNewton* solver = new g2o::OptimizationAlgorithmGaussNewton( solver_ptr );//对初始值要求低，精度最高，但只能处理单LED情况
//   //g2o::OptimizationAlgorithmDogleg* solver = new g2o::OptimizationAlgorithmDogleg( solver_ptr );
//   g2o::SparseOptimizer optimizer;
//   optimizer.setAlgorithm ( solver );
//  
//   //Adding vertexs
//   vector<VertexSE3LieAlgebra*> vectices;//T_wc情况
//    //vector<g2o::VertexSE3Expmap*> vectices;//T_cw情况
//   
//   cerr << "Optimization: Adding robot poses ... ";
//   for (int i=0; i<simulator.init_poses().size(); ++i) {
//     const Simulator::InitPose& p =simulator.init_poses()[i];//初始pose
//     const Simulator::Pose& ptrue =simulator.poses()[i];//pose
//     
//      VertexSE3LieAlgebra* v = new VertexSE3LieAlgebra;
//     
//     //v->setId(ptrue.id);//初始点采用正确值
//     //v->setEstimate(ptrue.simulatePose);//假设simulatePose是通过单应矩阵H求得的
//     
//     v->setId(p.id);
//    v->setEstimate(p.initPose );
//    //v->setEstimate(simulator.init_poses()[0].initPose);//实验1：统一用第一个frame的坐标作初始值 实验2：初始值未知，只知道最近的灯的坐标
//     if ( i==0 ) {
//       //v->setFixed(true);
// //       v->setEstimate(noiseInit);
//     }
//     
//     
// //     g2o::VertexSE3Expmap* v = new g2o::VertexSE3Expmap;
// //     
// //     //v->setId(ptrue.id);//初始点采用正确值
// //     v->setEstimate(    g2o::SE3Quat (
// // //         ptrue.simulatePose.rotation_matrix(), 
// // //         ptrue.simulatePose.translation()
// //       
// // // 	simulator.init_poses()[0].initPose.rotation_matrix(), 
// // //         simulator.init_poses()[0].initPose.translation()
// //       
// //       p.initPose.rotation_matrix(), 
// //       p.initPose.translation()
// //       
// // 	
// //     ) );//假设simulatePose是通过单应矩阵H求得的
// // 
// //     v->setId(p.id);
// // //      if ( i==0 ) {
// // //       v->setFixed(true);
// // //     }
//     
//     vectices.push_back(v);
//     optimizer.addVertex(v);
//   }
//   cerr << "done." << endl;
//   
//   //Adding odometry edges
//   cerr << " Adding odometry measurements ...";
//   for (int i=0; i<simulator.odometry().size(); ++i) {
//     const Simulator::PoseEdge& poseEdge =simulator.odometry()[i];
//     
//     EdgeSE3LieAlgebra* odometry = new EdgeSE3LieAlgebra;
//     odometry->vertices()[0] = optimizer.vertex(poseEdge.from);
//     odometry->vertices()[1] = optimizer.vertex(poseEdge.to);
//     odometry->setMeasurement(poseEdge.noiseTransf);//用noise就对了！
//     //cout<<"带噪声测量值："<<poseEdge.noiseTransf.translation().transpose()<<endl;
//     odometry->setInformation(poseEdge.information);
//     optimizer.addEdge(odometry);
//   }
//   cerr << "done." << endl;
// 
//   int index = 0;
//   for (int i=0; i<simulator.landmarkObservations().size(); ++i) 
//   {
//         EdgeSE3ProjectXYZ2UVPoseOnly* edge  = new  EdgeSE3ProjectXYZ2UVPoseOnly();
//         edge->setId(index);
// 	edge->vertices()[0] = optimizer.vertex(simulator.landmarkObservations()[i].to);
//         edge->camera_ = &cam;
//         edge->point_ = simulator.landmarks()[ simulator.landmarkObservations()[i].from].xyzLED;
//         edge->setMeasurement( simulator.landmarkObservations()[i].simulateMeas );
//         edge->setInformation( simulator.landmarkObservations()[i].information );
// 	//edge->setInformation( Eigen::Matrix2d::Identity() );
// 	//cout<<"检验pose "<<simulator.landmarkObservations()[i].to<<", LED "<<simulator.landmarkObservations()[i].from<<" : "<<simulator.landmarkObservations()[i].simulateMeas.transpose() <<endl;
//         optimizer.addEdge( edge );
// 	index++;
//     }
//    
  //Start optimizing**********************************************************************复制部分
  
  
  //显示数据
  //   cout<<"PnP estimated Pose: "<<statePnPVector.back().transpose()<<endl;
//   cout<<"Batch estimated Pose: "<<stateBatchVector.back().transpose()<<endl;
//   cout<<"Window estimated Pose: "<<stateWindowVector.back().transpose()<<endl;
// 
//   cout<<"PnP  method's solve time cost = "<<timePnP<<" seconds. "<<endl;
//   cout<<"Batch method's solve time cost = "<<timeBatch<<" seconds. "<<endl;
//   cout<<"Window  method's solve time cost = "<<timeWindow<<" seconds. "<<endl;

//   cout<<"\nwindow error X: "<<var_window[0]<<"  vs  batch error X: "<<var_batch[0]<<endl;
//   cout<<"window error Y: "<<var_window[1]<<"  vs  batch error Y: "<<var_batch[1]<<endl;
//   cout<<"window error Z: "<<var_window[2]<<"  vs  batch error Z: "<<var_batch[2]<<"\n"<<endl;
//   
//   cout<<"  vs  PnP error X: "<<var_PnP[0]<<"  vs  VO error X: "<<var_init[0]<<endl;
//   cout<<"  vs  PnP error Y: "<<var_PnP[1]<<"  vs  VO error Y: "<<var_init[1]<<endl;
//   cout<<"  vs  PnP error Z: "<<var_PnP[2]<<"  vs  VO error Z: "<<var_init[2]<<endl;
  
  
//作图
//   Mat plot_image( 504, 1004, CV_8UC3 );  
//   cv::namedWindow("Test"); 
//   for (size_t y=0; y<plot_image.rows; y++)        //遍历每一行每一列并设置其像素值  
//   {  
//     for (size_t x=0; x<plot_image.cols; x++)  
//     {
//       unsigned char* row_ptr = plot_image.ptr<unsigned char> ( y );  // row_ptr是第y行的头指针
//       unsigned char* data_ptr = &row_ptr[ x*plot_image.channels() ]; // data_ptr 指向待访问的像素数据，这里的&是指c语言的取地址
//       for ( int c = 0; c != plot_image.channels(); c++ )
//       {
// 	  data_ptr[c] = 255; // data为I(x,y)第c个通道的值
//       }
//     }  
//   } 
//     
//   cv::rectangle( plot_image,  cv::Rect(2,2,1000,500), cv::Scalar( 0, 0, 0 ), 1, 8 );  
//   
//   for ( size_t i; i<simulator.landmarks().size(); ++i) {
//     cv::Point p(2+simulator.landmarks()[i].xyzLED[0], 2+simulator.landmarks()[i].xyzLED[1]);
//     cv::circle(plot_image, p, 2, cv::Scalar(0, 0, 255));
//   }
//   
//   for ( size_t i; i<simulator.landmarkObservations().size(); ++i ) {
//     cv::Point p(2+simulator.landmarks()[ simulator.landmarkObservations()[i].from ].xyzLED[0], 2+simulator.landmarks()[ simulator.landmarkObservations()[i].from ].xyzLED[1]);
//     cv::circle(plot_image, p, 3, cv::Scalar(0, 0, 255), -1);
//   }
//   
//   //True
//   cv::Point pre_point, cur_point;
//   pre_point = cv::Point(simulator.poses()[0].truePose.translation()[0], simulator.poses()[0].truePose.translation()[1]);
//   for (size_t i=1; i<simulator.poses().size(); ++i) {
//     cur_point = cv::Point(simulator.poses()[i].truePose.translation()[0], simulator.poses()[i].truePose.translation()[1]);
//     cv::line(plot_image, pre_point, cur_point, cv::Scalar(0, 0, 0));  
//     pre_point = cur_point;
//   }
//   //VO
//   pre_point = cv::Point(simulator.init_poses()[0].initPose.translation()[0], simulator.init_poses()[0].initPose.translation()[1]);
//   for (size_t i=1; i<simulator.poses().size(); ++i) {
//     cur_point = cv::Point(simulator.init_poses()[i].initPose.translation()[0], simulator.init_poses()[i].initPose.translation()[1]);
//     cv::line(plot_image, pre_point, cur_point, cv::Scalar(0, 255, 0));  
//     pre_point = cur_point;
//   }
//   
//   //Window
//   pre_point = cv::Point(stateWindowVector[0].transpose()[0], stateWindowVector[0].transpose()[1]);
//   for (size_t i=1; i<simulator.poses().size(); ++i) {
//     cur_point = cv::Point(stateWindowVector[i].transpose()[0], stateWindowVector[i].transpose()[1]);
//     cv::line(plot_image, pre_point, cur_point, cv::Scalar(0, 0, 255));  
//     pre_point = cur_point;
//   }
//   //Batch
//   pre_point = cv::Point(stateBatchVector[0].transpose()[0], stateBatchVector[0].transpose()[1]);
//   for (size_t i=1; i<simulator.poses().size(); ++i) {
//     cur_point = cv::Point(stateBatchVector[i].transpose()[0], stateBatchVector[i].transpose()[1]);
//     cv::line(plot_image, pre_point, cur_point, cv::Scalar(255, 0, 255));  
//     pre_point = cur_point;
//   }
//   //PnP
//   pre_point = cv::Point(statePnPVector[0].transpose()[0], statePnPVector[0].transpose()[1]);
//   for (size_t i=1; i<simulator.poses().size(); ++i) {
//     cur_point = cv::Point(statePnPVector[i].transpose()[0], statePnPVector[i].transpose()[1]);
//     cv::line(plot_image, pre_point, cur_point, cv::Scalar(227,207,87));  
//     pre_point = cur_point;
//   }
//   
//   cv::imshow( "Test", plot_image );  
//   cv::waitKey( 0 );  
//   ***************************************************************
 //LED 3D 布局示意图
//   char filenametemp_LED[]="LED3D布局示意图.txt";
//   fout.open( filenametemp_LED );           //不能有空格
//   for (size_t i=0; i<simulator.landmarks().size(); ++i) {
//     fout.precision(10); //设置有效数字10位
//     fout<<simulator.landmarks()[i].xyzLED[0]<<" "
// 	   <<simulator.landmarks()[i].xyzLED[1]<<" "
// 	   <<simulator.landmarks()[i].xyzLED[2]<<" ";
//     fout<<endl;
//   }
//   fout<<flush;
//   fout.close();
//   printf("LED3D布局示意图文件输出完毕\r\n");
// 
//   char filenametemp_truepose[]="truepose.txt";
//   fout.open(filenametemp_truepose);           
//   for (size_t i=0; i<numOfTotalFrame; ++i) {
//     fout.precision(10); //设置有效数字10位
//     fout<<simulator.poses()[i].truePose.translation()[0]<<" "
// 	  <<simulator.poses()[i].truePose.translation()[1]<<" "
// 	  <<simulator.poses()[i].truePose.translation()[2]<<" ";
//     fout<<endl;
//   }
//   fout<<flush;
//   fout.close();
//   printf("truepose文件输出完毕\r\n");
//   
//   char filenametemp_window[]="window.txt";
//   fout.open(filenametemp_window);           
//   for (size_t i=0; i<numOfTotalFrame; ++i) {
//     fout.precision(10); //设置有效数字10位
//     fout<<stateWindowVector[i][0]<<" "
// 	   <<stateWindowVector[i][1]<<" "
// 	   <<stateWindowVector[i][2]<<" ";
//     fout<<endl;
//   }
//   fout<<flush;
//   fout.close();
//   printf("window文件输出完毕\r\n");
// 
//   char filenametemp_batch[]="batch.txt";
//   fout.open(filenametemp_batch);           
//   for (size_t i=0; i<numOfTotalFrame; ++i) {
//     fout.precision(10); //设置有效数字10位
//     fout<<stateBatchVector[i][0]<<" "
// 	   <<stateBatchVector[i][1]<<" "
// 	   <<stateBatchVector[i][2]<<" ";
//     fout<<endl;
//   }
//   fout<<flush;
//   fout.close();
//   printf("batch文件输出完毕\r\n");
//   
//   char filenametemp_pnp[]="pnp.txt";
//   fout.open(filenametemp_pnp);           
//   for (size_t i=0; i<numOfTotalFrame; ++i) {
//     fout.precision(10); //设置有效数字10位
//     fout<<statePnPVector[i][0]<<" "
// 	   <<statePnPVector[i][1]<<" "
// 	   <<statePnPVector[i][2]<<" ";
//     fout<<endl;
//   }
//   fout<<flush;
//   fout.close();
//   printf("pnp文件输出完毕\r\n");
//   ×××××××××××××××××××××××××××××××××××××××××××××××××××××××××

 
//   for ( size_t i=0; i<states.size(); ++i ) {
//     int num_of_visible_led = generateCameraMeas( states[i],  ledLayout, measNoise, cam);
//     cout<<"状态： "<<states[i].pose.translation().transpose();
//     cout<<" 可见到 "<<num_of_visible_led<<"个LED";
//     cout<<" 初始状态： "<<states[i].init_pose.translation().transpose();
//     cout<<" 里程计： "<<states[i].odometry.translation().transpose()<<endl;
//   }
  