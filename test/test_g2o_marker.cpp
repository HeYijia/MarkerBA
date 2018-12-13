#include"g2otypes_marker.h"
#include "camera_model.h"
#include "frame.h"

#include <g2o/core/sparse_optimizer.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/solver.h>
#include <g2o/core/robust_kernel_impl.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <g2o/types/sba/types_six_dof_expmap.h>


void optWithMarker()
{
  Camera_Model::PinholeCamera* pinhole;

  pinhole = new Camera_Model::PinholeCamera(640,480,525.0,525.0, 275.0, 228.0);
  double tag_size_ = 0.2;  // 20cm

  // create a marker
  markerslam::Marker* tag = new markerslam::Marker(tag_size_);
  // set marker in world frame
  tag->T_w_m = Eigen::Matrix4d::Identity();
  tag->T_w_m(1,1) = -1;
  tag->T_w_m(2,2) = -1;
  tag->T_w_m(2,3) = 2;

  //create two keyframe
  markerslam::FramePtr kf1( new markerslam::Frame(1));
  kf1->T_w_c_ = Eigen::Matrix4d::Identity();

  markerslam::FramePtr kf2( new markerslam::Frame(2));
  kf2->T_w_c_ = Eigen::Matrix4d::Identity();
  kf2->T_w_c_(0,3) = 1;


  // setup g2o
  g2o::SparseOptimizer optimizer;
  optimizer.setVerbose(false);

  g2o::BlockSolver_6_3::LinearSolverType * linearSolver;
  linearSolver = new g2o::LinearSolverCholmod< g2o::BlockSolver_6_3::PoseMatrixType >();
  g2o::BlockSolver_6_3* solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

  g2o::OptimizationAlgorithmLevenberg * solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);

  solver->setMaxTrialsAfterFailure(5);
  optimizer.setAlgorithm(solver);

  //set keyframe vertices , kf1
  g2o::VertexSE3Expmap * vSE3_1 = new g2o::VertexSE3Expmap();
  vSE3_1->setEstimate( g2o::SE3Quat(kf1->TcwRotation(),kf1->TcwPosition()) );
  vSE3_1->setId(1);
  vSE3_1->setFixed(true);
  optimizer.addVertex(vSE3_1);

  // kf2 vertices
  g2o::VertexSE3Expmap * vSE3_2 = new g2o::VertexSE3Expmap();
  //Eigen::Vector3d kf2_position_noise(0.5,0.03,0.02);
  Eigen::Vector3d kf2_position_noise(0.2,0.3,0.1);
  vSE3_2->setEstimate( g2o::SE3Quat(kf2->TcwRotation(),kf2->TcwPosition()+kf2_position_noise) );
  vSE3_2->setId(2);
  vSE3_2->setFixed(false);
  optimizer.addVertex(vSE3_2);

  // set marker vertices
  g2o::VertexSE3Expmap * vSE3_marker = new g2o::VertexSE3Expmap();
  Eigen::Vector3d marker_position_noise(0.5,0.12,0.02);

  Eigen::Matrix3d Rnoise( Eigen::AngleAxisd(0.1, Eigen::Vector3d::UnitZ()) );
//  std::cout << "R noise \n" << Rnoise * tag->TmwRotation() << std::endl;

  vSE3_marker->setEstimate( g2o::SE3Quat( Rnoise * tag->TmwRotation(),tag->TmwPosition() + marker_position_noise ) );
  vSE3_marker->setId(3);
  vSE3_marker->setFixed(false);
  optimizer.addVertex(vSE3_marker);

  // set g2o edges
  // every marker have four conner
  for(int i = 0; i < 4; i++)
  {

    // Transform the conner coordinate to camera frame
    Eigen::Vector4d pm( tag->pxy_m.col(i));
    Eigen::Vector4d pw = tag->T_w_m * pm;
    Eigen::Vector4d pc1 = kf1->T_w_c_.inverse() * pw;
    Eigen::Vector4d pc2 = kf2->T_w_c_.inverse() * pw;

    // projectio to camera image plane, create pixel observation
    Eigen::Vector2d obs1 = pinhole->world2cam(Eigen::Vector3d(pc1(0),pc1(1),pc1(2)));  // pixel coordinate in kf1 img
    Eigen::Vector2d obs2 = pinhole->world2cam(Eigen::Vector3d(pc2(0),pc2(1),pc2(2)));  // pixel coordinate in kf2 img

    // create Edge between marker and kf1
    g2o::EdgeSE3ProjectMarker *e1 = new g2o::EdgeSE3ProjectMarker();
    e1->fx = pinhole->fx();
    e1->fy = pinhole->fy();
    e1->cx = pinhole->cx();
    e1->cy = pinhole->cy();

    e1->setVertex(0, dynamic_cast< g2o::OptimizableGraph::Vertex* >(vSE3_1));
    e1->setVertex(1,dynamic_cast<g2o::OptimizableGraph::Vertex* >(vSE3_marker) );

    e1->setMeasurement(obs1);
    e1->SetMarkerPoint(Eigen::Vector3d(pm(0),pm(1),pm(2)));
    e1->setInformation(Eigen::Matrix2d::Identity());
    g2o::RobustKernelHuber *rk = new g2o::RobustKernelHuber;
    rk->setDelta(2.0);   // 2 pixel
    e1->setRobustKernel(rk);

    optimizer.addEdge(e1);

    // edge between marker and kf2
    g2o::EdgeSE3ProjectMarker *e2 = new g2o::EdgeSE3ProjectMarker();
    e2->fx = pinhole->fx();
    e2->fy = pinhole->fy();
    e2->cx = pinhole->cx();
    e2->cy = pinhole->cy();

    e2->setVertex(0, dynamic_cast< g2o::OptimizableGraph::Vertex* >(vSE3_2));
    e2->setVertex(1,dynamic_cast<g2o::OptimizableGraph::Vertex* >(vSE3_marker) );
    e2->setMeasurement(obs2);
    e2->SetMarkerPoint(Eigen::Vector3d(pm(0),pm(1),pm(2)));
    e2->setInformation(Eigen::Matrix2d::Identity());
    g2o::RobustKernelHuber *rk2 = new g2o::RobustKernelHuber;
    rk2->setDelta(2.0);   // 2 pixel
    e2->setRobustKernel(rk2);

    optimizer.addEdge(e2);

  }


  //g2o::VertexSE3Expmap* vSE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(3));
  //std::cout<<" before kf1 pose: "<<std::endl<<vSE3_1->estimate()<<std::endl;
  std::cout<<"before optimizer, kf2 pose: "<<std::endl<<vSE3_2->estimate().inverse()<<std::endl;
  std::cout<<"before optimizer, marker pose: "<<std::endl<<vSE3_marker->estimate().inverse() << std::endl;

  optimizer.initializeOptimization();
  optimizer.optimize(3);

  //g2o::VertexSE3Expmap* vSE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(3));
  //std::cout<<"kf1 pose: "<<std::endl<<vSE3_1->estimate()<<std::endl;
  std::cout<<"after optimizer, kf2 pose: "<<std::endl<<vSE3_2->estimate().inverse()<<std::endl;
  std::cout<<"the ground truth kf2 pose: "<<std::endl<<kf2->T_w_c_<<std::endl;

  std::cout<<"after optimizer, marker pose: "<<std::endl<<vSE3_marker->estimate().inverse() << std::endl;
  std::cout<<"the ground truth marker pose: "<<std::endl<<tag->T_w_m << std::endl;
  //std::cout<<"marker pose, after optimizer: "<<std::endl<<vSE3->estimate() << std::endl;
}

int main()
{

  std::cout<<"hello ba"<<std::endl;

  optWithMarker();

  return 0;
}
