#include"g2otypes_marker.h"
#include "g2o/core/factory.h"
#include "g2o/stuff/macros.h"

namespace g2o {

using namespace std;


Vector2d project2d(const Vector3d& v)  {
  Vector2d res;
  res(0) = v(0)/v(2);
  res(1) = v(1)/v(2);
  return res;
}

Vector3d unproject2d(const Vector2d& v)  {
  Vector3d res;
  res(0) = v(0);
  res(1) = v(1);
  res(2) = 1;
  return res;
}
/*
VertexSE3Expmap::VertexSE3Expmap() : BaseVertex<6, SE3Quat>() {
}

bool VertexSE3Expmap::read(std::istream& is) {
  Vector7d est;
  for (int i=0; i<7; i++)
    is  >> est[i];
  SE3Quat cam2world;
  cam2world.fromVector(est);
  setEstimate(cam2world.inverse());
  return true;
}

bool VertexSE3Expmap::write(std::ostream& os) const {
  SE3Quat cam2world(estimate().inverse());
  for (int i=0; i<7; i++)
    os << cam2world[i] << " ";
  return os.good();
}
*/

EdgeSE3ProjectMarker::EdgeSE3ProjectMarker() : BaseBinaryEdge<2, Vector2d, VertexSE3Expmap, VertexSE3Expmap>() {
}

bool EdgeSE3ProjectMarker::read(std::istream& is){
  for (int i=0; i<2; i++){
    is >> _measurement[i];
  }
  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }
  return true;
}

bool EdgeSE3ProjectMarker::write(std::ostream& os) const {

  for (int i=0; i<2; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}


void EdgeSE3ProjectMarker::linearizeOplus() {
  VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[1]);
  SE3Quat Tmw(vj->estimate());
  VertexSE3Expmap* vi = static_cast<VertexSE3Expmap*>(_vertices[0]);
  SE3Quat Tcw(vi->estimate());

  SE3Quat Tcm(Tcw * Tmw.inverse());
  //Vector3d xyz = point3d_marker;
  Vector3d point3d_camera = Tcm.map(point3d_marker);    //  transfrom 3d point from marker frame to camera frame

  double x = point3d_camera[0];
  double y = point3d_camera[1];
  double z = point3d_camera[2];
  double z_2 = z*z;

  // jacobian used to update marker pose
  Matrix<double,2,3> tmp;
  tmp(0,0) = fx;
  tmp(0,1) = 0;
  tmp(0,2) = -x/z*fx;

  tmp(1,0) = 0;
  tmp(1,1) = fy;
  tmp(1,2) = -y/z*fy;

  Matrix<double,3,6> partital_PctoLie;
  double xm = point3d_marker[0];
  double ym = point3d_marker[1];
  double zm = point3d_marker[2];
  partital_PctoLie(0,0)=0;  partital_PctoLie(0,1)=-zm;  partital_PctoLie(0,2)=ym;  partital_PctoLie(0,3)=-1;  partital_PctoLie(0,4)=0;  partital_PctoLie(0,5)=0;
  partital_PctoLie(1,0)=zm;  partital_PctoLie(1,1)=0;  partital_PctoLie(1,2)=-xm;  partital_PctoLie(1,3)=0;  partital_PctoLie(1,4)=-1;  partital_PctoLie(1,5)=0;
  partital_PctoLie(2,0)=-ym;  partital_PctoLie(2,1)=xm;  partital_PctoLie(2,2)=0;  partital_PctoLie(2,3)=0;  partital_PctoLie(2,4)=0;  partital_PctoLie(2,5)=-1;

  _jacobianOplusXj =  -1./z * tmp * Tcm.rotation().toRotationMatrix()*partital_PctoLie;

  // jacobian used to update  camera pose
  _jacobianOplusXi(0,0) =  x*y/z_2 *fx;
  _jacobianOplusXi(0,1) = -(1+(x*x/z_2)) *fx;
  _jacobianOplusXi(0,2) = y/z *fx;
  _jacobianOplusXi(0,3) = -1./z *fx;
  _jacobianOplusXi(0,4) = 0;
  _jacobianOplusXi(0,5) = x/z_2 *fx;

  _jacobianOplusXi(1,0) = (1+y*y/z_2) *fy;
  _jacobianOplusXi(1,1) = -x*y/z_2 *fy;
  _jacobianOplusXi(1,2) = -x/z *fy;
  _jacobianOplusXi(1,3) = 0;
  _jacobianOplusXi(1,4) = -1./z *fy;
  _jacobianOplusXi(1,5) = y/z_2 *fy;

  //std::cout<<"_jacobianOplusXi: "<<_jacobianOplusXi;
}

Vector2d EdgeSE3ProjectMarker::cam_project(const Vector3d & trans_xyz) const{
  Vector2d proj = project2d(trans_xyz);
  Vector2d res;
  res[0] = proj[0]*fx + cx;
  res[1] = proj[1]*fy + cy;
  return res;
}

EdgeSE3ProjectXYZ::EdgeSE3ProjectXYZ() : BaseBinaryEdge<2, Vector2d, VertexSBAPointXYZ, VertexSE3Expmap>() {
}

bool EdgeSE3ProjectXYZ::read(std::istream& is){
  for (int i=0; i<2; i++){
    is >> _measurement[i];
  }
  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }
  return true;
}

bool EdgeSE3ProjectXYZ::write(std::ostream& os) const {

  for (int i=0; i<2; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}


void EdgeSE3ProjectXYZ::linearizeOplus() {
  VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[1]);
  SE3Quat T(vj->estimate());
  VertexSBAPointXYZ* vi = static_cast<VertexSBAPointXYZ*>(_vertices[0]);
  Vector3d xyz = vi->estimate();
  Vector3d xyz_trans = T.map(xyz);

  double x = xyz_trans[0];
  double y = xyz_trans[1];
  double z = xyz_trans[2];
  double z_2 = z*z;

  Matrix<double,2,3> tmp;
  tmp(0,0) = fx;
  tmp(0,1) = 0;
  tmp(0,2) = -x/z*fx;

  tmp(1,0) = 0;
  tmp(1,1) = fy;
  tmp(1,2) = -y/z*fy;

  _jacobianOplusXi =  -1./z * tmp * T.rotation().toRotationMatrix();

  _jacobianOplusXj(0,0) =  x*y/z_2 *fx;
  _jacobianOplusXj(0,1) = -(1+(x*x/z_2)) *fx;
  _jacobianOplusXj(0,2) = y/z *fx;
  _jacobianOplusXj(0,3) = -1./z *fx;
  _jacobianOplusXj(0,4) = 0;
  _jacobianOplusXj(0,5) = x/z_2 *fx;

  _jacobianOplusXj(1,0) = (1+y*y/z_2) *fy;
  _jacobianOplusXj(1,1) = -x*y/z_2 *fy;
  _jacobianOplusXj(1,2) = -x/z *fy;
  _jacobianOplusXj(1,3) = 0;
  _jacobianOplusXj(1,4) = -1./z *fy;
  _jacobianOplusXj(1,5) = y/z_2 *fy;
}

Vector2d EdgeSE3ProjectXYZ::cam_project(const Vector3d & trans_xyz) const{
  Vector2d proj = project2d(trans_xyz);
  Vector2d res;
  res[0] = proj[0]*fx + cx;
  res[1] = proj[1]*fy + cy;
  return res;
}

} // end namespace
