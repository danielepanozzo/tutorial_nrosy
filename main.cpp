#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/local_basis.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/viewer/Viewer.h>
#include <igl/triangle_triangle_adjacency.h>

#include "nrosy.h"

// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// Triangle-triangle adjacency
Eigen::MatrixXi TT;
Eigen::MatrixXi TTi;

// Constrained faces id
Eigen::VectorXi b;

// Cosntrained faces representative vector
Eigen::MatrixXd bc;

// Degree of the N-RoSy field
int N = 4;

// Converts a representative vector per face in the full set of vectors that describe
// an N-RoSy field
void representative_to_nrosy(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& R,
  const int N,
  Eigen::MatrixXd& Y)
{
  using namespace Eigen;
  using namespace std;
  MatrixXd B1, B2, B3;

  igl::local_basis(V,F,B1,B2,B3);

  Y.resize(F.rows()*N,3);
  for (unsigned i=0;i<F.rows();++i)
  {
    double x = R.row(i) * B1.row(i).transpose();
    double y = R.row(i) * B2.row(i).transpose();
    double angle = atan2(y,x);

    for (unsigned j=0; j<N;++j)
    {
      double anglej = angle + 2*M_PI*double(j)/double(N);
      double xj = cos(anglej);
      double yj = sin(anglej);
      Y.row(i*N+j) = xj * B1.row(i) + yj * B2.row(i);
    }
  }
}

// Plots the mesh with an N-RoSy field
// The constrained faces (b) are colored in red.
void plot_mesh_nrosy(
  igl::viewer::Viewer& viewer,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  int N,
  Eigen::MatrixXd& PD1,
  Eigen::VectorXi& b)
{
  using namespace Eigen;
  using namespace std;
  // Clear the mesh
  viewer.data.clear();
  viewer.data.set_mesh(V,F);

  // Expand the representative vectors in the full vector set and plot them as lines
  double avg = igl::avg_edge_length(V, F);
  MatrixXd Y;
  representative_to_nrosy(V, F, PD1, N, Y);

  MatrixXd B;
  igl::barycenter(V,F,B);

  MatrixXd Be(B.rows()*N,3);
  for(unsigned i=0; i<B.rows();++i)
    for(unsigned j=0; j<N; ++j)
      Be.row(i*N+j) = B.row(i);

  viewer.data.add_edges(Be,Be+Y*(avg/2),RowVector3d(0,0,1));

  // Highlight in red the constrained faces
  MatrixXd C = MatrixXd::Constant(F.rows(),3,1);
  for (unsigned i=0; i<b.size();++i)
    C.row(b(i)) << 1, 0, 0;
  viewer.data.set_colors(C);
}

  // It allows to change the degree of the field when a number is pressed
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
  using namespace Eigen;
  using namespace std;
  if (key >= '1' && key <= '9')
    N = key - '0';

  MatrixXd R = nrosy(V,F,TT,TTi,b,bc,N);


  plot_mesh_nrosy(viewer,V,F,N,R,b);

  return false;
}

int main(int argc, char *argv[])
{
  using namespace std;
  using namespace Eigen;

  // Load a mesh in OBJ format
  igl::readOBJ("../2triangles.obj", V, F);

  // Triangle-triangle adjacency
  igl::triangle_triangle_adjacency(F,TT,TTi);

  // Simple contraint
  b.resize(1);
  b << 0;
  bc.resize(1,3);
  bc << 1,1,1;

  igl::viewer::Viewer viewer;

  // Interpolate the field and plot
  key_down(viewer, '1', 0);

  // Plot the mesh
  viewer.data.set_mesh(V, F);
  viewer.callback_key_down = &key_down;

  // Disable wireframe
  viewer.core.show_lines = false;

  // Launch the viewer
  viewer.launch();
}
