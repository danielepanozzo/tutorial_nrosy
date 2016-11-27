#include "nrosy.h"
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <iostream>

using namespace Eigen;

MatrixXd nrosy
        (
        const MatrixXd& V,          // Vertices of the mesh
        const MatrixXi& F,          // Faces
        const MatrixXi& TT,         // Adjacency triangle-triangle
        const MatrixXi& TTi,        // Adjacency triangle-triangle-index
        const VectorXi& soft_id,    // Soft constraints face ids
        const MatrixXd& soft_value, // Soft constraints 3d vectors
        const int n                 // Degree of the n-rosy field
        )
{
  assert(soft_id.size() > 0); // One constraint is necessary to make the solution unique

  std::cerr << V << std::endl;
  std::cerr << F << std::endl;

  Matrix<double,Dynamic,3> T1(F.rows(),3), T2(F.rows(),3), T3(F.rows(),3);

  for (unsigned i=0;i<F.rows();++i)
  {
    // Computes the local reference systems for each face
    Vector3d e1 =  V.row(F(i, 1)) - V.row(F(i, 0));
    Vector3d e2 =  V.row(F(i, 2)) - V.row(F(i, 0));
    T1.row(i) = e1.normalized();
    T3.row(i) = T1.row(i).cross(e2).normalized();
    T2.row(i) = T1.row(i).cross(T3.row(i)).normalized();
  }

  // Build the sparse matrix, with an energy term for each edge
  std::vector< Triplet<std::complex<double> > > t;
  std::vector< Triplet<std::complex<double> > > tb;
  t.reserve(V.rows()*10);
  tb.reserve(soft_id.size());

  int count = 0;

  for (int fi=0;fi<F.rows();++fi)
  {
    for (int ei=0;ei<F.cols();++ei)
    {
      // Look up the opposite face
      int fj = TT(fi,ei);
      int ej = TTi(fi,ei);

      // If it is a boundary edge, it does not contribute to the energy
      if ((fi == -1) || (fj == -1)) continue;

      // Compute the complex representation of the common edge
      Vector3d vi = V.row(F(fi,(ei+1)%3)) - V.row(F(fi,ei));
      std::complex<double> ci(vi.dot(T1.row(fi)),vi.dot(T2.row(fi)));
      // Note that the edge is flipped for fj to account for the different sorting of the vertices
      Vector3d vj = V.row(F(fj,ej)) - V.row(F(fj,(ej+1)%3));
      std::complex<double> cj(vj.dot(T1.row(fj)),vj.dot(T2.row(fj)));

      // Energy term is ui^n*ei^n - uj^n*ej^n
      t.push_back(Triplet<std::complex<double> >(count,fi, std::pow(std::conj(ci),n)));
      t.push_back(Triplet<std::complex<double> >(count,fj,-std::pow(std::conj(cj),n)));

      ++count;
    }
  }

  // Convert the constraints into the complex polynomial coefficients and add them as soft constraints
  for (int r=0; r<soft_id.size(); ++r)
  {
    int f = soft_id(r);
    Vector3d v = soft_value.row(r);
    std::complex<double> c(v.dot(T1.row(f)),v.dot(T2.row(f)));
    std::cerr << "Const: " << std::pow(c,n) << std::endl;
    t.push_back(Triplet<std::complex<double> >(count,f, 100));
    tb.push_back(Triplet<std::complex<double> >(count,0, std::pow(c,n) * std::complex<double>(100,0)));
    count++;
  }

  for(int i=0; i<t.size(); ++i)
    std::cerr << t[i].row() << "\t\t" << t[i].col() << "\t\t" << t[i].value() << std::endl;

  // Solve the linear system
  SparseMatrix<std::complex<double>,RowMajor> A(count,F.rows());
  A.setFromTriplets(t.begin(), t.end());
  SparseMatrix<std::complex<double>,RowMajor> b(count,1);
  b.setFromTriplets(tb.begin(), tb.end());

  SparseLU< SparseMatrix<std::complex<double> > > solver;
  solver.compute(A.transpose()*A);
  assert(solver.info()==Success);

  Matrix<std::complex<double>, Dynamic, Dynamic> x = solver.solve(A.transpose()*MatrixXcd(b));
  assert(solver.info()==Success);

  std::cerr << "soft_id:" << std::endl << soft_id << std::endl;
  //std::cerr << "soft_value:" << std::endl << soft_value << std::endl;


  std::cerr << "A:" << std::endl << MatrixXcd(A) << std::endl;

  std::cerr << "b:" << std::endl << MatrixXcd(b) << std::endl;

  std::cerr << "x:" << std::endl << MatrixXcd(x) << std::endl;

  // Convert the interpolated polyvector into Euclidean vectors
  MatrixXd R(F.rows(),3);
  for (int f=0; f<F.rows(); ++f)
  {
    std::complex<double> c = std::pow(x(f),-n);
    R.row(f) = T1.row(f) * c.real() + T2.row(f) * c.imag();
  }

  std::cerr << "R:" << std::endl << R << std::endl;

  return R;


}

