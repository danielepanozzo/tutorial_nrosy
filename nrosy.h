#ifndef NROSY
#define NROSY

#include <Eigen/Core>

Eigen::MatrixXd nrosy
        (
        const Eigen::MatrixXd& V,          // Vertices of the mesh
        const Eigen::MatrixXi& F,          // Faces
        const Eigen::MatrixXi& TT,         // Adjacency triangle-triangle
        const Eigen::VectorXi& soft_id,    // Soft constraints face ids
        const Eigen::MatrixXd& soft_value, // Soft constraints 3d vectors
        const int n                        // Degree of the n-rosy field
        );

#endif
