/*
*   delete vertices from mesh and update the faces
*/

#ifndef DELETE_VERTICES_H
#define DELETE_VERTICES_H

#include <Eigen/Core>

#include <vector>
#include <algorithm>

inline std::vector<int> vector_complement(std::vector<int> v, int size) {

    std::vector<int> complement(size-v.size());

    if (v.size() > 0) {
        // sort and remove duplicates
        std::sort (v.begin(), v.end());
        v.erase( unique( v.begin(), v.end() ), v.end() );

        // populate the complement
        int counter = 0;
        for (int j=0; j<v[0]; j++) {
            complement[counter] = j;
            counter ++;
        }
        if (v.size() > 1) {
            for (int i=0; i<v.size()-1; i++) {
                for (int j=v[i]+1; j<v[i+1]; j++) {
                    complement[counter] = j;
                    counter ++;
                }
            }
        }
        for (int j=v.back()+1; j<size; j++) {
            complement[counter] = j;
            counter ++;
        }
    } else {
        std::iota(complement.begin(), complement.end(), 0);
    }

    return complement;
}

inline void deleteFacesConnectedToVertices(Eigen::MatrixXi& F, 
                                           std::vector<int>& vertices_to_delete,
                                           int number_of_vertices) {
    // build adjacency of faces
    std::vector < std::vector <int> > face_adjacency(number_of_vertices);
    for (int i=0; i<F.cols(); i++){
        for (int j=0; j<F.rows(); j++) {
            face_adjacency.at( F(j,i) ).push_back(i);
        }
    }

    // build list of faces to delete
    std::vector<int> faces_to_delete;
    for (auto vertex_id : vertices_to_delete) {
        for (auto face_id : face_adjacency.at(vertex_id)) {
            faces_to_delete.push_back(face_id);
        }
    }

    // make unique
    sort( faces_to_delete.begin(), faces_to_delete.end() );
    faces_to_delete.erase( unique( faces_to_delete.begin(), faces_to_delete.end() ), faces_to_delete.end() );

    std::vector<int> faces_to_keep = vector_complement(faces_to_delete, F.cols());

    // update the faces
    Eigen::MatrixXi F_temp;
    F_temp.resize(3, faces_to_keep.size());
    for (int i=0; i<faces_to_keep.size(); i++) {
        F_temp.col(i) = F.col(faces_to_keep[i]);
    }
    F = F_temp;

    // update the faces indexes
    std::vector<int> vertices_to_keep = vector_complement(vertices_to_delete, number_of_vertices);
    for (int i=0; i<F.cols(); i++) {
        for (int j=0; j<F.rows(); j++) {
            std::vector<int>::iterator itr = std::find(vertices_to_keep.begin(), vertices_to_keep.end(), F(j,i) );
            F(j,i) = std::distance(vertices_to_keep.begin(), itr);
        }
    }
}

inline void deleteVertices(Eigen::MatrixXd& V,
                           Eigen::MatrixXi& F, 
                           Eigen::MatrixXd& N, 
                           Eigen::MatrixXi& C, 
                           std::vector<int> vertices_to_delete) {

    if (F.cols() != 0)
        deleteFacesConnectedToVertices(F, vertices_to_delete, V.cols());

    std::vector<int> vertices_to_keep = vector_complement(vertices_to_delete, V.cols());

    // update the normals
    if (N.cols() == V.cols()) {
        Eigen::MatrixXd N_temp;
        N_temp.resize(3, vertices_to_keep.size());
        for (int i=0; i<vertices_to_keep.size(); i++) {
            N_temp.col(i) = N.col(vertices_to_keep[i]);
        }
        N = N_temp;
    }

    // update the colors
    if (C.cols() == V.cols()) {
        Eigen::MatrixXi C_temp;
        C_temp.resize(3, vertices_to_keep.size());
        for (int i=0; i<vertices_to_keep.size(); i++) {
            C_temp.col(i) = C.col(vertices_to_keep[i]);
        }
        C = C_temp;
    }

    // update the vertices
    Eigen::MatrixXd V_temp;
    V_temp.resize(3, vertices_to_keep.size());
    for (int i=0; i<vertices_to_keep.size(); i++) {
        V_temp.col(i) = V.col(vertices_to_keep[i]);
    }
    V = V_temp;

}

inline void deleteVertices(Eigen::MatrixXd& V,
                           Eigen::MatrixXi& F, 
                           std::vector<int> vertices_to_delete) {
    Eigen::MatrixXd N;
    Eigen::MatrixXi C;
    deleteVertices(V, F, N, C, vertices_to_delete);
}

void deleteVerticesWithoutColor(Eigen::Vector3i color, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &N, Eigen::MatrixXi &C) {
    std::vector<int> vertices_to_delete;
    for (int i=0; i<V.cols(); i++) {
        if (C.col(i) != color) 
            vertices_to_delete.push_back(i);
    }
    deleteVertices(V, F, N, C, vertices_to_delete);
}


#endif
