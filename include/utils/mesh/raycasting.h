/*
 *   raycasting on mesh - wrapper towards the CGAL
 */

#ifndef RAYCASTING_H
#define RAYCASTING_H

#include <Eigen/Dense>

#include "utils/convertors/igl_to_cgal.h"

// headers specific to raycasting
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Ray_3 Ray;
typedef K::Point_3 Point;
typedef CGAL::AABB_face_graph_triangle_primitive<CGAL_Mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

class Raycasting {
    private: 
        Eigen::MatrixXd V_;
        Eigen::MatrixXi F_;
        CGAL_Mesh mesh;
        std::shared_ptr<Tree> tree;
        
    public: 
        Raycasting(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
        bool cast_ray(Eigen::Vector3d source, Eigen::Vector3d vector, Eigen::Vector3d & hit);
        bool cast_ray(Eigen::Vector3d source, Eigen::Vector3d vector, double & distance);
};


Raycasting::Raycasting(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) : V_(V), F_(F)
{
    convert_mesh_from_IGL_to_CGAL(V_, F_, mesh);
    this->tree = std::make_shared< Tree >(faces(mesh).first, faces(mesh).second, mesh);

    // call to random raycasting to force the initialization of the AABB tree
    Eigen::Vector3d a, b, c;
    a << 0, 0, 0;
    b << 1, 1, 1;
    cast_ray(a, b, c);
}

bool Raycasting::cast_ray(Eigen::Vector3d source, Eigen::Vector3d vector, Eigen::Vector3d & hit)
{
    // output as tuple? https://www.educative.io/edpresso/how-to-return-multiple-values-from-a-function-in-cpp17
    bool success = false;

    Point a(source(0), source(1), source(2));
    Point b(source(0) + vector(0), source(1) + vector(1), source(2) + vector(2));
    Ray ray_query(a, b);

    auto intersection = this->tree->first_intersection(ray_query);

    if(intersection)
    {
        success = true;
        if(boost::get<Point>(&(intersection->first)))
        {
            const Point* p =  boost::get<Point>(&(intersection->first));
            hit(0) = p[0][0];
            hit(1) = p[0][1];
            hit(2) = p[0][2];
        }
    }

    return success;
}

bool Raycasting::cast_ray(Eigen::Vector3d source, Eigen::Vector3d vector, double & distance)
{
    Eigen::Vector3d hit;
    bool success = cast_ray(source, vector, hit);
    if (success)
        distance = (source - hit).norm();
    else
        distance = -1;

    return success;
}

#endif
