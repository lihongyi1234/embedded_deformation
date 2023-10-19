#ifndef TO_PCL
#define TO_PCL

/*
 *   convertion towards PCL
 */

#include <Eigen/Core>
#include <pcl/point_types.h>
#include <pcl/common/common.h>

pcl::PointCloud<pcl::PointXYZ>::Ptr to_PCL(Eigen::MatrixXd& V)
{
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);

    cloud->points.resize( V.cols() );
    for(size_t i=0; i<V.cols(); ++i)
        cloud->points[i].getVector3fMap() = V.row(i).cast<float>();

    return cloud;
}

pcl::PointCloud<pcl::Normal>::Ptr to_PCL_normals(Eigen::MatrixXd& N)
{
    pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);

    // set points
    normals->points.resize( N.cols() );
    for (int i = 0; i < N.cols(); i++)
    {
        normals->points[i].normal_x = N(i,0);
        normals->points[i].normal_y = N(i,1);
        normals->points[i].normal_z = N(i,2);
    }

    return normals;
}

pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr to_PCL(Eigen::MatrixXd& V, Eigen::MatrixXd& N, Eigen::MatrixXi& C)
{
    pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);

    // set points
    cloud->points.resize( V.cols() );
    for (int i = 0; i < V.cols(); i++)
    {
        cloud->points[i].x = V(i,0);
        cloud->points[i].y = V(i,1);
        cloud->points[i].z = V(i,2);

        cloud->points[i].r = C(i,0);
        cloud->points[i].g = C(i,1);
        cloud->points[i].b = C(i,2);

        cloud->points[i].curvature = 0;

        cloud->points[i].normal_x = N(i,0);
        cloud->points[i].normal_y = N(i,1);
        cloud->points[i].normal_z = N(i,2);
    }

    return cloud;
}

#endif
