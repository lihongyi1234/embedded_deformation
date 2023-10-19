

int convert_mesh_from_PCL_to_CGAL(pcl::PolygonMesh::Ptr PCL_mesh, Eigen::MatrixXd& V, Eigen::MatrixXd& F)
{
    pcl::PointCloud<pcl::PointXYZ>::Ptr mesh_cloud (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::fromPCLPointCloud2( PCL_mesh->cloud, *mesh_cloud );


    // get the vertices
    Eigen::MatrixXf M = mesh_cloud->getMatrixXfMap();

    // clear and reserve the memory
    CGAL_mesh.clear();
    int n = mesh_cloud->size();
    int f = PCL_mesh->polygons.size();
    int e = 0;
    CGAL_mesh.reserve(n, 2*f, e);

    //copy the vertices
    double x, y, z;
    for (int i=0; i<mesh_cloud->size(); i++)
    {
        Point p;
        x = mesh_cloud->points[i].x;
        y = mesh_cloud->points[i].y;
        z = mesh_cloud->points[i].z;
        p = Point(x, y, z);
        CGAL_mesh.add_vertex(p);
    }

    // copy the faces
    std::vector <int> vertices;
    for(int i=0; i<PCL_mesh->polygons.size(); i++)
    {
        vertices.resize(3);
        vertices[0] = PCL_mesh->polygons[i].vertices[0];
        vertices[1] = PCL_mesh->polygons[i].vertices[1];
        vertices[2] = PCL_mesh->polygons[i].vertices[2];
        CGAL_mesh.add_face(CGAL_Mesh::Vertex_index (vertices[0]), 
                           CGAL_Mesh::Vertex_index (vertices[1]),
                           CGAL_Mesh::Vertex_index (vertices[2]));
    }

    return 0;
}