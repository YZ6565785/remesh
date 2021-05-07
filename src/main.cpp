#include <string>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/barycentric_coordinates.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Voronoi_diagram_2.h>

// #include <boost/geometry.hpp>
// #include <boost/geometry/geometries/point_xy.hpp>
// #include <boost/geometry/geometries/polygon.hpp>




#include <typeinfo>


// Delaunay Triangulation
/*
A triangulation is a Delaunay triangulation, if the circumscribing circle of any
facet of the triangulation contains no vertex in its interior
*/
typedef CGAL::Exact_predicates_inexact_constructions_kernel              Kernel;
typedef Kernel::Point_2                                                   Point;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, Kernel>    Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                            Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>                    Delaunay;
typedef Delaunay::Vertex_handle                                   Vertex_handle;
typedef Delaunay::Face_handle                                       Face_handle;

typedef CGAL::Delaunay_triangulation_adaptation_traits_2<Delaunay>           AT;
typedef CGAL::Voronoi_diagram_2<Delaunay, AT>                                VD;

// typedef boost::geometry::model::d2::point_xy<double>                 boostPoint;
// typedef boost::geometry::model::polygon<boostPoint>            boostPolygon;



struct Mesh
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd bnd_uv; // boundary point in uv plane
  Eigen::MatrixXi Fuv; // Faces in uv plane
  Eigen::MatrixXd Parameterized_uv; // points in uv plane
  // after remeshing
  Eigen::MatrixXd Parameterized_uv_new; // remeshed points in uv plane
  Eigen::MatrixXd V_remeshed; // remeshed vertices in 3D
};
class Context
{
  public:
    //magic Eigen3 macro : https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Mesh mesh;

    // GUI states
    int    num_vex;
    float  nv_len;
    float  point_size;
    float  line_width;
    int    sel_vidx;
    int    mode;


    bool   show_para = false;
    bool   show_texture = false;
    bool   show_original = true;

	// bool show_mesh;
	// bool show_normals;
};

Context g_ctx;


int g_cur_ind = -0;





////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS
// callback function of key input for igl viewer
bool
keyDown (igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier);

// reset the mesh in viewer
void
resetDisplay (igl::opengl::glfw::Viewer& viewer);

// print key input instructions
void
printHelp ();

////////////////////////////////////////////////////////////////////////////////
// ALGORITHMS FUNCTIONS
// remeshing operation
void
doRemeshing (igl::opengl::glfw::Viewer& viewer, Context & g_ctx);

// parameterization
void
doParameterization (Context & g_ctx);

// Lloyd relaxation
void
iterativeLloydRelaxation (
  igl::opengl::glfw::Viewer& viewer,
  Context & g_ctx, 
  Delaunay & dt
);

// Lloyd relaxation
void
getBarcentricCoordinates (Context & g_ctx, Delaunay & dt);




struct Point2D
{
    double x;
    double y;
};

Point2D compute2DPolygonCentroid(const Point2D* vertices, int vertexCount)
{
    Point2D centroid = {0, 0};
    double signedArea = 0.0;
    double x0 = 0.0; // Current vertex X
    double y0 = 0.0; // Current vertex Y
    double x1 = 0.0; // Next vertex X
    double y1 = 0.0; // Next vertex Y
    double a = 0.0;  // Partial signed area

    int lastdex = vertexCount-1;
    const Point2D* prev = &(vertices[lastdex]);
    const Point2D* next;

    // For all vertices in a loop
    for (int i=0; i<vertexCount; ++i)
    {
        next = &(vertices[i]);
        x0 = prev->x;
        y0 = prev->y;
        x1 = next->x;
        y1 = next->y;
        a = x0*y1 - x1*y0;
        signedArea += a;
        centroid.x += (x0 + x1)*a;
        centroid.y += (y0 + y1)*a;
        prev = next;
    }

    signedArea *= 0.5;
    centroid.x /= (6.0*signedArea);
    centroid.y /= (6.0*signedArea);

    return centroid;
}



//##############################################################################
//##############################################################################

void
add_points (
  igl::opengl::glfw::Viewer& viewer,
  Eigen::MatrixXd const  & pts_n3,
  Eigen::RowVector3d const  & color
)
{
    //mark points
    viewer.data().add_points(pts_n3, color);

}


////////////////////////////////////////////////////////////////////////////////
void
getBarcentricCoordinates (Context & g_ctx, Delaunay & dt)
{
  std::cout << "Updating the mesh..." << std::endl;
  // prepare new vertices for remeshing in 3D
  g_ctx.mesh.V_remeshed.resize(g_ctx.mesh.V.rows(), 3);

  for (size_t ind = 0; ind < g_ctx.mesh.Parameterized_uv.rows(); ind++)
  {
    // centroids
    // g_ctx.mesh.Parameterized_uv_new(ind);

    auto fh = dt.locate(Point(
      g_ctx.mesh.Parameterized_uv_new(ind,0),
      g_ctx.mesh.Parameterized_uv_new(ind,1)
    ));

    // std::cout << "centroid [" << ind+1 << "] in Face Handle: " 
    //           <<  fh->vertex(0)->info() << " " << fh->vertex(1)->info() << " " 
    //           <<  fh->vertex(2)->info() << std::endl;

    // define the trinagle from the located face
    Eigen::RowVector2d a, b, c;
    a = g_ctx.mesh.Parameterized_uv.row(fh->vertex(0)->info());
    b = g_ctx.mesh.Parameterized_uv.row(fh->vertex(1)->info());
    c = g_ctx.mesh.Parameterized_uv.row(fh->vertex(2)->info());

    Eigen::RowVector3d coords;
    igl::barycentric_coordinates(
      g_ctx.mesh.Parameterized_uv_new.row(ind), // centroids
      a,
      b,
      c,
      coords 
    );
    // std::cout << "###########################" << std::endl;
    // std::cout << coords << std::endl;


    // Eigen::MatrixXd A(
    //   g_ctx.mesh.V.row(fh->vertex(0)->info()),
    //   g_ctx.mesh.V.row(fh->vertex(1)->info()),
    //   g_ctx.mesh.V.row(fh->vertex(2)->info())
    // );
    // std::cout << g_ctx.mesh.V.row(fh->vertex(0)->info()).transpose() << std::endl;
    // std::cout << g_ctx.mesh.V.row(fh->vertex(1)->info()).transpose() << std::endl;
    // std::cout << g_ctx.mesh.V.row(fh->vertex(2)->info()).transpose() << std::endl;
    
    Eigen::MatrixXd A(3,3);
    A.row(0) = g_ctx.mesh.V.row(fh->vertex(0)->info());
    A.row(1) = g_ctx.mesh.V.row(fh->vertex(1)->info());
    A.row(2) = g_ctx.mesh.V.row(fh->vertex(2)->info());

    // std::cout << A << "\n\n";
    // std::cout << A.transpose() << "\n\n";
    // std::cout << coords.transpose() << "\n\n";
    // std::cout << (A.transpose()*coords.transpose()).transpose() << "\n\n";

    // update the vertices in 3D
    g_ctx.mesh.V_remeshed.row(ind) = \
      (A.transpose()*coords.transpose()).transpose();

  }
}

////////////////////////////////////////////////////////////////////////////////
void
iterativeLloydRelaxation (
  igl::opengl::glfw::Viewer& viewer,
  Context & g_ctx, 
  Delaunay & dt
)
{ 

  Eigen::MatrixXd uv_temp = g_ctx.mesh.Parameterized_uv;
  std::vector<double> error_hist;

  for (size_t iter = 0; iter < 10; iter++)
  {
    std::cout << "ITERATION >> " << iter+1 << "      ########################\n";
    std::cout << "ITERATION >> " << iter+1 << "      ########################\n";
    std::cout << "ITERATION >> " << iter+1 << "      ########################\n";
    std::cout << "ITERATION >> " << iter+1 << "      ########################\n";

  
    Delaunay::Finite_vertices_iterator v_itr = dt.finite_vertices_begin();

    std::cout << "\tTotal points " << uv_temp.rows() << std::endl;
    std::cout << "\tThere are " << g_ctx.mesh.bnd_uv.rows() << " boundary points."
              << std::endl;
    std::cout << "\tThere are " 
              << uv_temp.rows()-g_ctx.mesh.bnd_uv.rows()
              << " interior points." << std::endl;

    VD vd(dt);
    assert(vd.is_valid());

    // std::cout << uv_temp << std::endl;

    int count = 0;

    // remeshed uv plane
    g_ctx.mesh.Parameterized_uv_new.resize(uv_temp.rows(), 2);
    int ind = 0;
    bool is_boundary;
    int check_point = ((int)uv_temp.rows()/20==0)?1:(int)uv_temp.rows()/20;
    
    for (; v_itr != dt.finite_vertices_end(); v_itr++)
    { 
      count++;
      if (count%check_point==0)
      {
        // round to two decimals
        float progress = (int)(( (float)count*100./uv_temp.rows() ) * 100 + 0.5);
        progress = (float)progress/100;
        std::cout << progress << "%" <<  std::endl;
      }


      ind = v_itr->info();
      // std::cout << "x=" << v_itr->point().x() << ", y=" << v_itr->point().y()
      //             << std::endl;


      is_boundary = false;
      for (size_t j = 0; j < g_ctx.mesh.bnd_uv.rows(); j++)
      {
        if (v_itr->point().x() == g_ctx.mesh.bnd_uv(j,0) &&
            v_itr->point().y() == g_ctx.mesh.bnd_uv(j,1)
        )
        {
          is_boundary = true;
          break;
        }
      }
      
      if (!is_boundary)
      {
        // std::cout << "Found a interior point." << std::endl;
        auto cell = vd.dual(v_itr);
        auto ec_begin = cell->ccb();

        // std::vector<boostPoint> pts;
        std::vector<Point2D> pts_2;


        auto ec = ec_begin;
        do
        {
          // std::cout << "(" << ec->source()->point().x() << ","
          //           << ec->source()->point().y() << ")" << std::endl;
          {
            Eigen::MatrixXd p (1,2);
            p(0, 0) = ec->source()->point().x();
            p(0, 1) = ec->source()->point().y();
            add_points(
              viewer,
              p,
              Eigen::RowVector3d(1, 0, 0)
            );
          }

          pts_2.push_back(
            (Point2D)
            {
              ec->source()->point().x(),
              ec->source()->point().y()
            }
          );
        } while (++ec != ec_begin);

        {
          Point2D polygon[pts_2.size()];
          int i = 0;
          for (Point2D p: pts_2)
          {
            polygon[i] = {p.x, p.y};
            i++;
          }
          Point2D centroid = compute2DPolygonCentroid(polygon, pts_2.size());
          // std::cout << "Centroid is (" << centroid.x << ", " << centroid.y 
          //           << ")\n";
        
          g_ctx.mesh.Parameterized_uv_new(ind,0) = centroid.x;
          g_ctx.mesh.Parameterized_uv_new(ind,1) = centroid.y;


          // add the centroid to the viewer
          {
            Eigen::MatrixXd p (1,2);
            p(0, 0) = centroid.x;
            p(0, 1) = centroid.y;
            add_points(
              viewer,
              p,
              Eigen::RowVector3d(0, 0, 1)
            );
          }

        }
      }
      else
      {
        g_ctx.mesh.Parameterized_uv_new(ind,0) = v_itr->point().x();
        g_ctx.mesh.Parameterized_uv_new(ind,1) = v_itr->point().y();
      }
    }

    double error = (uv_temp-g_ctx.mesh.Parameterized_uv_new).sum();
    error_hist.push_back(error);
    if (error<1e-10)
    {
      break;
    }

    uv_temp = g_ctx.mesh.Parameterized_uv_new;
    

  }
  
  std::cout << "Llyod done." << std::endl;
  std::cout << "\tErrors:" << std::endl;

  for(double err: error_hist)
  {
    std::cout << "\t" << err << std::endl;
  }


  // option to show the remeshed plane
  g_ctx.show_original = false;
  resetDisplay(viewer);
  
  

}

////////////////////////////////////////////////////////////////////////////////
// parameterization
void
doParameterization (Context & g_ctx)
{
  Eigen::VectorXi boundary;
  igl::boundary_loop(g_ctx.mesh.F, boundary);

  // compute bundary points in uv plane
  igl::map_vertices_to_circle(g_ctx.mesh.V, boundary, g_ctx.mesh.bnd_uv);

  // Compute the initial solution for ARAP (harmonic parametrization)
  igl::harmonic(
    g_ctx.mesh.V,
    g_ctx.mesh.F,
    boundary,
    g_ctx.mesh.bnd_uv,
    1,
    g_ctx.mesh.Parameterized_uv
  );

  std::cout << "\tboundary: " << boundary.size() << std::endl;
  std::cout << "Parameterized_uv: "
    << g_ctx.mesh.Parameterized_uv.rows() << ","
    << g_ctx.mesh.Parameterized_uv.cols() << std::endl;
    
}

////////////////////////////////////////////////////////////////////////////////
void
doRemeshing (igl::opengl::glfw::Viewer& viewer, Context & g_ctx)
{
  // parameterization
  doParameterization(g_ctx);

  // delaunay triangulation
  // reference: https://github.com/libigl/libigl/blob/main/include/igl/copyleft/cgal/delaunay_triangulation.cpp
  std::vector< std::pair<Point,unsigned int> > p;
  for (size_t i = 0; i < g_ctx.mesh.Parameterized_uv.rows(); i++)
  {
    p.push_back(std::make_pair(
      Point(
        g_ctx.mesh.Parameterized_uv(i, 0),
        g_ctx.mesh.Parameterized_uv(i, 1)),
      i // info to be stored
    ));
  }

  Delaunay dt;
  dt.insert(p.begin(), p.end());

  g_ctx.mesh.Fuv.resize(g_ctx.mesh.F.rows(), 3);
  Delaunay::Finite_faces_iterator f_itr = dt.finite_faces_begin();
  size_t ind = 0;
  for (; f_itr != dt.finite_faces_end(); f_itr++)
  {
    g_ctx.mesh.Fuv(ind,0) = f_itr->vertex(0)->info();
    g_ctx.mesh.Fuv(ind,1) = f_itr->vertex(1)->info();
    g_ctx.mesh.Fuv(ind,2) = f_itr->vertex(2)->info();


    ind++;
  }

  // std::cout << g_ctx.mesh.Parameterized_uv << std::endl;
  // std::cout << g_ctx.mesh.Fuv << std::endl;
  // std::cout << "original: " << std::endl;
  // std::cout << g_ctx.mesh.V << std::endl;
  // std::cout << g_ctx.mesh.F << std::endl;

  {
    if (g_ctx.show_para)
    {
      add_points(
        viewer,
        g_ctx.mesh.bnd_uv,
        Eigen::RowVector3d(0, 1, 0)
      );
      viewer.data().point_size = 10;
    }
  }

  // iterative lloyd relaxation algorithm
  iterativeLloydRelaxation(viewer, g_ctx, dt);

  // compute the barcentric coordinates
  getBarcentricCoordinates(g_ctx, dt);

  std::cout << "Remeshing Finished!" << std::endl;
  
}

////////////////////////////////////////////////////////////////////////////////
void
resetDisplay (igl::opengl::glfw::Viewer& viewer)
{
  viewer.data().clear();

  // reset all options
  viewer.data().show_texture = g_ctx.show_texture;


  // if option is true: show the parameterized uv plane
  if (g_ctx.show_para)
  {
    if (g_ctx.mesh.Parameterized_uv.size() <= 0)
    {
      doParameterization(g_ctx);
    }

    Eigen::MatrixXd uv;
    if (g_ctx.show_original)
    {
      uv = g_ctx.mesh.Parameterized_uv;

      // show centroids if done remeshing
      if (g_ctx.mesh.Parameterized_uv_new.rows()>0)
      {
        add_points(
          viewer,
          g_ctx.mesh.Parameterized_uv_new,
          Eigen::RowVector3d(1, 0, 0)
        );
      }
    }
    else
    {
      uv = g_ctx.mesh.Parameterized_uv_new;
    }

    
    viewer.data().set_mesh(uv,g_ctx.mesh.F);
    viewer.core().align_camera_center(uv,g_ctx.mesh.F);
  }
  else
  { 
    // if optioned to show remeshed model, check if done remeshing 
    if (!g_ctx.show_original)
    {
      // if done for remeshing, show the result
      if (g_ctx.mesh.V_remeshed.size() > 0)
      {
        viewer.data().set_mesh(g_ctx.mesh.V_remeshed,g_ctx.mesh.F);
        viewer.core().align_camera_center(g_ctx.mesh.V_remeshed,g_ctx.mesh.F);
        return;
      }
    }
    
    viewer.data().set_mesh(g_ctx.mesh.V,g_ctx.mesh.F);
    viewer.core().align_camera_center(g_ctx.mesh.V,g_ctx.mesh.F);
  }
}

////////////////////////////////////////////////////////////////////////////////
bool
keyDown (igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{

    std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;
    if (key=='q' || key=='Q')
    {
        exit(0);
    }
    else if (key=='h' || key=='H')
    {
        printHelp();
    }
    else if (key=='r' || key=='R')
    {
        doRemeshing(viewer, g_ctx);
    }
    else if (key=='p' || key=='P')
    {
      g_ctx.show_para = !g_ctx.show_para;
      resetDisplay(viewer);
    }
    else if (key=='s' || key=='S')
    {
      g_ctx.show_texture = !g_ctx.show_texture;
      resetDisplay(viewer);
    }
    else if (key=='b' || key=='B')
    {
      if (g_ctx.mesh.Parameterized_uv_new.rows()>0)
      {
        g_ctx.show_original = !g_ctx.show_original;
        std::cout << (g_ctx.show_original?"Original":"Remeshed") << std::endl;
        resetDisplay(viewer);
      }
      
    }
    else if (key=='d' || key=='D')
    {
      if (g_ctx.mesh.Parameterized_uv.rows() > 0)
      {
        g_cur_ind++;
        if (g_cur_ind>=g_ctx.mesh.Parameterized_uv.rows())
        {
          g_cur_ind = 0;
        }
        std::cout << "VERTEX: " << g_cur_ind+1 << "\t (" 
                  << g_ctx.mesh.Parameterized_uv.row(g_cur_ind) 
                  << ")" << std::endl;
        resetDisplay(viewer);
        add_points(
          viewer,
          g_ctx.mesh.Parameterized_uv.row(g_cur_ind),
          Eigen::RowVector3d(0, 1, 0)
        );
      }
      
    }

    viewer.data().compute_normals();
    return false;
}

////////////////////////////////////////////////////////////////////////////////
void
printHelp ()
{
  std::cout
    << "#####################################################################"
    << std::endl;
  std::cout << "Help manu:" << std::endl;
  std::cout << "  Key Inputs:\n";
  std::cout << "\tq, Q:  Quit.\n";
  std::cout << "\tr, R:  Process Remeshing.\n";
  std::cout << "\ts, S:  Switch for Showing Texture.\n";
  std::cout << "\tp, P:  Switch for Displaying Mesh or Parameterized UV plane.\n";
  std::cout << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
int main (int argc, char const *argv[])
{
  std::string data_path = "data/";
  std::string mesh_file_name;
  mesh_file_name = "camelhead.obj";
  // mesh_file_name = "terrain.obj";
  // mesh_file_name = "plane_simple.obj";

  Mesh m;
  bool file_found = igl::readOBJ(data_path+mesh_file_name, m.V, m.F);

  if (file_found)
  {
    g_ctx.mesh = m;
    std::cout << "Mesh file loaded at: "<< data_path+mesh_file_name<< std::endl;
    std::cout << "Vertices: " << m.V.rows() << std::endl;
    std::cout << "Faces: " << m.F.rows() << std::endl;

    // init viewer
    igl::opengl::glfw::Viewer viewer;
    viewer.callback_key_down = &keyDown;

    // display the mesh
    viewer.data().set_mesh(m.V, m.F);
    viewer.launch();


  }
  else
  {
    std::cout << "Mesh file " << data_path+mesh_file_name
              << " NOT found!"<< std::endl;
  }

  return 0;
}
