#include <string>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Voronoi_diagram_2.h>



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

struct Mesh
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd bnd_uv; // boundary point in uv plane
  Eigen::MatrixXi Fuv; // Faces in uv plane
  Eigen::MatrixXd Parameterized_uv; // points in uv plane
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
iterativeLloydRelaxation (Context & g_ctx, Delaunay & dt);






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
iterativeLloydRelaxation (Context & g_ctx, Delaunay & dt)
{
  Delaunay::Finite_vertices_iterator v_itr = dt.finite_vertices_begin();

  std::cout << "vertex: " << v_itr->info() << std::endl;
  std::cout << "Total points " << dt.number_of_vertices() << std::endl;
  std::cout << "There are " << g_ctx.mesh.bnd_uv.rows() << " boundary points."
            << std::endl;

  VD vd(dt);
  assert(vd.is_valid());

  bool is_boundary;
  for (; v_itr != dt.finite_vertices_end(); v_itr++)
  {
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
      std::cout << "x=" << v_itr->point().x() << ", y=" << v_itr->point().y()
                << std::endl;

                
      auto cell = vd.dual(v_itr);
      auto ec_begin = cell->ccb();

      auto ec = ec_begin;
      do
      {
        std::cout << ec->source()->point().x() << std::endl;
        std::cout << ec->source()->point().y() << std::endl;
      } while (++ec != ec_begin);
    }

  }

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

  std::cout << "boundary: " << boundary.size() << std::endl;
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
      i
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
  {
    if (g_ctx.show_para)
    {
      // add_points(
      //   viewer,
      //   g_ctx.mesh.Parameterized_uv,
      //   Eigen::RowVector3d(1, 0, 0)
      // );
      add_points(
        viewer,
        g_ctx.mesh.bnd_uv,
        Eigen::RowVector3d(0, 1, 0)
      );
      viewer.data().point_size = 10;
    }
  }

  // iterative lloyd relaxation algorithm
  iterativeLloydRelaxation(g_ctx, dt);
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
    viewer.data().set_mesh(g_ctx.mesh.Parameterized_uv,g_ctx.mesh.F);
    viewer.core().align_camera_center(g_ctx.mesh.Parameterized_uv,g_ctx.mesh.F);
  }
  else{
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
    else if (key=='d' || key=='D')
    {
      g_cur_ind++;
      if (g_cur_ind>=g_ctx.mesh.V.rows())
      {
        g_cur_ind = 0;
      }
      std::cout << "FACE: " << g_cur_ind+1 << std::endl;
      resetDisplay(viewer);
      add_points(
        viewer,
        g_ctx.mesh.Parameterized_uv.row(g_cur_ind),
        Eigen::RowVector3d(0, 1, 0)
      );
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
  std::string mesh_file_name = "plane_simple.obj";

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
