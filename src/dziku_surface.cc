#include <config.h>

#include <type_traits>
#include <vector>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/curvedgrid/geometries/implicitsurface.hh>
#include <dune/geometry/quadraturerules/compositequadraturerule.hh>
#include <dune/functions/common/differentiablefunctionfromcallables.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvedgrid/grid.hh>
#include <dune/gmsh4/gmsh4reader.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/vtk/vtkwriter.hh>
#include <dune/vtk/datacollectors/lagrangedatacollector.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <math.h>


  using namespace Dune;
  
const int quad_order = 14;
const int num_levels = 1;
double f(const FieldVector<double,3>& x)
{
    auto numerator = pow(x[1], 2) + pow(x[2], 2) - (x[0] - pow(x[2], 2)) * (x[0] * (2 * x[0] - 1) + 2 * pow(x[1], 2) + pow(x[2], 2) - 4 * x[0] * pow(x[2], 2) + 2 * pow(x[2], 4));
    auto denominator = pow(x[1], 2) + pow(x[2], 2) + (x[0] - pow(x[2], 2)) * (x[0] + (4 * x[0] - 5) * pow(x[2], 2) - 4 * pow(x[2], 4));
    return numerator / pow(denominator, 2);
      
  };


template <class Grid>
typename Grid::ctype L_error (const Grid& grid)
{
  std::cout << "  L_error..." << std::endl;
   using QuadProvider = Dune::QuadratureRules<typename Grid::ctype, 2>;
  typename Grid::ctype result=0;
 for (const auto& element : elements(grid.leafGridView()))
  {
 auto geo=element.geometry();
const auto& quadRule = QuadProvider::rule(element.type(), quad_order);
    for (const auto& quadPoint : quadRule) 
    {

result+=f(geo.global(quadPoint.position()))*geo.integrationElement(quadPoint.position())*quadPoint.weight();

}
}
return std::abs(4*M_PI-result);
}

template <class Grid>
typename Grid::ctype edge_length (const Grid& grid)
{
  typename Grid::ctype h = 0;
  for (const auto& e : edges(grid.hostGrid().leafGridView()))
    h = std::max(h, e.geometry().volume());

  return h;
}
    
int main(int argc, char** argv)
{   using namespace Dune;


  MPIHelper::instance(argc, argv);

  // Construct a reference grid
//vogitMesh_N=3144_d=0.5_c=0.375.msh
auto refGrid = Gmsh4Reader<FoamGrid<2,3>>::createGridFromFile(DUNE_GRID_PATH "dziku_surface_N=8088.msh");

    // Define the levelset function
  using Signature = double(FieldVector<double,3>);
  auto psi = Functions::makeDifferentiableFunctionFromCallables(
    Functions::SignatureTag<Signature>{},
    [](auto const& x) -> double {
      double d=0.8,c=-0.934;
      return pow((x[0] - pow(x[2], 2)), 2) + pow(x[1], 2) + pow(x[2], 2) - 1;
    },
    [](auto const& x) -> FieldVector<double,3> {
      return {2*(x[0]-pow(x[2],2)), 2*x[1], 2*(-2*x[0]*x[2]+2*pow(x[2],3)+x[2])};
    }
  );

// Construct a projection by an iteration scheme
auto surfaceGF1 = SimpleImplicitSurfaceProjection{psi};
  using std::log;
  using std::abs;
  CurvedGrid grid{*refGrid, surfaceGF1,4};
  using Grid = decltype(grid);
  std::vector<typename Grid::ctype> L_errors, edge_lengths;
  for (int i = 0; i < num_levels; ++i) {
    if (i > 0)
      grid.globalRefine(1);

    std::cout << "level = " << i << std::endl;
    L_errors.push_back(L_error(grid));
    edge_lengths.push_back( edge_length(grid));}
   // longest edge in the host grid

      auto u = Dune::Functions::makeAnalyticGridViewFunction(f, grid.leafGridView());
  Vtk::LagrangeDataCollector dataCollector{grid.leafGridView(), 3};
  VtkUnstructuredGridWriter writer2{dataCollector};
  writer2.write("dziku_surface.vtu");
    writer2.addPointData(u, Dune::Vtk::FieldInfo{"K\_Gauss", 1, Dune::Vtk::RangeTypes::SCALAR});
    writer2.write("dziku_surface_scalar.vtu");
  
    std::vector<typename Grid::ctype>  eocL(num_levels, 0) ;
  for (int i = 1; i < num_levels; ++i) {
    auto ratio = edge_lengths[i]/edge_lengths[i-1];
    eocL[i] = log(L_errors[i]/L_errors[i-1]) / log(ratio);

  }

auto print_line = [](auto i, const auto& data) {
    std::cout << i;
    for (const auto& d : data) {
      std::cout << '\t' << "| ";
      std::cout << d;
    }
    std::cout << std::endl;
  };

  auto print_break = [] {
    std::cout.width(8 + 9*16);
    std::cout.fill('-');
    std::cout << '-' << std::endl;
    std::cout.fill(' ');
  };

  std::cout.setf(std::ios::scientific);
  print_line("level", std::vector<std::string>{"h       ",
  "Error"  "(eoc)   "});
  print_break();

  for (int i = 0; i < num_levels; ++i)
    print_line(i,std::vector<typename Grid::ctype>{
      edge_lengths[i],
      L_errors[i], eocL[i]});
}
