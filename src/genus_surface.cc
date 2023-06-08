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

auto x0=x[0], x1=x[1], x2=x[2];

return (4*1*(-900*(pow(x[0], 2) + pow(x[1], 2))*pow(x[2], 2) + 45*1*(pow(x[0], 2) + pow(x[1], 2))*pow((-3*pow(x[0], 2)*x[1] + pow(x[1], 3)),2)*pow(x[2], 6) - 6*1*x[1]*(-3*pow(x[0], 2) + pow(x[1], 2))*pow(x[2], 2)*(159*(pow(x[0], 2) + pow(x[1], 2)) - 460*pow(x[2], 2))+\
                                             15*1*x[1]*(-3*pow(x[0], 2) + pow(x[1], 2))*(9*(pow(x[0], 2) + pow(x[1], 2)) - 40*pow(x[2], 2)) + 15*1*x[1]*(-3*pow(x[0], 2) + pow(x[1], 2))*pow(x[2], 4)*(3*pow(x[0], 6) - 9*pow(x[0], 4)*pow(x[1], 2) + 21*pow(x[0], 2)*pow(x[1], 4) + pow(x[1], 6) + 27*(pow(x[0], 2) + pow(x[1], 2))*pow(x[2], 4)) + \
                                             15*1*(3*(pow(x[0], 6) + 21*pow(x[0], 4)*pow(x[1], 2) - 9*pow(x[0], 2)*pow(x[1], 4) + 3*pow(x[1], 6)) + 20*pow((pow(x[0], 2) + pow(x[1], 2)),2)*pow(x[2], 2) + 336*(pow(x[0], 2) + pow(x[1], 2))*pow(x[2], 4)) + 9*1*x[1]*(-3*pow(x[0], 2) + pow(x[1], 2))*(pow(x[0], 6) + 21*pow(x[0], 4)*pow(x[1], 2) - 9*pow(x[0], 2)*pow(x[1], 4) + 3*pow(x[1], 6) + 212*(pow(x[0], 2) + pow(x[1], 2))*pow(x[2], 4) - 456*pow(x[2], 6)) + \
                                             1*(-20*pow(x[0], 8) + 163*pow(x[0], 6)*pow(x[1], 2) - 39*pow(x[0], 4)*pow(x[1], 4) - 215*pow(x[0], 2)*pow(x[1], 6) + 7*pow(x[1], 8) - 3*(171*pow(x[0], 6) + 2151*pow(x[0], 4)*pow(x[1], 2) - 579*pow(x[0], 2)*pow(x[1], 4) + 353*pow(x[1], 6))*pow(x[2], 2) - 1080*pow((pow(x[0], 2) + pow(x[1], 2)),2)*pow(x[2], 4) - 10296*(pow(x[0], 2) + pow(x[1], 2))*pow(x[2], 6)) + \
                                             3*1*pow(x[2], 2)*(3*(pow(x[0], 2) + pow(x[1], 2))*(12*pow(x[0], 6) + 27*pow(x[0], 4)*pow(x[1], 2) + 42*pow(x[0], 2)*pow(x[1], 4) + 11*pow(x[1], 6)) + (345*pow(x[0], 6) + 3213*pow(x[0], 4)*pow(x[1], 2) - 417*pow(x[0], 2)*pow(x[1], 4) + 587*pow(x[1], 6))*pow(x[2], 2) + 324*pow((pow(x[0], 2) + pow(x[1], 2)),2)*pow(x[2], 4) + 3024*(pow(x[0], 2) + pow(x[1], 2))*pow(x[2], 6)) - \
                                             2*1*x[1]*(-3*pow(x[0], 2) + pow(x[1], 2))*(2*pow((pow(x[0], 2) + pow(x[1], 2)),4) + 3*(9*pow(x[0], 6) + 9*pow(x[0], 4)*pow(x[1], 2) + 39*pow(x[0], 2)*pow(x[1], 4) + 7*pow(x[1], 6))*pow(x[2], 2) + 747*(pow(x[0], 2) + pow(x[1], 2))*pow(x[2], 6) - 972*pow(x[2], 8)) + 3*1*pow(x[2], 2)*(-4*pow((pow(x[0], 2) + pow(x[1], 2)),2)*(3*pow(x[0], 6) - 9*pow(x[0], 4)*pow(x[1], 2) + \
                                             21*pow(x[0], 2)*pow(x[1], 4) + pow(x[1], 6)) - 21*pow(x[1], 2)*pow((-3*pow(x[0], 2) + pow(x[1], 2)),2)*(pow(x[0], 2) + pow(x[1], 2))*pow(x[2], 2) - 9*(21*pow(x[0], 6) + 153*pow(x[0], 4)*pow(x[1], 2) + 3*pow(x[0], 2)*pow(x[1], 4) + 31*pow(x[1], 6))*pow(x[2], 4) - 972*(pow(x[0], 2) + pow(x[1], 2))*pow(x[2], 8))))/pow((100*pow(x[2],2) - 12*1*x[1]*(-3*pow(x[0],2) + pow(x[1],2))*pow(x[2],2)*(pow(x[0],2) + pow(x[1],2) + 6*pow(x[2],2)) + \
                                             4*1*x[1]*(-3*pow(x[0],2) + pow(x[1],2))*(3*(pow(x[0],2) + pow(x[1],2)) + 10*pow(x[2],2)) + 1*pow(x[2],2)*(4*pow((-3*pow(x[0],2)*x[1] + pow(x[1],3)),2) + 9*pow((pow(x[0],2) + pow(x[1],2)),2)*pow(x[2],2)) + 9*1*(pow((pow(x[0],2) + pow(x[1],2)),2) - 40*pow(x[2],4)) + 2*1*(2*pow((pow(x[0],2) + pow(x[1],2)),3) - 9*pow((pow(x[0],2) + pow(x[1],2)),2)*pow(x[2],2) + 162*pow(x[2],6))),2)
;
          
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
      //std::cout<< "p=" <<std:: setprecision (15)<<geo.global(quadPoint.position())<<std::endl;
result+=f(geo.global(quadPoint.position()))*geo.integrationElement(quadPoint.position())*quadPoint.weight();

}
}
return std::abs(-4*M_PI-result);
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
auto refGrid = Gmsh4Reader<FoamGrid<2,3>>::createGridFromFile(DUNE_GRID_PATH "genus_N=15632.msh");

  // Define the levelset function
  using Signature = double(FieldVector<double,3>);
  auto psi = Functions::makeDifferentiableFunctionFromCallables(
    Functions::SignatureTag<Signature>{},
    [](auto const& x) -> double {
      auto x2 = x[0]*x[0], y2 = x[1]*x[1], z2 = x[2]*x[2];
      return 2*x[1]*(y2 - 3*x2)*(1 - z2)+power(x2 + y2, 2)-(9*z2 - 1)*(1 - z2);
    },
    [](auto const& x) -> FieldVector<double,3> {
      auto x2 = x[0]*x[0], y2 = x[1]*x[1], z2 = x[2]*x[2];
      return {
        4*x[0]*(x2 + y2 + 3*x[1]*(z2 - 1)) ,
        4*x[1]*(x2 + y2) + 4*y2*(1 - z2) + 2*(3*x2 - y2)*(z2 - 1) ,
        4*x[2]*(x[1]*(3*x2 - y2) + 9*z2 - 5)
      };
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
  writer2.write("genus_surface.vtu");
    writer2.addPointData(u, Dune::Vtk::FieldInfo{"K\_Gauss", 1, Dune::Vtk::RangeTypes::SCALAR});
    writer2.write("genus_surface_scalar.vtu");
  
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
