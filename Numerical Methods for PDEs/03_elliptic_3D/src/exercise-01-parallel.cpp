#include "DiffusionReaction.hpp"

static constexpr unsigned int dim = DiffusionReaction::dim;

// Main function.
int
main(int /*argc*/, char * /*argv*/[])
{

  const std::vector<unsigned int> N_el_values = {5, 10, 20, 40};
  const unsigned int              r           = 2;

  const auto mu = [](const Point<dim> & p) {
    if (p[0] <= 0.5){
      return 100.;
    } else {
      return 1.;
    }
  };

  const auto sigma = [](const Point<dim> & /*p*/) { return 1.0; };

  const auto f = [](const Point<dim> & /*p*/) {return 1.;};

  for (const auto &N_el : N_el_values)
    {
      const std::string mesh_file_name =
        "../mesh/mesh-cube-" + std::to_string(N_el) + ".msh";

        DiffusionReaction problem(mesh_file_name, r, mu, sigma, f);

      problem.setup();
      problem.assemble();
      problem.solve();
      problem.output();

    }

  return 0;
}