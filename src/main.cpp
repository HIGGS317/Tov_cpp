#include "TOV.hpp"
#include "constants.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char *argv[]) {

  try {
    int a, b,num;
    char love;

    std::ifstream file(argv[1]);
    std::string output_file;

    // if (!file.is_open()) {
    //   std::cerr << "Cannot open file: " << argv[1] << "\n";
    //   return 1;
    // }

    std::cout << "Enter the coloumn number for pressure value\n";
    std::cin >> a;

    std::cout << "Enter the coloumn number for density value\n";
    std::cin >> b;

    std::cout<<"Do you want to calculate love number (Y/N) \n";
    std::cin>>love;

    std::cout << "Enter the number MR points \n";
    std::cin >> num;

    std::cout << "Enter the name of output file\n";
    std::cin >> output_file;

    auto table = read_csv(argv[1], b, a);

    std::vector<double> pressure, density;
    for (const auto &row : table) {
      auto pressure1 = row.radius * Mev_fm3_to_GU;
      pressure.push_back(pressure1);
      // std::cout << "pressure = " << pressure1 << "\n";
    }

    for (const auto &row : table) {
      auto density1 = row.mass * Mev_fm3_to_GU;
      density.push_back(density1);
      // std::cout << "density = " << density1 << "\n";
    }

    double max_density = *std::max_element(density.begin(), density.end());
    // std::cout << "Maximum density = " << max_density << "\n";

    // Step 2: prepare limits
    double upper = std::log10(max_density * 0.999);
    double lower = std::log10(2.072392843084521e-10);

    std::vector<double> central_energy = logspace(upper, lower, num);

    // std::cout << "Minimum Central_energy" << lower << "\n";
    // std::cout << "Maximum Central_energy" << upper << "\n";

    double min_pressure = *std::min_element(pressure.begin(), pressure.end());

    double minimum_pressure = min_pressure * 1.1;
    // std::cout << "minimum pressure = " << minimum_pressure << "\n";
if (love == 'Y') {
  std::vector<vec4> MR =
      ToV_love(central_energy, pressure, density, minimum_pressure);
  for (size_t i = 0; i < 50; ++i) {
    std::cout << "Mass=" << MR[i].x << " radius=" << MR[i].y << "Compactness = "<<MR[i].z << "Love Number = "<<MR[i].w <<'\n';
  }

  std::ofstream out_file(output_file);
  if (!out_file) {
    std::cerr << "Failed to open output file.\n";
    return 1;
  }

  // Write header
  out_file << "M,R,Compactness,Love Number\n";

  // Write data
  for (size_t i = 0; i < central_energy.size(); ++i) {
    out_file << "Mass=" << MR[i].x << " radius=" << MR[i].y << "Compactness = "<<MR[i].z << "Love Number = "<<MR[i].w <<'\n';
  }
}else {
  std::vector<data> MR =
      ToV(central_energy, pressure, density, minimum_pressure);

    // Print first 10 results
    for (size_t i = 0; i < 50; ++i) {
      std::cout << "Mass=" << MR[i].mass << " radius=" << MR[i].radius << '\n';
    }

    std::ofstream out_file(output_file);
    if (!out_file) {
      std::cerr << "Failed to open output file.\n";
      return 1;
    }

    // Write header
    out_file << "M,R\n";

    // Write data
    for (size_t i = 0; i < central_energy.size(); ++i) {
      out_file << MR[i].mass << ',' << MR[i].radius << '\n';
    }
    }
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}