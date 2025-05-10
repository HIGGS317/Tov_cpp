#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "TOV.hpp"
#include "constants.hpp"




int find_ind(const std::vector<double> &arr, double val) {
  for (size_t i = 0; i < arr.size(); ++i) {
    if (val > arr[i]) {
      continue;
    } else {
      return static_cast<int>(i);
    }
  }
  return static_cast<int>(arr.size());
}

double ene_interp(const std::vector<double> &pre_arr,
                  const std::vector<double> &ene_arr, double pressure) {
  if (pressure < *std::min_element(pre_arr.begin(), pre_arr.end()) or
      pressure > *std::max_element(pre_arr.begin(), pre_arr.end())) {

    throw std::runtime_error("pressure out of range");
  } else {
    const int ind = find_ind(pre_arr, pressure);
    double left_p = pre_arr[ind - 1];
    double right_p = pre_arr[ind];
    double left_e = ene_arr[ind - 1];
    double right_e = ene_arr[ind];
    double ene_val =
        (pressure - left_p) * (right_e - left_e) / (right_p - left_p) + left_e;
    return ene_val;
  }
}

double pre_interp(const std::vector<double> &pre_arr,
                  const std::vector<double> &ene_arr, double energy) {
  if (energy < *std::min_element(ene_arr.begin(), ene_arr.end()) or
      energy > *std::max_element(ene_arr.begin(), ene_arr.end())) {

    throw std::runtime_error("pressure out of range");
  } else {
    int ind = find_ind(ene_arr, energy);
    double left_p = pre_arr[ind - 1];
    double right_p = pre_arr[ind];
    double left_e = ene_arr[ind - 1];
    double right_e = ene_arr[ind];
    double pre_val =
        (energy - left_e) * (right_p - left_p) / (right_e - left_e) + left_p;
    return pre_val;
  }
}



double en_dens(const std::vector<double> &parr, const std::vector<double> &earr,double P) {
  double e;
  if (P < *std::min_element(parr.begin(), parr.end()) or P > *std::max_element(parr.begin(), parr.end())) {
    e = 0;
  }
  else{
    e = ene_interp(parr, earr, P);

  }
  return e;
}

data beta_and_H(const double r, const double p, const double H , const double m, const double beta, const std::vector<double> &parr,const std::vector<double> &earr) {
  const double dp = p*0.005;
  const double el_3 = en_dens(parr, earr, p - 3 * dp);
  const double el_2 = en_dens(parr, earr, p - 2 * dp);
  const double el_1 = en_dens(parr, earr, p - 1 * dp);
  const double er_3 = en_dens(parr, earr, p + 3 * dp);
  const double er_2 = en_dens(parr, earr, p + 2 * dp);
  const double er_1 = en_dens(parr, earr, p + 1 * dp);

  const double de_dp =  (-1.0 / 60 * el_3 + 3.0 / 20 * el_2 - 3.0 / 4 * el_1 + 3.0 / 4 * er_1 - 3.0 / 20 * er_2 + 1.0 / 60 * er_3) / dp;

  const double e = en_dens(parr, earr, p);

  double dbeta_dr = 2 * std::pow((1 - 2 * m / r ) ,-1) *H * std::pow(-2 * pi * G / std::pow(c,2) * (5 * e + 9 * p / std::pow(c,2) + de_dp * std::pow(c,2) * (e + p / std::pow(c,2))) + 3 / std::pow(r,2) + 2 * std::pow((1 - 2 * m / r ) , -1) * ( m / std::pow(r,2) + G / std::pow(c,4) * 4 * pi * r * p), 2) + 2 * std::pow((1 - 2 * m / r ) , -1) *beta / r * (-1 + m / r + 2 * pi * std::pow(r,2) * G / std::pow(c,2) * (e - p / std::pow(c, 2)));

  const double dHdr = beta;

  return {dbeta_dr, dHdr};
}

inline double love_number(const double C, const double y) {
  const double k2 = 8.0 / 5 * std::pow(C , 5) * std::pow((1 - 2 * C) ,2) * (2 + 2 * C * (y - 1) - y) * (2 * C * (6 - 3 * y + 3 * C * (5 * y - 8)) + 4 * std::pow(C ,3) * (13 - 11 * y + C * (3 * y - 2) + 2 * std::pow(C , 2) * (1 + y)) + 3 * std::pow((1 - 2 * C), 2) * (2 - y + 2 * C * (y - 1)) * std::pow(((std::log(1 - 2 * C))), -1));
  return k2;

}

inline double Tov_eqn(const double P, const double r, const double m, const std::vector<double> &dens,
               const std::vector<double> &press, const double min_pressure) {
  if (P < min_pressure) {

    return 0.0;

  } else {

    const double eden = ene_interp(press, dens, P);
    // return -(G * ((P / std::pow(c, 2)) + eden)
    //          *(m + 4 * pi * std::pow(r, 3) * P / std::pow(c, 2))) /
    //        (r * (r - 2 * G * m / std::pow(c, 2)));

    // return (-2*G*eden*m/(std::pow(r,2)))/(1-(2*G*m/r));

    return (-2*G*eden*m/std::pow(r,2))*(1+P/(eden*std::pow(c,2)))*(1+(4*pi*std::pow(r,3)*P)/(m*std::pow(c,2)));


  }
}

inline double mass_eqn(double r, double ene) { return 4 * pi * std::pow(r, 2) * ene; }



std::vector<data> ToV(const std::vector<double> &cen_dens,
                      const std::vector<double> &pressure,
                      const std::vector<double> &density, double exit_Pre) {

  double P_exit = exit_Pre;
  std::vector<data> MR;
  for (size_t i = 0; i < cen_dens.size(); i++) {
    data info{0,0};

    const std::vector<double> &press = pressure;
    const std::vector<double>& dens = density;
    const double d = cen_dens[i];
    const double P0 = pre_interp(press, dens, d);
    double r = 10;
    double P = P0;
    double m = mass_eqn(r, d);
    const double min_pressure = *std::min_element(press.begin(), press.end());
    std::cout << "Iteration Number = " << i << "\n";

    while (P > P_exit) {
      constexpr double h = 1;
      // std::cout << "Pressure in the loop = " << P << "\n";
      // std::cout << "Exit pressure = " << P_exit << "\n";
      double k1_m = mass_eqn(r, ene_interp(press, dens, P));
      double k2_m = mass_eqn(r + h / 2, ene_interp(press, dens, P));
      double k3_m = mass_eqn(r + h / 2, ene_interp(press, dens, P));
      double k4_m = mass_eqn(r + h, ene_interp(press, dens, P));

      double k1_p = Tov_eqn(P, r, m, dens, press, min_pressure);
      double k2_p = Tov_eqn(P + k1_p * h / 2, r + h / 2, m + k1_m * h / 2, dens,
                            press, min_pressure);
      double k3_p = Tov_eqn(P + k2_p * h / 2, r + h / 2, m + k2_m * h / 2, dens,
                            press, min_pressure);
      double k4_p =
          Tov_eqn(P + k3_p * h, r + h, m + k3_m * h, dens, press, min_pressure);

      P += h * (k1_p + 2 * k2_p + 2 * k3_p + k4_p) / 6;
      m += h * (k1_m + 2 * k2_m + 2 * k3_m + k4_m) / 6;
      r += h;
    }

    info.mass = m / M01;
    info.radius = r / 1000;

    MR.push_back(info);
  }
  // std::vector<data> MR;
  // MR.push_back(info);

  return MR;
}

std::vector<vec4> ToV_love(const std::vector<double> &cen_dens,
                      const std::vector<double> &pressure,
                      const std::vector<double> &density, double exit_Pre) {
  double P_exit = exit_Pre;
  std::vector<vec4> MR;
  for (size_t i = 0; i < cen_dens.size(); i++) {
    vec4 info{0,0,0,0};

    const std::vector<double> &press = pressure;
    const std::vector<double>& dens = density;
    const double d = cen_dens[i];
    const double P0 = pre_interp(press, dens, d);
    double r = 10;
    double P = P0;
    double m = mass_eqn(r, d);
    const double min_pressure = *std::min_element(press.begin(), press.end());
    constexpr double a0 = 1;
    const double  H0 = a0*std::pow(r,2);
    const double beta0 = 2*a0*r;
    double beta = beta0;
    double H = H0;
    std::cout << "Iteration Number = " << i << "\n";

    while (P > P_exit) {
      constexpr double h = 1;
      // std::cout << "Pressure in the loop = " << P << "\n";
      // std::cout << "Exit pressure = " << P_exit << "\n";
      double k1_m = mass_eqn(r, ene_interp(press, dens, P));
      double k2_m = mass_eqn(r + h / 2, ene_interp(press, dens, P));
      double k3_m = mass_eqn(r + h / 2, ene_interp(press, dens, P));
      double k4_m = mass_eqn(r + h, ene_interp(press, dens, P));

      double k1_p = Tov_eqn(P, r, m, dens, press, min_pressure);
      double k2_p = Tov_eqn(P + k1_p * h / 2, r + h / 2, m + k1_m * h / 2, dens,
                            press, min_pressure);
      double k3_p = Tov_eqn(P + k2_p * h / 2, r + h / 2, m + k2_m * h / 2, dens,
                            press, min_pressure);
      double k4_p =
          Tov_eqn(P + k3_p * h, r + h, m + k3_m * h, dens, press, min_pressure);


      auto  k1_dHdr = beta_and_H(r, P, H, m, beta, press, dens);
      auto  k2_dHdr = beta_and_H(r + 0.5 * h, P + 0.5 * h * k1_p, H + 0.5 * h * k1_dHdr.radius, m + 0.5*h*k1_m, beta + 0.5*h*k1_dHdr.mass, press, dens);
      auto  k3_dHdr = beta_and_H(r + 0.5 * h, P + 0.5 * h * k2_p, H + 0.5 * h * k2_dHdr.radius, m + 0.5*h*k2_m, beta + 0.5*h*k2_dHdr.mass, press, dens);
      auto  k4_dHdr = beta_and_H(r + h, P + h * k3_p, H + h * k3_dHdr.radius, m + h*k3_m, beta + h*k3_dHdr.mass, press, dens);

      beta = beta + (h / 6.0) * (k1_dHdr.mass + 2 * k2_dHdr.mass + 2 * k3_dHdr.mass + k4_dHdr.mass);
      H = H + (h / 6.0) * (k1_dHdr.radius + 2 * k2_dHdr.radius + 2 * k3_dHdr.radius + k4_dHdr.radius);

      P += h * (k1_p + 2 * k2_p + 2 * k3_p + k4_p) / 6;
      m += h * (k1_m + 2 * k2_m + 2 * k3_m + k4_m) / 6;
      r += h;
    }

    info.x = m / M01;
    info.y = r / 1000;
    double y = r*beta/H;
    double C = m/r;
    double k2 = love_number(C, y);
    info.z = k2;
    info.w = beta;
    MR.push_back(info);
  }
  return MR;
}

std::vector<double> logspace(const double start_log10, const double end_log10,
                             const int num) {
  std::vector<double> result;
  result.reserve(num);
  double step = (end_log10 - start_log10) / (num - 1);

  for (int i = 0; i < num; ++i) {
    result.emplace_back(std::pow(10, start_log10 + i * step));
  }
  return result;
}

std::vector<data> read_csv(const std::string &filename, int density_col,
                           int pressure_col) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + filename);
  }

  std::vector<data> table;
  std::string line;

  // Skip header
  std::getline(file, line);

  // Process data rows
  while (std::getline(file, line)) {
    std::istringstream ss(line);
    std::string cell;
    int col_index = 1; // 1-based column index
    double density = 0.0, pressure = 0.0;

    // Parse CSV by splitting on commas
    while (std::getline(ss, cell, ',')) {
      try {
        if (col_index == density_col) {
          density = std::stod(cell);
        } else if (col_index == pressure_col) {
          pressure = std::stod(cell);
        }
      } catch (const std::invalid_argument &) {
        std::cerr << "Error converting value: " << cell << " in line: " << line
                  << "\n";
        throw;
      } catch (const std::out_of_range &) {
        std::cerr << "Value out of range: " << cell << " in line: " << line
                  << "\n";
        throw;
      }
      ++col_index;
    }

    // Add the parsed density and pressure to the table
    if (col_index > std::max(density_col, pressure_col)) {
      table.push_back({density, pressure});
    }
  }
  return table;
}