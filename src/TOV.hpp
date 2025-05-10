#ifndef TOV_H
#define TOV_H

#include <vector>



struct data {
    double mass;
    double radius;
};

struct vec4 {
    double x, y, z, w;
};

int find_ind(const std::vector<double> &arr, double val);

double ene_interp(const std::vector<double> &pre_arr,
                  const std::vector<double> &ene_arr, double pressure);

double pre_interp(const std::vector<double> &pre_arr,
                  const std::vector<double> &ene_arr, double energy);

double en_dens(const std::vector<double> &parr, const std::vector<double> &earr,double P);

data beta_and_H( double r,  double p,  double H ,  double m,  double beta, const std::vector<double> &parr,const std::vector<double> &earr);

inline double love_number( double C,  double y);

inline double Tov_eqn( double P, double r,  double m, const std::vector<double> &dens,
               const std::vector<double> &press,  double min_pressure);

inline double mass_eqn(double r, double ene);

std::vector<data> ToV(const std::vector<double> &cen_dens,
                      const std::vector<double> &pressure,
                      const std::vector<double> &density, double exit_Pre);

std::vector<vec4> ToV_love(const std::vector<double> &cen_dens,
                      const std::vector<double> &pressure,
                      const std::vector<double> &density, double exit_Pre);

std::vector<double> logspace( double start_log10, double end_log10,
                              int num);

std::vector<data> read_csv(const std::string &filename, int density_col,
                           int pressure_col);

#endif