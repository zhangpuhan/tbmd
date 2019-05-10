//
//  analysis.hpp
//  GQMD-Hubbard-liquid
//
//

#ifndef analysis_hpp
#define analysis_hpp

#include <stdio.h>
#include "util.hpp"
#include "potential.hpp"
#include "tbmd.hpp"

using namespace std;

class MD_param {
public:

    double kT;

    double atom_mass;
    double Langevin_damping;

    double filling_fraction;
    int N_atoms;

    double l_box;

    double dlt_t;
    int N_steps;
};


class MD_Data {
public:

    double M1_ee, M2_ee;
    double M1_ep, M2_ep;
    double M1_ek, M2_ek;

    double M1_pr, M2_pr;        // pressure
    double M1_vl, M2_vl;

    double M1_vol, M2_vol;

    Vec<double> g_r;
    int num_r;
    double dr;

    int numSteps;
    double dt;


    MD_Data(void) {
        reset();
    };

    ~MD_Data(void) {

    };

    void reset(void) {
        M1_ee = M2_ee = M1_ep = M2_ep = M1_ek = M2_ek = 0;
        M1_pr = M2_pr = M1_vl = M2_vl = 0;
        M1_vol = M2_vol = 0;
    };

    void init_g_r(int size, double dlt_r) {

        g_r.resize(size);
        num_r = size;
        dr = dlt_r;

        for(int i=0; i<size; i++) g_r[i] = 0;
    };

    void print_g_r(string const filename) {

        std::ofstream fs;

        fs.open(filename.c_str(), ios::out);
        fs.precision(10);
        for(unsigned m=0; m<g_r.size(); m++) fs << (m + 0.5) * dr << '\t' <<  g_r[m] << endl;
        fs.close();

    };

    void print_data(string const filename, double param) {

        std::ofstream fs;

        fs.open(filename.c_str(), ios::out);
        fs.precision(12);

        fs << param << '\t';
        fs << M1_ee << '\t' << sqrt(M2_ee - pow(M1_ee, 2)) << '\t';
        fs << M1_ep << '\t' << sqrt(M2_ep - pow(M1_ep, 2)) << '\t';
        fs << M1_ek << '\t' << sqrt(M2_ek - pow(M1_ek, 2)) << '\t';
        fs << M1_pr << '\t' << sqrt(M2_pr - pow(M1_pr, 2)) << '\t';
        fs << M1_vol << '\t' << sqrt(M2_vol - pow(M1_vol, 2)) << '\t';
        fs << endl;

        fs.close();

    };

    void print_data_append(string const filename, double param) {

        std::ofstream fs;

        fs.open(filename.c_str(), ios::app);
        fs.precision(12);

        fs << param << '\t';
        fs << M1_ee << '\t' << sqrt(M2_ee - pow(M1_ee, 2)) << '\t';
        fs << M1_ep << '\t' << sqrt(M2_ep - pow(M1_ep, 2)) << '\t';
        fs << M1_ek << '\t' << sqrt(M2_ek - pow(M1_ek, 2)) << '\t';
        fs << M1_pr << '\t' << sqrt(M2_pr - pow(M1_pr, 2)) << '\t';
        fs << M1_vol << '\t' << sqrt(M2_vol - pow(M1_vol, 2)) << '\t';

        fs << endl;

        fs.close();

    };


    void basic_measurements(SOrbitalSystem & system);

    void compute_g_r(Vec<AtomPair> pairs, int Nm, double volume, double r_min, double r_max);

};


void average_MD_Data(MD_Data & accu_data, MD_Data & mddata, int nc);

#endif /* analysis_hpp */
