//
//  tbmd.hpp
//  GQMD-Hubbard-liquid
//
//

#ifndef tbmd_hpp
#define tbmd_hpp

#include <iostream>
#include <fstream>

#include "util.hpp"
#include "potential.hpp"


class SOrbitalSystem {
public:
    int numAtoms;

    Vec<vec3> position;
    Vec<vec3> velocity;
    Vec<vec3> force;

    Vec<vec3> pos_init;
    Vec<vec3> v_init;

    Vec<double> mu_n;
    Vec<double> n0;

    vec3 cm_velocity;

    double mass;
    double kT;
    double gamma;   // Langevin damping

    double filling_fraction;
    int num_filling;
    double mu;              // chemical potential

    double kT_elec;

    arma::vec eigE;
    arma::cx_mat eig_mode;

    ToyModelSOribtal *potential;

    RNG rng;
    std::random_device seed;

    int dimH;

    arma::cx_mat Hamiltonian;
    arma::cx_mat rho;

    vec3 system_box;
    double volume;

    // constructor:
    SOrbitalSystem(int nAtoms, double atomMass, double filling, double temperature, double gm);

    // ============================================================================
    // initialization of various lattices:
    void init_bcc(double rs, int M, vec3& bdsLo, vec3& bdsHi);
    void init_sc(double rs, int M, vec3& bdsLo, vec3& bdsHi);
    void init_linear_chain(double rnn, int M, vec3& bdsLo, vec3& bdsHi);
    void init_random(double Lt, double Lz, vec3& bdsLo, vec3& bdsHi);
    void init_random_dimers(double Lt, double Lz, vec3& bdsLo, vec3& bdsHi);

    // ============================================================================
    // set parameters:
    int boundary_type;   // 0 : Open BC,  1 : Periodic BC
    int langevin_dym;

    void set_BD_type(int t) {boundary_type = t;};   // default = 1
    void set_LD_dynm(int t) {langevin_dym = t;};    // default = 1

    // ============================================================================
    Vec<AtomPair> pairs;

    void find_all_pairs(void);
    void find_all_pairs(double cutoff);

    // ============================================================================
    void built_Hamiltonian(void);

    //fkpm::EnergyScale computeEnergyScale(double extra, double tolerance);

    // ============================================================================
    // for single-particle density matrix:
    double compute_chemical_potential(void);
    double compute_avg_density(double x);
    arma::cx_mat compute_density_matrix(arma::cx_mat & H_elec);
    void compute_density_matrix(void);
    arma::vec fd_factor;

    // ============================================================================
    // Solving bare Hamiltonian
    void solve_TB_Hamiltonian(void);

    // ============================================================================
    // MD-related routines
    void move_atoms(double dt);
    void compute_forces(void);
    void integrate_forces(double dt);
    void integrate_Langevin(double dt);
    void step_NVT(double dt);


    // For machine-learning training

    void compute_forces_elec(void);

    // ============================================================================
    // compute energies
    double ee, ep, ek;

    double e_kin(void);
    double e_elec(void);
    double e_pair(void);

    // ============================================================================
    // measurements:

    void compute_cm_velocity(void);
    void subtract_cm_velocity(void);

    double virial;
    double pressure;
    arma::mat P_tensor;
    void compute_pressure_tensor(void);
    double compute_virial(void);
    double compute_pressure(void);

    // ============================================================================

    void save_ML_snapshots(string const filename, string const filename_2);
    void save_configuration(string const filename);
    void read_configuration(double Lt, double Lz, vec3& bdsLo, vec3& bdsHi, string const filename);

    // ============================================================================

    void plot_binding_curve(string const filename);

    // ============================================================================

    void print_spectrum(const string filename);

    void print_ga_data(const string filename);

};

#endif /* tbmd_hpp */
