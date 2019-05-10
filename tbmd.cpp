//
//  tbmd.cpp
//  GQMD-Hubbard-liquid
//
//

#include <iostream>
#include <fstream>
#include <cassert>

#include "util.hpp"
#include "potential.hpp"
#include "tbmd.hpp"

SOrbitalSystem::SOrbitalSystem(int nAtoms, double atomMass, double filling, double temperature, double gm) {
    numAtoms = nAtoms;
    position = Vec<vec3>(nAtoms);
    velocity = Vec<vec3>(nAtoms);
    force = Vec<vec3>(nAtoms);

    mu_n = Vec<double>(nAtoms);
    n0 = Vec<double>(nAtoms);

    mass = atomMass;
    kT = temperature;
    gamma = gm;

    ee = ep = ek = 0;

    kT_elec = kT;

    boundary_type = 1;      // default: 1 corresponds to periodic BC
    langevin_dym = 1;

    filling_fraction = filling;
    num_filling = (int) (filling_fraction * nAtoms);

    dimH = nAtoms;

    eigE = arma::vec(dimH);
    eig_mode = arma::cx_mat(dimH, dimH);
    fd_factor = arma::vec(dimH);

    //    Hamiltonian = arma::sp_cx_mat(nAtoms, nAtoms);
    //    rho = arma::sp_cx_mat(nAtoms, nAtoms);


    potential = new(ToyModelSOribtal);

    rng = RNG(seed());

    P_tensor.resize(3,3);
}

void SOrbitalSystem::init_bcc(double rnn, int M, vec3& bdsLo, vec3& bdsHi) {
    position.resize(numAtoms);

    double a0 = (2./sqrt(3.)) * rnn;

    bdsLo = {0, 0, 0};
    bdsHi = {M*a0, M*a0, M*a0};

    int ir = 0;
    for(int x=0; x<M; x++) {
        for(int y=0; y<M; y++) {
            for(int z=0; z<M; z++) {

                position[ir] = vec3 {x*a0, y*a0, z*a0};
                position[ir+1] = position[ir] + vec3 {0.5*a0, 0.5*a0, 0.5*a0};

                ir+=2;
            }
        }
    }


    system_box = bdsHi - bdsLo;
    cout << "system box = " << system_box << endl;
    volume = system_box.x * system_box.y * system_box.z;
    cout << "volumn = " << volume << endl;

    for(int i=0; i<numAtoms; i++) {
        velocity[i] = vec3 {0, 0, 0};
        force[i] = vec3 {0, 0, 0};
    }

    double sgm_v = sqrt(kT / mass);
    std::normal_distribution<double> rn;    // default mean = 0, var = 1
    for(int i=0; i<numAtoms; i++) {
        for(int k=0; k<3; k++) velocity[i](k) = sgm_v * rn(rng);
    }

    pos_init = Vec<vec3>(numAtoms);
    for(int i=0; i<numAtoms; i++) pos_init[i] = position[i];
}

void SOrbitalSystem::init_sc(double rnn, int M, vec3& bdsLo, vec3& bdsHi) {
    position.resize(numAtoms);

    double a0 = rnn;

    bdsLo = {0, 0, 0};
    bdsHi = {M*a0, M*a0, M*a0};

    int ir = 0;
    for(int x=0; x<M; x++) {
        for(int y=0; y<M; y++) {
            for(int z=0; z<M; z++) {

                position[ir] = vec3 {x*a0, y*a0, z*a0};

                ir++;
            }
        }
    }


    system_box = bdsHi - bdsLo;
    cout << "system box = " << system_box << endl;
    volume = system_box(0) * system_box(1) * system_box(2);
    cout << "volumn = " << volume << endl;

    for(int i=0; i<numAtoms; i++) {
        velocity[i] = vec3 {0, 0, 0};
        force[i] = vec3 {0, 0, 0};
    }

    double sgm_v = sqrt(kT / mass);
    std::normal_distribution<double> rn;    // default mean = 0, var = 1
    for(int i=0; i<numAtoms; i++) {
        for(int k=0; k<3; k++) velocity[i](k) = sgm_v * rn(rng);
    }

    pos_init = Vec<vec3>(numAtoms);
    for(int i=0; i<numAtoms; i++) pos_init[i] = position[i];
}

void SOrbitalSystem::init_linear_chain(double rnn, int L, vec3& bdsLo, vec3& bdsHi) {
    position.resize(numAtoms);

    double a0 = rnn;

    bdsLo = {0, 0, 0};
    bdsHi = {L*a0, a0, a0};

    for(int x=0; x<L; x++) position[x] = vec3 {x * a0, 0, 0};

    std::normal_distribution<double> rd;

    double dm = 0.0001 * a0;
    for(int i=0; i<L; i++) position[i] += vec3 {dm * rd(rng), dm * rd(rng), dm * rd(rng) };

    system_box = bdsHi - bdsLo;
    cout << "system box = " << system_box << endl;
    volume = system_box(0) * system_box(1) * system_box(2);
    cout << "volumn = " << volume << endl;

    for(int i=0; i<numAtoms; i++) {
        velocity[i] = vec3 {0, 0, 0};
        force[i] = vec3 {0, 0, 0};
    }

    double sgm_v = sqrt(kT / mass);
    std::normal_distribution<double> rn;    // default mean = 0, var = 1
    for(int i=0; i<numAtoms; i++) {
        for(int k=0; k<3; k++) velocity[i](k) = sgm_v * rn(rng);
    }

    pos_init = Vec<vec3>(numAtoms);
    for(int i=0; i<numAtoms; i++) pos_init[i] = position[i];
}

void SOrbitalSystem::init_random(double Lt, double Lz, vec3& bdsLo, vec3& bdsHi) {
    //RNG rg(0);
    std::uniform_real_distribution<double> rd;

    bdsLo = {0, 0, 0};
    bdsHi = {Lt, Lt, Lz};
    position.resize(numAtoms);

    for(int i=0; i<numAtoms; i++) {
        position[i] = vec3 {rd(rng) * Lt, rd(rng) * Lt, rd(rng) * Lz};
    }

    //        position[0] = vec3 {1*bohr, 0, 0};
    //        position[1] = vec3 {3.*bohr, 0, 0};

    system_box = bdsHi - bdsLo;
    cout << "system box = " << system_box << endl;
    volume = system_box(0) * system_box(1) * system_box(2);
    cout << "volumn = " << volume << endl;

    for(int i=0; i<numAtoms; i++) {
        velocity[i] = vec3 {0, 0, 0};
        force[i] = vec3 {0, 0, 0};
    }

//    double sgm_v = sqrt(kT / mass);
//    std::normal_distribution<double> rn;    // default mean = 0, var = 1
//    for(int i=0; i<numAtoms; i++) {
//        for(int k=0; k<3; k++) velocity[i](k) = sgm_v * rn(rng);
//    }

    pos_init = Vec<vec3>(numAtoms);
    for(int i=0; i<numAtoms; i++) pos_init[i] = position[i];
}

void SOrbitalSystem::init_random_dimers(double Lt, double Lz, vec3& bdsLo, vec3& bdsHi) {
    //RNG rg(0);
    std::uniform_real_distribution<double> rd;

    bdsLo = {0, 0, 0};
    bdsHi = {Lt, Lt, Lz};
    position.resize(numAtoms);

    double d_dm = 1.0 * bohr;
    for(int i=0; i<numAtoms; i+=2) {
        position[i] = vec3 {rd(rng) * Lt, rd(rng) * Lt, rd(rng) * Lz};
        position[i+1] = position[i] + vec3 {rd(rng) * d_dm, rd(rng) * d_dm, rd(rng) * d_dm};
    }

    //        position[0] = vec3 {1*bohr, 0, 0};
    //        position[1] = vec3 {3.*bohr, 0, 0};

    system_box = bdsHi - bdsLo;
    cout << "system box = " << system_box << endl;
    volume = system_box(0) * system_box(1) * system_box(2);

    for(int i=0; i<numAtoms; i++) {
        velocity[i] = vec3 {0, 0, 0};
        force[i] = vec3 {0, 0, 0};
    }

    double sgm_v = sqrt(kT / mass);
    std::normal_distribution<double> rn;    // default mean = 0, var = 1
    for(int i=0; i<numAtoms; i++) {
        for(int k=0; k<3; k++) velocity[i](k) = sgm_v * rn(rng);
    }

    pos_init = Vec<vec3>(numAtoms);
    for(int i=0; i<numAtoms; i++) pos_init[i] = position[i];
}

void SOrbitalSystem::find_all_pairs(void) {
    if(boundary_type == 0)
        pairs = allPairs(position, potential->rcut());
    else
        pairs = allPairsPeriodic(position, potential->rcut(), system_box);
}

void SOrbitalSystem::find_all_pairs(double cutoff) {
    if(boundary_type == 0)
        pairs = allPairs(position, cutoff);
    else
        pairs = allPairsPeriodic(position, cutoff, system_box);
}

void SOrbitalSystem::built_Hamiltonian(void) {
    Hamiltonian = potential->build_Hamiltonian(numAtoms, pairs);
    //Hamiltonian = potential->buildHamiltonian(numAtoms, pairs, system_box);
}

double SOrbitalSystem::compute_avg_density(double x) {
    double sum = 0;
    for(int i=0; i<eigE.size(); i++) {
        sum += fermi_density(eigE(i), kT_elec, x);
    }
    return sum/((double) numAtoms);
}

double SOrbitalSystem::compute_chemical_potential(void) {

    double x1 = eigE(0);
    double x2 = eigE(eigE.size()-1);

    int max_bisection = 100;
    double eps_bisection = 1.e-12;

    int iter = 0;
    while(iter < max_bisection || fabs(x2 - x1) > eps_bisection) {

        double xm = 0.5*(x1 + x2);
        double density = compute_avg_density(xm);

        if(density <= filling_fraction) x1 = xm;
        else x2 = xm;

        iter++;
    }

    return 0.5*(x1 + x2);
}

//kpm::EnergyScale SOrbitalSystem::computeEnergyScale(double extra, double tolerance) {
//    arma::cx_vec eigval;
//    eigs_gen(eigval, Hamiltonian, 1, "sr", tolerance);
//    double eig_min = eigval(0).real();
//    eigs_gen(eigval, Hamiltonian, 1, "lr", tolerance);
//    double eig_max = eigval(0).real();
//    double slack = extra * (eig_max - eig_min);
//    e_min = eig_min;
//    return {eig_min-slack, eig_max+slack};
//}

arma::cx_mat SOrbitalSystem::compute_density_matrix(arma::cx_mat & H_elec) {

    arma::cx_mat rho_elec(dimH, dimH);

    arma::cx_mat eigvec;
    auto Hd = H_elec; //fkpm::sparseToDense(H_elec);
    arma::eig_sym(eigE, eigvec, Hd);
    eig_mode = eigvec;

    mu = compute_chemical_potential();

    //cout << "Hd: " << endl;
    //cout << Hd << endl;

    for(int i=0; i<dimH; i++) fd_factor(i) = fermi_density(eigE(i), kT_elec, mu);

    rho_elec.zeros();
    for(int a=0; a<dimH; a++)
        for(int b=a; b<dimH; b++) {

            cx_double sum = 0;
            for(int m=0; m<eigE.size(); m++) {
                sum += fd_factor(m) * conj(eigvec(a, m)) * eigvec(b, m);
            }
            rho_elec(a, b) = sum;
            if(a != b) rho_elec(b, a) = conj(sum);
        }

    //cout << "rho: " << endl;
    //cout << rho_elec << endl;

    return rho_elec;
}


void SOrbitalSystem::compute_density_matrix(void) {
    rho = compute_density_matrix(Hamiltonian);
}

void SOrbitalSystem::compute_forces(void) {

    potential->force(rho, pairs, force, virial);

    //potential->force(potential->numSpins() * rho, pairs, force, virial, system_box);  // old version CHECK !!!!!!
    // the numSpins() is absorbed into the calculation in potential->force(...) !!!

}

void SOrbitalSystem::compute_forces_elec(void) {

    potential->force_elec(rho, pairs, force, virial);

}

void SOrbitalSystem::move_atoms(double dt) {
    for(int i=0; i<numAtoms; i++) {
        vec3 dlt = 0.5 * (force[i]/mass) * dt;
        dlt += velocity[i];
        velocity[i] = dlt;  // velecity at t + 0.5*dt
        dlt *= dt;
        position[i] += dlt;
    }
}

void SOrbitalSystem::integrate_forces(double dt) {
    for(int i=0; i<numAtoms; i++) {
        velocity[i] += 0.5 * (force[i]/mass) * dt;
    }
}

void SOrbitalSystem::integrate_Langevin(double dt) {
    //RNG rng(seed());
    std::normal_distribution<double> rd;    // default mean = 0, var = 1

    double alpha2 = exp(-gamma * dt);
    double sigma2 = sqrt((1-pow(alpha2, 2))*kT/mass);

    for(int i=0; i<numAtoms; i++) {
        velocity[i] = alpha2 * velocity[i] + sigma2 * vec3(rd(rng), rd(rng), rd(rng));
    }
}

void SOrbitalSystem::solve_TB_Hamiltonian(void) {

    built_Hamiltonian();
    compute_density_matrix();

}

void SOrbitalSystem::step_NVT(double dt) {

    // velocity verlet integration with Langevin damping
    move_atoms(dt);

    find_all_pairs();

    solve_TB_Hamiltonian();

    compute_forces();

    integrate_forces(dt);

    if(langevin_dym == 1)
        integrate_Langevin(dt);
}

double SOrbitalSystem::e_elec(void) {
    double sum1 = 0;

    //for(int i=0; i<numFilling; i++) sum1 += eigE(i);

    //cout << "mu = " << mu << endl;

    //        for(int i=0; i<dimH; i++) {
    //            sum1 += eigE(i) * fermiDensity(eigE(i), kT, mu);
    //        }

    //auto rr = Hamiltonian.to_arma() * rho.to_arma();
    //cout << trace(rr).real() << endl;

    auto H_elec = Hamiltonian;
    auto rho_elec = rho;

    auto rr = H_elec * rho_elec;

    //    cout << H_elec << endl;
    //    cout << rho_elec << endl;
    //    cout << rr << endl;

    sum1 = trace(rr).real();

    //sum1 = electronicGrandEnergy(gm, energyScale, kT_electron, mu) + Hamiltonian.n_rows * fillingFraction * mu;

    // ------- on-site Coloumb repulsion -------------

    ee = potential->numSpins() * sum1 / ((double) numAtoms);

    return ee;
}

double SOrbitalSystem::e_kin(void) {
    double sum = 0;
    for(int i=0; i<numAtoms; i++) sum += velocity[i].norm2();
    ek = (0.5 * mass * sum) / ((double) numAtoms);
    return ek;
}

double SOrbitalSystem::e_pair(void) {
    ep = potential->pair_energy(pairs) / ((double) numAtoms);
    return ep;
}

void SOrbitalSystem::compute_cm_velocity(void) {
    cm_velocity = {0, 0, 0};
    for(int i=0; i<numAtoms; i++) cm_velocity += velocity[i];
    cm_velocity /= ((double) numAtoms);
}

void SOrbitalSystem::subtract_cm_velocity(void) {
    for(int i=0; i<numAtoms; i++) velocity[i] -= cm_velocity;
}

void SOrbitalSystem::compute_pressure_tensor(void) {

    P_tensor.zeros();

    for(int i=0; i<numAtoms; i++) {
        for(int x=0; x<3; x++)
            for(int y=0; y<3; y++) P_tensor(x, y) += mass * velocity[i](x) * velocity[i](y);
    }

    compute_forces();

    for (auto & p : pairs) {
        if (p.delta.norm2() < pow(potential->rcut(),2)) {
            int i = p.index1;
            int j = p.index2;

            for(int x=0; x<3; x++)
                for(int y=0; y<3; y++) P_tensor(x, y) += 0.5 * p.delta(x) * (force[j](y) - force[i](y));

        }
    }

    P_tensor /= volume;
}

double SOrbitalSystem::compute_virial(void) {

    double EPS_V = 1.e-3;
    Vec<AtomPair> pairs_x = pairs;

    for(int i=0; i<pairs.size(); i++)
        pairs_x[i].delta = (1. + EPS_V) * pairs[i].delta;

    Hamiltonian = potential->build_Hamiltonian(numAtoms, pairs_x);
    compute_density_matrix();

    double U1 = e_elec() * numAtoms + potential->pair_energy(pairs_x);

    for(int i=0; i<pairs.size(); i++)
        pairs_x[i].delta = (1. - EPS_V) * pairs[i].delta;

    Hamiltonian = potential->build_Hamiltonian(numAtoms, pairs_x);
    compute_density_matrix();

    double U2 = e_elec() * numAtoms + potential->pair_energy(pairs_x);

    virial = -(U1 - U2) / (pow(1.+EPS_V, 3) - pow(1.-EPS_V, 3));

    return virial;
}

double SOrbitalSystem::compute_pressure(void) {

    //        compute_virial();
    //        pressure = (virial + (2./3.) * ek * numAtoms) / volume;     // 3-dimensions

    pressure = (virial/3. + (2./3.) * ek * numAtoms) / volume;     // 3-dimensions

    return pressure;
}

void SOrbitalSystem::save_ML_snapshots(string const filename, string const filename_2) {

    Vec<vec3> force_elec = Vec<vec3>(numAtoms);

    potential->force_elec(rho, pairs, force_elec, virial);

    std::ofstream fs;
    std::ofstream fs_2;

    fs.open(filename.c_str(), ios::out);
    fs.precision(12);
    fs_2.open(filename_2.c_str(), ios::out);
    fs_2.precision(12);


    vec3 bdsLo = {0, 0, 0};
    vec3 bdsHi = {system_box.x, system_box.y, system_box.z};


    double e_electron = e_elec();
    fs_2 << e_electron << endl;
    fs_2.close();

    for(int i=0; i<numAtoms; i++) {

        vec3 pos = wrapPosition(position[i], bdsLo, bdsHi);

        fs << pos.x << "," << pos.y << "," << pos.z << ",";
        fs << force_elec[i].x << "," << force_elec[i].y << "," << force_elec[i].z;
        //fs << real(rho(i, i)) << '\t';
        fs << endl;
    }

    fs.close();
}

void SOrbitalSystem::save_configuration(string const filename) {

    std::ofstream fs;

    fs.open(filename.c_str(), ios::out);
    fs.precision(12);

    for(int i=0; i<numAtoms; i++) {

        fs << position[i].x << '\t' << position[i].y << '\t' << position[i].z << '\t';
        fs << velocity[i].x << '\t' << velocity[i].y << '\t' << velocity[i].z << endl;
    }
    fs.close();
}

void SOrbitalSystem::read_configuration(double Lt, double Lz, vec3& bdsLo, vec3& bdsHi, const string filename) {

    position.resize(numAtoms);

    std::ifstream fs;
    fs.open(filename.c_str(), ios::in);
    for(int i=0; i<numAtoms; i++) {

        // ====== for reading from config.dat
        fs >> position[i].x >> position[i].y >> position[i].z >> velocity[i].x >> velocity[i].y >> velocity[i].z;

        cout << position[i].x << '\t' << position[i].y << '\t' << position[i].z << endl;
    }

    bdsLo = {0, 0, 0};
    bdsHi = {Lt, Lt, Lz};

    system_box = bdsHi - bdsLo;
    cout << "system box = " << system_box << endl;
    volume = system_box(0) * system_box(1) * system_box(2);
    cout << "volumn = " << volume << endl;

    for(int i=0; i<numAtoms; i++) {
        force[i] = vec3 {0, 0, 0};
    }


    // IF read from "config.dat", then comment off the following initialization code for velocity
    /*
    double sgm_v = sqrt(kT / mass);
    std::normal_distribution<double> rn;    // default mean = 0, var = 1
    for(int i=0; i<numAtoms; i++) {
        for(int k=0; k<3; k++) velocity[i](k) = sgm_v * rn(rng);
    }
    */


    pos_init = Vec<vec3>(numAtoms);
    for(int i=0; i<numAtoms; i++) pos_init[i] = position[i];

}

void SOrbitalSystem::plot_binding_curve(const string filename) {

    std::ofstream fs;

    fs.open(filename.c_str(), ios::out);
    fs.precision(12);

    double dr = 0.01;
    for(double r = 1.e-5; r<=potential->rcut(); r+= dr) {
        double e_t = potential->hopping(r)(0,0);
        double e_p = potential->phi(r);
        fs << r << '\t' << e_t << '\t' << e_p << '\t' << 2 * e_t + e_p << endl;
    }

    fs.close();
}

void SOrbitalSystem::print_spectrum(const string filename) {

    std::ofstream fs;

    fs.open(filename.c_str(), ios::out);
    fs.precision(12);

    for(int m=0; m<dimH; m++) {
        fs << eigE(m) << '\t' << fd_factor(m) << endl;
    }

    fs.close();
}



