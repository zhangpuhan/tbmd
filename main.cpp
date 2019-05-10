//
//  main.cpp
//  TB_MD_s_orbital
//
//

#include <iostream>
#include "util.hpp"
#include "potential.hpp"
#include "tbmd.hpp"
#include "analysis.hpp"
#include "test.hpp"

using namespace std;

void GQMD_simulation(MD_param param, MD_Data &accu_data) {

    double l_box = param.l_box;
    double kT = param.kT;
    double mass = param.atom_mass;
    double gamma = param.Langevin_damping;
    int N_atoms = param.N_atoms;
    double dt = param.dlt_t;

    SOrbitalSystem atoms(N_atoms, mass, param.filling_fraction, kT, gamma);
    atoms.set_BD_type(1);

    //atoms.plot_binding_curve("./data/eb.dat");

    double r_s = pow(3. / (4.*M_PI) / param.N_atoms, 1./3.) * l_box;
    cout << "r_s = " << r_s << endl;

    double lt = l_box;
    double lz = l_box;
    vec3 bd_lo, bd_hi;

    atoms.init_random(lt, lz, bd_lo, bd_hi);

    MD_Data mddata;
    double r_min = 0.000001;
    double r_max = 0.25 * (lt + lz);
    double dr = 0.01;
    int num_r = (int) ((r_max - r_min)/dr);
    mddata.init_g_r(num_r, dr);


    //ofstream fs("./data/r.dat");

    int n_save = 40;
    int nx = 0;

    int n_save_ML = 500;
    int nc = 0;

    for(int r=0; r < param.N_steps; r++) {

        cout << "r = " << r << endl;

        atoms.step_NVT(dt);

        double ee = atoms.e_elec();
        double ep = atoms.e_pair();
        double ek = atoms.e_kin();

        //fs << r * dt << '\t' << ee << '\t' << ep << '\t' << ek << '\t';
        //fs << endl;

        if(r % n_save == 0) {

            cout << "r = " << r << endl;

            mddata.basic_measurements(atoms);
            mddata.compute_g_r(atoms.pairs, atoms.numAtoms, atoms.volume, r_min, r_max);


            average_MD_Data(accu_data, mddata, nx);
            nx++;

            //accu_data.print_g_r("./data/gr.dat");

            //atoms.save_configuration("./data/c.dat");

            //atoms.print_spectrum("./data/em.dat");
        }

        if(r % n_save_ML == 0) {

            atoms.save_ML_snapshots("./testdata/g" + std::to_string(nc) + ".csv",
                                    "./testdata/ge" + std::to_string(nc) + ".csv");
            nc++;
        }

    }

    //fs.close();
}


int main(int argc, const char * argv[]) {

    //test_dimer2();

    double k_T          = argc > 1 ? atof(argv[1]) : 0.15; //0.15;
    double num_atoms    = argc > 2 ? atoi(argv[2]) : 50;
    double l_box        = argc > 3 ? atof(argv[3]) : 6.5;

    cout << "kT = " << k_T << endl;
    cout << "N_atoms = " << num_atoms << endl;
    cout << "l_box = " << l_box << endl;

    MD_param sim1;

    sim1.l_box = l_box;
    sim1.kT = k_T;
    sim1.N_atoms = num_atoms;
    sim1.filling_fraction = 0.5;
    sim1.N_steps = 500010;
    sim1.atom_mass = 1;
    sim1.Langevin_damping = 0.05;
    sim1.dlt_t = 0.02;

    MD_Data mdt;

    GQMD_simulation(sim1, mdt);


    return 0;
}
