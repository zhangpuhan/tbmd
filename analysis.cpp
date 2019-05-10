//
//  analysis.cpp
//  GQMD-Hubbard-liquid
//
//

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "util.hpp"
#include "potential.hpp"
#include "tbmd.hpp"
#include "analysis.hpp"

void average_MD_Data(MD_Data & accu_data, MD_Data & mddata, int nc) {

    if(nc == 0) {
        cout << "resetting accumulated MD_data." << endl;
        accu_data.reset();

        accu_data.init_g_r(mddata.num_r, mddata.dr);

    }



    accu_data.M1_ee = (accu_data.M1_ee * nc + mddata.M1_ee) / (nc + 1.);
    accu_data.M2_ee = (accu_data.M2_ee * nc + mddata.M2_ee) / (nc + 1.);

    accu_data.M1_ep = (accu_data.M1_ep * nc + mddata.M1_ep) / (nc + 1.);
    accu_data.M2_ep = (accu_data.M2_ep * nc + mddata.M2_ep) / (nc + 1.);

    accu_data.M1_ek = (accu_data.M1_ek * nc + mddata.M1_ek) / (nc + 1.);
    accu_data.M2_ek = (accu_data.M2_ek * nc + mddata.M2_ek) / (nc + 1.);

    accu_data.M1_pr = (accu_data.M1_pr * nc + mddata.M1_pr) / (nc + 1.);
    accu_data.M2_pr = (accu_data.M2_pr * nc + mddata.M2_pr) / (nc + 1.);

    accu_data.M1_vl = (accu_data.M1_vl * nc + mddata.M1_vl) / (nc + 1.);
    accu_data.M2_vl = (accu_data.M2_vl * nc + mddata.M2_vl) / (nc + 1.);

    accu_data.M1_vol = (accu_data.M1_vol * nc + mddata.M1_vol) / (nc + 1.);
    accu_data.M2_vol = (accu_data.M2_vol * nc + mddata.M2_vol) / (nc + 1.);

    for(unsigned m=0; m<accu_data.g_r.size(); m++) {
        accu_data.g_r[m] = (accu_data.g_r[m] * nc + mddata.g_r[m]) / (nc + 1.);
    }

}

void MD_Data::basic_measurements(SOrbitalSystem & system) {

    double ee = system.e_elec();
    double ep = system.e_pair();
    double ek = system.e_kin();

    double pr = system.compute_pressure();
    double vl = system.virial;

    double vol = system.volume;

    M1_ee = ee;
    M2_ee = pow(ee, 2);

    M1_ep = ep;
    M2_ep = pow(ep, 2);

    M1_ek = ek;
    M2_ek = pow(ek, 2);

    M1_pr = pr;
    M2_pr = pow(pr, 2);

    M1_vl = vl;
    M2_vl = pow(vl, 2);

    M1_vol = vol;
    M2_vol = pow(vol, 2);

}

void MD_Data::compute_g_r(Vec<AtomPair> pairs, int Nm, double volume, double r_min, double r_max) {

    for(int m=0; m<num_r; m++) g_r[m] = 0;

    for(int i=0; i<pairs.size(); i++) {
        double rij = pairs[i].delta.norm();
        if(rij < r_max && rij > r_min) {
            int m = (int) ((rij - r_min)/dr);

            if(m<num_r && m>=0)
                g_r[m]++;
        }
    }

    for(int m=0; m<num_r; m++) {
        double rm = (m + 0.5) * dr;
        g_r[m] *= volume / (2. * 3.1415926 * dr * pow(Nm * rm, 2));
    }
}

