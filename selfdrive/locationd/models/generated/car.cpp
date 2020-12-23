
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3941128655911992447) {
   out_3941128655911992447[0] = delta_x[0] + nom_x[0];
   out_3941128655911992447[1] = delta_x[1] + nom_x[1];
   out_3941128655911992447[2] = delta_x[2] + nom_x[2];
   out_3941128655911992447[3] = delta_x[3] + nom_x[3];
   out_3941128655911992447[4] = delta_x[4] + nom_x[4];
   out_3941128655911992447[5] = delta_x[5] + nom_x[5];
   out_3941128655911992447[6] = delta_x[6] + nom_x[6];
   out_3941128655911992447[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_1859391793331179110) {
   out_1859391793331179110[0] = -nom_x[0] + true_x[0];
   out_1859391793331179110[1] = -nom_x[1] + true_x[1];
   out_1859391793331179110[2] = -nom_x[2] + true_x[2];
   out_1859391793331179110[3] = -nom_x[3] + true_x[3];
   out_1859391793331179110[4] = -nom_x[4] + true_x[4];
   out_1859391793331179110[5] = -nom_x[5] + true_x[5];
   out_1859391793331179110[6] = -nom_x[6] + true_x[6];
   out_1859391793331179110[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_5793256231659498220) {
   out_5793256231659498220[0] = 1.0;
   out_5793256231659498220[1] = 0.0;
   out_5793256231659498220[2] = 0.0;
   out_5793256231659498220[3] = 0.0;
   out_5793256231659498220[4] = 0.0;
   out_5793256231659498220[5] = 0.0;
   out_5793256231659498220[6] = 0.0;
   out_5793256231659498220[7] = 0.0;
   out_5793256231659498220[8] = 0.0;
   out_5793256231659498220[9] = 1.0;
   out_5793256231659498220[10] = 0.0;
   out_5793256231659498220[11] = 0.0;
   out_5793256231659498220[12] = 0.0;
   out_5793256231659498220[13] = 0.0;
   out_5793256231659498220[14] = 0.0;
   out_5793256231659498220[15] = 0.0;
   out_5793256231659498220[16] = 0.0;
   out_5793256231659498220[17] = 0.0;
   out_5793256231659498220[18] = 1.0;
   out_5793256231659498220[19] = 0.0;
   out_5793256231659498220[20] = 0.0;
   out_5793256231659498220[21] = 0.0;
   out_5793256231659498220[22] = 0.0;
   out_5793256231659498220[23] = 0.0;
   out_5793256231659498220[24] = 0.0;
   out_5793256231659498220[25] = 0.0;
   out_5793256231659498220[26] = 0.0;
   out_5793256231659498220[27] = 1.0;
   out_5793256231659498220[28] = 0.0;
   out_5793256231659498220[29] = 0.0;
   out_5793256231659498220[30] = 0.0;
   out_5793256231659498220[31] = 0.0;
   out_5793256231659498220[32] = 0.0;
   out_5793256231659498220[33] = 0.0;
   out_5793256231659498220[34] = 0.0;
   out_5793256231659498220[35] = 0.0;
   out_5793256231659498220[36] = 1.0;
   out_5793256231659498220[37] = 0.0;
   out_5793256231659498220[38] = 0.0;
   out_5793256231659498220[39] = 0.0;
   out_5793256231659498220[40] = 0.0;
   out_5793256231659498220[41] = 0.0;
   out_5793256231659498220[42] = 0.0;
   out_5793256231659498220[43] = 0.0;
   out_5793256231659498220[44] = 0.0;
   out_5793256231659498220[45] = 1.0;
   out_5793256231659498220[46] = 0.0;
   out_5793256231659498220[47] = 0.0;
   out_5793256231659498220[48] = 0.0;
   out_5793256231659498220[49] = 0.0;
   out_5793256231659498220[50] = 0.0;
   out_5793256231659498220[51] = 0.0;
   out_5793256231659498220[52] = 0.0;
   out_5793256231659498220[53] = 0.0;
   out_5793256231659498220[54] = 1.0;
   out_5793256231659498220[55] = 0.0;
   out_5793256231659498220[56] = 0.0;
   out_5793256231659498220[57] = 0.0;
   out_5793256231659498220[58] = 0.0;
   out_5793256231659498220[59] = 0.0;
   out_5793256231659498220[60] = 0.0;
   out_5793256231659498220[61] = 0.0;
   out_5793256231659498220[62] = 0.0;
   out_5793256231659498220[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_8865508588157271509) {
   out_8865508588157271509[0] = state[0];
   out_8865508588157271509[1] = state[1];
   out_8865508588157271509[2] = state[2];
   out_8865508588157271509[3] = state[3];
   out_8865508588157271509[4] = state[4];
   out_8865508588157271509[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_8865508588157271509[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_8865508588157271509[7] = state[7];
}
void F_fun(double *state, double dt, double *out_4204984779909159995) {
   out_4204984779909159995[0] = 1;
   out_4204984779909159995[1] = 0;
   out_4204984779909159995[2] = 0;
   out_4204984779909159995[3] = 0;
   out_4204984779909159995[4] = 0;
   out_4204984779909159995[5] = 0;
   out_4204984779909159995[6] = 0;
   out_4204984779909159995[7] = 0;
   out_4204984779909159995[8] = 0;
   out_4204984779909159995[9] = 1;
   out_4204984779909159995[10] = 0;
   out_4204984779909159995[11] = 0;
   out_4204984779909159995[12] = 0;
   out_4204984779909159995[13] = 0;
   out_4204984779909159995[14] = 0;
   out_4204984779909159995[15] = 0;
   out_4204984779909159995[16] = 0;
   out_4204984779909159995[17] = 0;
   out_4204984779909159995[18] = 1;
   out_4204984779909159995[19] = 0;
   out_4204984779909159995[20] = 0;
   out_4204984779909159995[21] = 0;
   out_4204984779909159995[22] = 0;
   out_4204984779909159995[23] = 0;
   out_4204984779909159995[24] = 0;
   out_4204984779909159995[25] = 0;
   out_4204984779909159995[26] = 0;
   out_4204984779909159995[27] = 1;
   out_4204984779909159995[28] = 0;
   out_4204984779909159995[29] = 0;
   out_4204984779909159995[30] = 0;
   out_4204984779909159995[31] = 0;
   out_4204984779909159995[32] = 0;
   out_4204984779909159995[33] = 0;
   out_4204984779909159995[34] = 0;
   out_4204984779909159995[35] = 0;
   out_4204984779909159995[36] = 1;
   out_4204984779909159995[37] = 0;
   out_4204984779909159995[38] = 0;
   out_4204984779909159995[39] = 0;
   out_4204984779909159995[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_4204984779909159995[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_4204984779909159995[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4204984779909159995[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4204984779909159995[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_4204984779909159995[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_4204984779909159995[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_4204984779909159995[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_4204984779909159995[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_4204984779909159995[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_4204984779909159995[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4204984779909159995[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4204984779909159995[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_4204984779909159995[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_4204984779909159995[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_4204984779909159995[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4204984779909159995[56] = 0;
   out_4204984779909159995[57] = 0;
   out_4204984779909159995[58] = 0;
   out_4204984779909159995[59] = 0;
   out_4204984779909159995[60] = 0;
   out_4204984779909159995[61] = 0;
   out_4204984779909159995[62] = 0;
   out_4204984779909159995[63] = 1;
}
void h_25(double *state, double *unused, double *out_3903272388947322771) {
   out_3903272388947322771[0] = state[6];
}
void H_25(double *state, double *unused, double *out_6002909593856940474) {
   out_6002909593856940474[0] = 0;
   out_6002909593856940474[1] = 0;
   out_6002909593856940474[2] = 0;
   out_6002909593856940474[3] = 0;
   out_6002909593856940474[4] = 0;
   out_6002909593856940474[5] = 0;
   out_6002909593856940474[6] = 1;
   out_6002909593856940474[7] = 0;
}
void h_24(double *state, double *unused, double *out_8819236755190351777) {
   out_8819236755190351777[0] = state[4];
   out_8819236755190351777[1] = state[5];
}
void H_24(double *state, double *unused, double *out_7453449519679828079) {
   out_7453449519679828079[0] = 0;
   out_7453449519679828079[1] = 0;
   out_7453449519679828079[2] = 0;
   out_7453449519679828079[3] = 0;
   out_7453449519679828079[4] = 1;
   out_7453449519679828079[5] = 0;
   out_7453449519679828079[6] = 0;
   out_7453449519679828079[7] = 0;
   out_7453449519679828079[8] = 0;
   out_7453449519679828079[9] = 0;
   out_7453449519679828079[10] = 0;
   out_7453449519679828079[11] = 0;
   out_7453449519679828079[12] = 0;
   out_7453449519679828079[13] = 1;
   out_7453449519679828079[14] = 0;
   out_7453449519679828079[15] = 0;
}
void h_30(double *state, double *unused, double *out_4504690546560686325) {
   out_4504690546560686325[0] = state[4];
}
void H_30(double *state, double *unused, double *out_2829935568903427862) {
   out_2829935568903427862[0] = 0;
   out_2829935568903427862[1] = 0;
   out_2829935568903427862[2] = 0;
   out_2829935568903427862[3] = 0;
   out_2829935568903427862[4] = 1;
   out_2829935568903427862[5] = 0;
   out_2829935568903427862[6] = 0;
   out_2829935568903427862[7] = 0;
}
void h_26(double *state, double *unused, double *out_8940380259488019646) {
   out_8940380259488019646[0] = state[7];
}
void H_26(double *state, double *unused, double *out_298316438917731201) {
   out_298316438917731201[0] = 0;
   out_298316438917731201[1] = 0;
   out_298316438917731201[2] = 0;
   out_298316438917731201[3] = 0;
   out_298316438917731201[4] = 0;
   out_298316438917731201[5] = 0;
   out_298316438917731201[6] = 0;
   out_298316438917731201[7] = 1;
}
void h_27(double *state, double *unused, double *out_2012358364786539565) {
   out_2012358364786539565[0] = state[3];
}
void H_27(double *state, double *unused, double *out_1542353581066802550) {
   out_1542353581066802550[0] = 0;
   out_1542353581066802550[1] = 0;
   out_1542353581066802550[2] = 0;
   out_1542353581066802550[3] = 1;
   out_1542353581066802550[4] = 0;
   out_1542353581066802550[5] = 0;
   out_1542353581066802550[6] = 0;
   out_1542353581066802550[7] = 0;
}
void h_29(double *state, double *unused, double *out_8375952343088355426) {
   out_8375952343088355426[0] = state[1];
}
void H_29(double *state, double *unused, double *out_3404494182758964037) {
   out_3404494182758964037[0] = 0;
   out_3404494182758964037[1] = 1;
   out_3404494182758964037[2] = 0;
   out_3404494182758964037[3] = 0;
   out_3404494182758964037[4] = 0;
   out_3404494182758964037[5] = 0;
   out_3404494182758964037[6] = 0;
   out_3404494182758964037[7] = 0;
}
void h_28(double *state, double *unused, double *out_6231938022780803683) {
   out_6231938022780803683[0] = state[5];
   out_6231938022780803683[1] = state[6];
}
void H_28(double *state, double *unused, double *out_435617041443333213) {
   out_435617041443333213[0] = 0;
   out_435617041443333213[1] = 0;
   out_435617041443333213[2] = 0;
   out_435617041443333213[3] = 0;
   out_435617041443333213[4] = 0;
   out_435617041443333213[5] = 1;
   out_435617041443333213[6] = 0;
   out_435617041443333213[7] = 0;
   out_435617041443333213[8] = 0;
   out_435617041443333213[9] = 0;
   out_435617041443333213[10] = 0;
   out_435617041443333213[11] = 0;
   out_435617041443333213[12] = 0;
   out_435617041443333213[13] = 0;
   out_435617041443333213[14] = 1;
   out_435617041443333213[15] = 0;
}
}

extern "C"{
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
