/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3941128655911992447);
void inv_err_fun(double *nom_x, double *true_x, double *out_1859391793331179110);
void H_mod_fun(double *state, double *out_5793256231659498220);
void f_fun(double *state, double dt, double *out_8865508588157271509);
void F_fun(double *state, double dt, double *out_4204984779909159995);
void h_25(double *state, double *unused, double *out_3903272388947322771);
void H_25(double *state, double *unused, double *out_6002909593856940474);
void h_24(double *state, double *unused, double *out_8819236755190351777);
void H_24(double *state, double *unused, double *out_7453449519679828079);
void h_30(double *state, double *unused, double *out_4504690546560686325);
void H_30(double *state, double *unused, double *out_2829935568903427862);
void h_26(double *state, double *unused, double *out_8940380259488019646);
void H_26(double *state, double *unused, double *out_298316438917731201);
void h_27(double *state, double *unused, double *out_2012358364786539565);
void H_27(double *state, double *unused, double *out_1542353581066802550);
void h_29(double *state, double *unused, double *out_8375952343088355426);
void H_29(double *state, double *unused, double *out_3404494182758964037);
void h_28(double *state, double *unused, double *out_6231938022780803683);
void H_28(double *state, double *unused, double *out_435617041443333213);
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
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
