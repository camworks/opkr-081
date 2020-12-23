/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_1705903580326446379);
void inv_err_fun(double *nom_x, double *true_x, double *out_1978619372336320749);
void H_mod_fun(double *state, double *out_9009401118396320287);
void f_fun(double *state, double dt, double *out_9182124003117273804);
void F_fun(double *state, double dt, double *out_6182448584763505904);
void h_3(double *state, double *unused, double *out_8858959669317795453);
void H_3(double *state, double *unused, double *out_6475888940787819655);
void h_4(double *state, double *unused, double *out_1593994730339534882);
void H_4(double *state, double *unused, double *out_959860441889961338);
void h_9(double *state, double *unused, double *out_4494595797978160276);
void H_9(double *state, double *unused, double *out_5761002484965262760);
void h_10(double *state, double *unused, double *out_2343020947428809476);
void H_10(double *state, double *unused, double *out_969346734059019354);
void h_12(double *state, double *unused, double *out_745762473140783671);
void H_12(double *state, double *unused, double *out_5545052452833616359);
void h_31(double *state, double *unused, double *out_7737512320488197638);
void H_31(double *state, double *unused, double *out_2869962677690196538);
void h_32(double *state, double *unused, double *out_4528292194023430765);
void H_32(double *state, double *unused, double *out_3430403066570204648);
void h_13(double *state, double *unused, double *out_5742034062832093046);
void H_13(double *state, double *unused, double *out_5953216333866922752);
void h_14(double *state, double *unused, double *out_4494595797978160276);
void H_14(double *state, double *unused, double *out_5761002484965262760);
void h_19(double *state, double *unused, double *out_2599082170334284050);
void H_19(double *state, double *unused, double *out_8253061737641857675);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);