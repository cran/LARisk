#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>


//--- Functions Declaration --------------------------------//
void call_lifeTable(double * dat_Prob_d_m, double * dat_Prob_d_f);
void call_incidTable(double * dat_Rate_m, double * dat_Rate_f);

double doseft(int * dist, double * info);
double latencyft(int * site);
double incidenceft(int age, int* sex, int* site);
double srateft(int e_age, int a_age, int* sex);
double DDREFft(int* DDREF_op, int* exposure_rate, double dose);

double approxft(double a, double b, double c);
double rtriangle(double a, double b, double c);
double vmvft(double* vec, double(*mat)[5], int column);

//--- Structure declaration for life & incidence table -----//
struct lifeTable
{
  double Prob_d_m;
  double Prob_d_f;
};
struct incidTable
{
  double Rate_m;
  double Rate_f;
};

//---Global Variables Declaration---------------------------//
struct lifeTable lifeTable[101];
struct incidTable incidTable[1919];

//--------- Main function-----------------------------------//
void larft(int * sex, int * age, int * exposeage, int * site, int * exposure_rate, int * dosedist, double * doseinfo, double * weight_err,
           int * DDREF_op, double * Prob_d_m, double * Prob_d_f, double * Rate_m, double * Rate_f,
           int * sim, double * ci, double * LAR_total, double * F_LAR_total, double * leukemia_result){
  // base data
  call_lifeTable(Prob_d_m, Prob_d_f);
  call_incidTable(Rate_m, Rate_f);

  // Variable Declaration
  int sum_el, len;
  double dose, latency;

  // Simulation
  if(*site==19){ // leukemia case
    // Initialization
    double SD_LAR_total[*sim];
    double F_SD_LAR_total[*sim];

    GetRNGstate();
    double z_value = qnorm(0.5+(*ci/2.0), 0.0, 1.0, 1, 0);
    PutRNGstate();


    // Parameters
    //double sigma_ERR_m[5][5] = {	{0.349923, -0.371562, 0.0279017, 0.0277847, 0.0206183},	{-0.371562, 0.542042, 0.0101649, 0.00150142, -0.00642578},	{0.0279017, 0.0101649, 0.0373066, 0.0169665, 0.0134684}, {0.0277847, 0.00150142, 0.0169665, 0.110026, 0.0578845}, {0.0206183, -0.00642578, 0.0134684, 0.0578845, 0.0696666}	};
    double sigma_EAR_m[5][5] = {	{0.677437, -0.545608, 0.0296734, 0, -0.00472918},	{-0.545608, 0.542886, 0.00169235, 0, 0.00737659},	{0.0296734, 0.00169235, 0.0217357, 0, 0.0103813},	{0, 0, 0, 0, 0},	{-0.00472918, 0.00737659, 0.0103813, 0, 0.0158421}		};
    //double sigma_ERR_f[5][5] = {	{0.420773, -0.399592, 0.0270153, 0.0244283, 0.0188047},	{-0.399592, 0.542042, 0.0101649, 0.00150142, -0.00642578},	{0.0270153, 0.0101649, 0.0373066, 0.0169665, 0.0134684},	{0.0244283, 0.00150142, 0.0169665, 0.110026, 0.0578845},	{0.0188047, -0.00642578, 0.0134684, 0.0578845, 0.0696666}	};
    double sigma_EAR_f[5][5] = {	{0.216888, -0.306358, 0.0150695, 0, -0.00201906},	{-0.306358, 0.542886, 0.00169235, 0, 0.00737659},	{0.0150695, 0.00169235, 0.0217357, 0, 0.0103813},	{0, 0, 0, 0, 0},		{-0.00201906, 0.00737659, 0.0103813, 0, 0.0158421}		};

    double gamma_ERR = -0.4, delta_ERR = -0.48, phi_ERR = 0.42, theta_ERR = 0.87;
    double gamma_EAR = 0.29, delta_EAR =  0.00, phi_EAR = 0.56, theta_EAR = 0.88;

    double beta_ERR, beta_EAR;

    if(*sex==1){ // male
      beta_ERR = 1.1, beta_EAR = 1.62;
    }else{ //female
      beta_ERR = 1.2, beta_EAR = 0.93;
    }

    if(*exposure_rate==2){ // chronic
      theta_ERR=0.0, theta_EAR=0.0;
    }

    // Simulation
    for(int s=0; s<*sim; s++){
      dose = doseft(dosedist, doseinfo); // exposure dose

      latency = latencyft(site); // latency
      sum_el = (int)(*exposeage + latency); // applying latency in exposed age
      len = (100 - sum_el) + 1; // length of sum...

      //compute ERR & EAR
      double * LAR_ERR = (double *)malloc(sizeof(double)* len);			// dynamic allocation
      double * LAR_EAR = (double *)malloc(sizeof(double)* len);			// dynamic allocation
      double * numerator_ERR = (double *)malloc(sizeof(double)* len);			// dynamic allocation
      double * numerator_EAR = (double *)malloc(sizeof(double)* len);			// dynamic allocation
      double * denominator_ERR = (double *)malloc(sizeof(double)* len);			// dynamic allocation
      double * denominator_EAR = (double *)malloc(sizeof(double)* len);			// dynamic allocation

      double estar=0.0, log_cal=0.0, incrate=0.0, srate=0.0, S_LAR_ERR=0.0, S_LAR_EAR=0.0;
      double S_numerator_ERR=0.0, S_numerator_EAR=0.0, S_denominator_ERR=0.0, S_denominator_EAR=0.0; // for integral (sum)...

      if(*exposeage < 30){
        estar = (*exposeage - (double)30) / (double)10;
      }else{
        estar = 0.0;
      }

      for(int a=sum_el; a<101; a++){
        log_cal = log((a - *exposeage) / (double)25);

        incrate = incidenceft(a, sex, site);
        srate = srateft(sum_el, a, sex);

        LAR_ERR[a - sum_el] = incrate * beta_ERR*dose * (1 + theta_ERR*dose) * exp(gamma_ERR*estar + delta_ERR*log_cal + phi_ERR*estar*log_cal) * srate;
        LAR_EAR[a - sum_el] = (double)10 * beta_EAR*dose * (1 + theta_EAR*dose) * exp(gamma_EAR*estar + delta_EAR*log_cal + phi_EAR*estar*log_cal) * srate;

        numerator_ERR[a - sum_el] = exp((delta_ERR + phi_ERR*estar) * log_cal) * incrate * srate * log_cal;
        numerator_EAR[a - sum_el] = exp((delta_EAR + phi_EAR*estar) * log_cal) * srate * log_cal;

        denominator_ERR[a - sum_el] = exp((delta_ERR + phi_ERR*estar) * log_cal) * incrate * srate;
        denominator_EAR[a - sum_el] = exp((delta_EAR + phi_EAR*estar) * log_cal) * srate;

        S_LAR_ERR += LAR_ERR[a - sum_el];
        S_LAR_EAR += LAR_EAR[a - sum_el];
        S_numerator_ERR += numerator_ERR[a - sum_el];
        S_numerator_EAR += numerator_EAR[a - sum_el];
        S_denominator_ERR += denominator_ERR[a - sum_el];
        S_denominator_EAR += denominator_EAR[a - sum_el];
      } // end integral (sum)...

      // Compute LAR
      LAR_total[s] = pow(S_LAR_ERR, *weight_err) * pow(S_LAR_EAR, (1 - *weight_err));

      // Compute Future LAR
      double F_S_LAR_ERR=0.0, F_S_LAR_EAR=0.0;
      double F_S_numerator_ERR=0.0, F_S_numerator_EAR=0.0, F_S_denominator_ERR=0.0, F_S_denominator_EAR=0.0; // for integral (sum)...

      if(*age<=sum_el){
        F_LAR_total[s] = LAR_total[s];
      }else{
        for(int a=*age; a<101; a++){
          log_cal = log((a - *exposeage) / (double)25);

          incrate = incidenceft(a, sex, site);
          srate = srateft(*age, a, sex);

          F_S_LAR_ERR += incrate * beta_ERR*dose * (1 + theta_ERR*dose) * exp(gamma_ERR*estar + delta_ERR*log_cal + phi_ERR*estar*log_cal) * srate;
          F_S_LAR_EAR += (double)10 * beta_EAR*dose * (1 + theta_EAR*dose) * exp(gamma_EAR*estar + delta_EAR*log_cal + phi_EAR*estar*log_cal) * srate;
          F_S_numerator_ERR += exp((delta_ERR + phi_ERR*estar) * log_cal) * incrate * srate * log_cal;
          F_S_numerator_EAR += exp((delta_EAR + phi_EAR*estar) * log_cal) * srate * log_cal;
          F_S_denominator_ERR += exp((delta_ERR + phi_ERR*estar) * log_cal) * incrate * srate;
          F_S_denominator_EAR += exp((delta_EAR + phi_EAR*estar) * log_cal) * srate;
        } // end integral (sum)...
        F_LAR_total[s] = pow(F_S_LAR_ERR, *weight_err) * pow(F_S_LAR_EAR, (1 - *weight_err));
      }

      // Compute variance of log(LAR)
      double S_nd_EAR=0, var_log_EAR=0; //S_nd_ERR=0, var_log_ERR
      //S_nd_ERR = S_numerator_ERR / S_denominator_ERR;
      S_nd_EAR = S_numerator_EAR / S_denominator_EAR;

      //double atrans_ERR[5] = { (double)1 / beta_ERR, dose / ((double)1 + theta_ERR * dose), estar, S_nd_ERR, S_nd_ERR*estar };
      double atrans_EAR[5] = { (double)1 / beta_EAR, dose / ((double)1 + theta_EAR * dose), estar, S_nd_EAR, S_nd_EAR*estar };

      if(*sex==1){
        //var_log_ERR = vmvft(atrans_ERR, sigma_ERR_m, sizeof(sigma_ERR_m)/sizeof(sigma_ERR_m[0]));
        var_log_EAR = vmvft(atrans_EAR, sigma_EAR_m, sizeof(sigma_EAR_m)/sizeof(sigma_EAR_m[0]));
      }else{
        //var_log_ERR = vmvft(atrans_ERR, sigma_ERR_f, sizeof(sigma_ERR_m)/sizeof(sigma_ERR_m[0]));
        var_log_EAR = vmvft(atrans_EAR, sigma_EAR_f, sizeof(sigma_EAR_m)/sizeof(sigma_EAR_m[0]));
      }

      SD_LAR_total[s] = sqrt( var_log_EAR + ( pow(log(S_LAR_ERR / S_LAR_EAR),2) * (*weight_err) * ((double)1 - *weight_err) ) );

      // Compute variance of Future log(LAR)
      if(*age<=sum_el){
        F_SD_LAR_total[s] = SD_LAR_total[s];
      }else{
        double F_S_nd_EAR=0, F_var_log_EAR=0; //F_S_nd_ERR=0, F_var_log_ERR
        //S_nd_ERR = F_S_numerator_ERR / F_S_denominator_ERR;
        S_nd_EAR = F_S_numerator_EAR / F_S_denominator_EAR;

        //double F_atrans_ERR[5] = { (double)1 / beta_ERR, dose / ((double)1 + theta_ERR * dose), estar, F_S_nd_ERR, F_S_nd_ERR*estar };
        double F_atrans_EAR[5] = { (double)1 / beta_EAR, dose / ((double)1 + theta_EAR * dose), estar, F_S_nd_EAR, F_S_nd_EAR*estar };

        if(*sex==1){
          //F_var_log_ERR = vmvft(F_atrans_ERR, sigma_ERR_m, sizeof(sigma_ERR_m)/sizeof(sigma_ERR_m[0]));
          F_var_log_EAR = vmvft(F_atrans_EAR, sigma_EAR_m, sizeof(sigma_EAR_m)/sizeof(sigma_EAR_m[0]));
        }else{
          //var_log_ERR = vmvft(F_atrans_ERR, sigma_ERR_f, sizeof(sigma_ERR_m)/sizeof(sigma_ERR_m[0]));
          var_log_EAR = vmvft(F_atrans_EAR, sigma_EAR_f, sizeof(sigma_EAR_m)/sizeof(sigma_EAR_m[0]));
        }

        F_SD_LAR_total[s] = sqrt( F_var_log_EAR + ( pow(log(F_S_LAR_ERR / F_S_LAR_EAR),2) * (*weight_err) * ((double)1 - *weight_err) ) );
      }

      free(LAR_ERR);
      free(LAR_EAR);
      free(numerator_ERR);
      free(numerator_EAR);
      free(denominator_ERR);
      free(denominator_EAR);

      leukemia_result[0] += exp( log(LAR_total[s]) - z_value * SD_LAR_total[s] );
      leukemia_result[1] += LAR_total[s];
      leukemia_result[2] += exp( log(LAR_total[s]) + z_value * SD_LAR_total[s] );

      leukemia_result[3] += exp( log(F_LAR_total[s]) - z_value * F_SD_LAR_total[s] );
      leukemia_result[4] += F_LAR_total[s];
      leukemia_result[5] += exp( log(F_LAR_total[s]) + z_value * F_SD_LAR_total[s] );
    }// end simulation

    if((*dosedist==1)&(*doseinfo==0.0)){
      leukemia_result[0] = (double) 0;
      leukemia_result[1] = (double) 0;
      leukemia_result[2] = (double) 0;

      leukemia_result[3] = (double) 0;
      leukemia_result[4] = (double) 0;
      leukemia_result[5] = (double) 0;

    }else{
      leukemia_result[0] = leukemia_result[0] / (double) *sim;
      leukemia_result[1] = leukemia_result[1] / (double) *sim;
      leukemia_result[2] = leukemia_result[2] / (double) *sim;

      leukemia_result[3] = leukemia_result[3] / (double) *sim;
      leukemia_result[4] = leukemia_result[4] / (double) *sim;
      leukemia_result[5] = leukemia_result[5] / (double) *sim;
    }



  }else{ // solid cancer case
    // initialization
    double gamma_ERR, eta_ERR, gamma_EAR, eta_EAR;
    double beta_ERR=0.0, beta_EAR=0.0;
    double DDREF;

    // parameters
    switch(*site){
    case 3: // liver
      gamma_ERR = -0.30; eta_ERR = -1.4;
      gamma_EAR = -0.41; eta_EAR =  4.1;
      break;
    case 4: // lung
      gamma_ERR = -0.30; eta_ERR = -1.4;
      gamma_EAR = -0.41; eta_EAR =  5.2;
      break;
    case 5: // breast
      gamma_ERR =  0.00; eta_ERR =  0.0;
      gamma_EAR = -0.50; eta_EAR =  3.5; // eta_EAR change
      break;
    case 9: // bladder
      gamma_ERR = -0.30; eta_ERR = -1.4;
      gamma_EAR = -0.41; eta_EAR =  6.0;
      break;
    case 10: // brain/cns
      gamma_ERR = -0.30; eta_ERR = -1.4;
      gamma_EAR =  0.00; eta_EAR =  0.0;
      break;
    case 11: // thyroid
      gamma_ERR = -0.83; eta_ERR =  0.0;
      gamma_EAR =  0.00; eta_EAR =  0.0;
      break;
    case 13: // oral
      gamma_ERR = -0.30; eta_ERR = -1.4;
      gamma_EAR = -0.41; eta_EAR =  0.5;
      break;
    case 16: // gallbladder
      gamma_ERR = -0.30; eta_ERR = -1.4;
      gamma_EAR =  0.00; eta_EAR =  0.0;
      break;
    default:
      gamma_ERR = -0.30; eta_ERR = -1.4;
      gamma_EAR = -0.41; eta_EAR =  2.8;
    }

    // Simulation
    for(int s=0; s<*sim; s++){
      dose = doseft(dosedist, doseinfo); // exposure dose

      latency = latencyft(site); // latency
      sum_el = (int)(*exposeage + latency); // applying latency in exposed age
      len = (100 - sum_el) + 1; // length of sum...

      switch(*site){
      case 1: //stomach
        if(*sex==1){
          GetRNGstate();
          beta_ERR = exp(rnorm(-1.5606477, 0.3293327));
          beta_EAR = exp(rnorm( 1.5892352, 0.3042856));
          PutRNGstate();
        }else{
          GetRNGstate();
          beta_ERR = exp(rnorm(-0.73396918, 0.21848782));
          beta_EAR = exp(rnorm( 1.58923521, 0.21038870));
          PutRNGstate();
        }
        break;
      case 2: // colon
        if(*sex==1){
          GetRNGstate();
          beta_ERR = exp(rnorm(-0.4620355, 0.2779496));
          beta_EAR = exp(rnorm( 1.1631508, 0.2895357));
          PutRNGstate();
        }else{
          GetRNGstate();
          beta_ERR = exp(rnorm( -0.84397007, 0.41324215));
          beta_EAR = exp(rnorm(  0.47000363, 0.35364650));
          PutRNGstate();
        }
        break;
      case 3: // liver
        if(*sex==1){
          GetRNGstate();
          beta_ERR = exp(rnorm(-1.1394343, 0.3536465));
          beta_EAR = exp(rnorm( 0.7884574, 0.4523131));
          PutRNGstate();
        }else{
          GetRNGstate();
          beta_ERR = exp(rnorm(-1.13943428, 0.58739416));
          beta_EAR = exp(rnorm( 0.00000000, 0.46749530));
          PutRNGstate();
        }
        break;
      case 4: // lung
        if(*sex==1){
          GetRNGstate();
          beta_ERR = exp(rnorm(-1.1394343, 0.3929707));
          beta_EAR = exp(rnorm( 0.8329091, 0.3862571));
          PutRNGstate();
        }else{
          GetRNGstate();
          beta_ERR = exp(rnorm(0.33647224, 0.20505427));
          beta_EAR = exp(rnorm(1.22377543, 0.19294030));
          PutRNGstate();
        }
        break;
      case 5: // breast
        beta_ERR = 0.0;
        GetRNGstate();
        beta_EAR = exp(rnorm(2.30258509, 0.1804418));
        PutRNGstate();
        break;
      case 6: // ovary
        GetRNGstate();
        beta_ERR = exp(rnorm(-0.96758403, 0.67322891));
        beta_EAR = exp(rnorm(-0.35667494, 0.59984060));
        PutRNGstate();
        break;
      case 7: // uterus
        GetRNGstate();
        beta_ERR = rnorm(0.055, 0.08418367);
        beta_EAR = rnorm(1.200, 0.71428570);
        PutRNGstate();
        break;
      case 8: // prostate
        GetRNGstate();
        beta_ERR = rnorm(0.12, 0.2908163);
        beta_EAR = rnorm(0.11, 0.4540816);
        PutRNGstate();
        break;
      case 9: // bladder
        if(*sex==1){
          GetRNGstate();
          beta_ERR = exp(rnorm(-0.6931472, 0.5232833));
          beta_EAR = exp(rnorm( 0.1823216, 0.5675060));
          PutRNGstate();
        }else{
          GetRNGstate();
          beta_ERR = exp(rnorm( 0.50077529, 0.44830562));
          beta_EAR = exp(rnorm(-0.28768207, 0.44250030));
          PutRNGstate();
        }
        break;
      case 10: // brain/cns
        if(*sex==1){
          GetRNGstate();
          beta_ERR = exp(rnorm(-0.3424903, 0.4183019));
          PutRNGstate();
          beta_EAR = 0.0;
        }else{
          GetRNGstate();
          beta_ERR = exp(rnorm(-1.42711636, 0.42166404));
          PutRNGstate();
          beta_EAR = 0.0;
        }
        break;
      case 11: // thyroid
        if(*sex==1){
          GetRNGstate();
          beta_ERR = exp(rnorm(-0.6348783, 0.6783827));
          PutRNGstate();
          beta_EAR = 0.0;
        }else{
          GetRNGstate();
          beta_ERR = exp(rnorm(0.04879016, 0.67192404));
          PutRNGstate();
          beta_EAR = 0.0;
        }
        break;
      case 12: // remainder
        if(*sex==1){
          GetRNGstate();
          beta_ERR = exp(rnorm(-0.1392621, 0.3375603));
          beta_EAR = exp(rnorm( 1.0043016, 0.2894181));
          PutRNGstate();
        }else{
          GetRNGstate();
          beta_ERR = exp(rnorm(-0.22314355, 0.45055679));
          beta_EAR = exp(rnorm(-0.05826891, 0.39222580));
          PutRNGstate();
        }
        break;
      case 13: // oral
        if(*sex==1){
          beta_ERR = approxft(-0.0047, 0.23, 0.66);
          beta_EAR = approxft(0.08, 0.44, 1.09);
        }else{
          beta_ERR = approxft(0.13, 0.53, 1.24);
          beta_EAR = approxft(0.06, 0.29, 0.66);
        }
        break;
      case 14: // oesophagus
        if(*sex==1){
          beta_ERR = approxft(0.093, 0.51, 1.13);
          beta_EAR = approxft(0.11, 0.88, 2.10);
        }else{
          beta_ERR = approxft(-0.10, 0.82, 3.09);
          beta_EAR = approxft(-0.08, 0.14, 0.63);
        }
        break;
      case 15: // rectum
        beta_ERR = approxft(-0.038, 0.12, 0.38);
        beta_EAR = approxft(-0.104, 0.34, 1.09);
        break;
      case 16: // gallbladder
        beta_ERR = approxft(-0.39, -0.018, 0.29);
        beta_EAR = 0.0;
        break;
      case 17: // pancreas
        beta_ERR = approxft(-0.0055, 0.36, 0.88);
        beta_EAR = approxft( 0.0904, 0.47, 1.08);
        break;
      case 18: // kidney
        beta_ERR = approxft(-0.012, 0.34, 1.00);
        beta_EAR = approxft( 0.082, 0.31, 0.68);
        break;
      } // parrameter generation

      //compute ERR & EAR
      double * LAR_ERR = (double *)malloc(sizeof(double)* len);			// dynamic allocation
      double * LAR_EAR = (double *)malloc(sizeof(double)* len);			// dynamic allocation

      double estar=0.0, astar=0.0, incrate=0.0, srate=0.0, S_LAR_ERR=0.0, S_LAR_EAR=0.0; // for integral (sum)...
      for(int a=sum_el; a<101; a++){

        if(*site==5){ //breast
          estar = (*exposeage - (double)25) / (double)10;
          astar = a / (double)50;

          if(a<=50){
            eta_EAR = 3.5;
          }else{
            eta_EAR = 1.0;
          }
        }else{ // others
          if(*exposeage < 30){
            estar = (*exposeage - (double)30) / (double)10;
          }else{
            estar = 0.0;
          }
          astar = a / (double)60;
        }

        incrate = incidenceft(a, sex, site);
        srate = srateft(sum_el, a, sex);

        LAR_ERR[a - sum_el] = incrate * beta_ERR*dose * exp(gamma_ERR*estar) * pow(astar, eta_ERR) * srate;
        LAR_EAR[a - sum_el] = (double)10 * beta_EAR*dose * exp(gamma_EAR*estar) * pow(astar, eta_EAR) * srate;

        S_LAR_ERR += LAR_ERR[a - sum_el];
        S_LAR_EAR += LAR_EAR[a - sum_el];
      } // end integral (sum)...

      //  compute LAR
      DDREF = DDREFft(DDREF_op, exposure_rate, dose);
      LAR_total[s] = (S_LAR_ERR * (*weight_err) + S_LAR_EAR * ((double)1 - *weight_err)) / DDREF;

      // compute Future LAR
      double F_S_LAR_ERR = 0.0, F_S_LAR_EAR = 0.0;

      if(*age<sum_el){
        F_LAR_total[s] = LAR_total[s];
      }else{
        for(int a=*age; a<101; a++){

          if(*site==5){ //breast
            estar = (*exposeage - (double)25) / (double)10;
            astar = a / (double)50;

            if(a<=50){
              eta_EAR = 3.5;
            }else{
              eta_EAR = 1.0;
            }
          }else{ // others
            if(*exposeage < 30){
              estar = (*exposeage - (double)30) / (double)10;
            }else{
              estar = 0.0;
            }
            astar = a / (double)60;
          }

          incrate = incidenceft(a, sex, site);
          srate = srateft(*age, a, sex);

          F_S_LAR_ERR += incrate * beta_ERR*dose * exp(gamma_ERR*estar) * pow(astar, eta_ERR) * srate;
          F_S_LAR_EAR += (double)10 * beta_EAR*dose * exp(gamma_EAR*estar) * pow(astar, eta_EAR) * srate;
        } // end integral (sum)...
        F_LAR_total[s] = (F_S_LAR_ERR* (*weight_err) + F_S_LAR_EAR*((double)1 - *weight_err)) / DDREF;
      }
    } // end simulation
  } // end solid cancer
}

void brft(int * sex, int * site, int * age, double * Prob_d_m, double * Prob_d_f, double * Rate_m, double * Rate_f, double * BR){
  // base data
  call_lifeTable(Prob_d_m, Prob_d_f);
  call_incidTable(Rate_m, Rate_f);

  if(*age>=100) *age=100;

  for(int a=*age; a<101; a++){
    *BR += srateft(*age, a, sex) * incidenceft(a, sex, site);
  }
}

//--------- Sub function -----------------------------------//
void call_lifeTable(double * dat_Prob_d_m, double * dat_Prob_d_f){
  // function to import NEW life table
  for(int i=0; i<101; i++){
    lifeTable[i].Prob_d_m = dat_Prob_d_m[i];
    lifeTable[i].Prob_d_f = dat_Prob_d_f[i];
  }
}

void call_incidTable(double * dat_Rate_m, double * dat_Rate_f){
  // function to import NEW incidence table
  for(int i=0; i<1919; i++){
    incidTable[i].Rate_m = dat_Rate_m[i];
    incidTable[i].Rate_f = dat_Rate_f[i];
  }
}


double doseft(int * dist, double * info){
  // function to calculate dose
  double doseresult=0.0;

  switch(*dist){
  case 1: // fixedvalue
    doseresult = info[0] / 1000.0;
    break;
  case 2: // lognormal
    GetRNGstate();
    doseresult = rlnorm(log(info[0]), log(info[1])) / 1000.0;
    PutRNGstate();
    break;
  case 3: // normal
    GetRNGstate();
    doseresult = rnorm(info[0], info[1]) / 1000.0;
    PutRNGstate();
    break;
  case 4: // triangular
    doseresult = rtriangle(info[0], info[1], info[2]) / 1000.0;
    break;
  case 5: // logtriangular
    doseresult = exp(rtriangle(log(info[0]), log(info[1]), log(info[2]))) / 1000.0;
    break;
  case 6: // uniform
    GetRNGstate();
    doseresult = runif(info[0],info[1]) / 1000.0;
    PutRNGstate();
    break;
  case 7: // loguniform
    GetRNGstate();
    doseresult = exp(runif(log(info[0]),log(info[1]))) / 1000.0;
    PutRNGstate();
    break;
  }

  return doseresult;
}

double latencyft(int * site){
  double latencyresult;

  switch(*site){
  case 19:
    while (1) {
      GetRNGstate();
      latencyresult = rtriangle(2.0, 2.25, 2.5) - 0.401*log((1.0 / runif(0.0,1.0)) - 1.0);
      PutRNGstate();
      if ((0.4 <= latencyresult) & (latencyresult <= 4.1)) break;
    }
  case 11:
    while (1) {
      GetRNGstate();
      latencyresult = rtriangle(2.0, 5.00, 7.0) - 0.544*log((1.0 / runif(0.0,1.0)) - 1.0);
      PutRNGstate();
      if ((2.5 <= latencyresult) & (latencyresult <= 7.6)) break;
    }
  default:
    while (1) {
      GetRNGstate();
      latencyresult = rtriangle(5.0, 7.5, 10.0) - 0.76*log((1.0 / runif(0.0,1.0)) - 1.0);
      PutRNGstate();
      if ((4.0 <= latencyresult) & (latencyresult <= 11.0)) break;
    }
  }

  return latencyresult;
}

double incidenceft(int age, int* sex, int* site){
  double result = 1.0;

  if(*sex==1){
    result = incidTable[(*site - 1) * 101 + (age)].Rate_m;
  }else{
    result = incidTable[(*site - 1) * 101 + (age)].Rate_f;
  }

  return result;
}

double srateft(int e_age, int a_age, int* sex){
  double product = 1.0;

  if(*sex==1){
    for(int i=e_age+1; i<a_age; i++){
      product *= 1.0 - lifeTable[i].Prob_d_m;
    }
  }else{
    for(int i=e_age+1; i<a_age; i++){
      product *= 1.0 - lifeTable[i].Prob_d_f;
    }
  }

  return product;
}

double DDREFft(int* DDREF_op, int* exposure_rate, double dose){
  double DDREF_tmp, Dlim;

  if(*DDREF_op==1){
    if(*exposure_rate==1){ //acute
      GetRNGstate();
      Dlim = exp(runif(log(0.03), log(0.2)));
      if(dose<Dlim){
        DDREF_tmp = exp(rnorm(log(1.5), log(1.35)));
      }else{
        DDREF_tmp = 1.0;
      }
      PutRNGstate();
    }else{ //chronic
      GetRNGstate();
      DDREF_tmp = exp(rnorm(log(1.5), log(1.35)));
      PutRNGstate();
    }
  }else{
    DDREF_tmp = 1.0;
  }
  return DDREF_tmp;
}

double approxft(double a, double b, double c){
  double result = 0.0;
  double r_unif = 0.0;

  GetRNGstate();
  r_unif = runif(0.0001, 0.9999);
  PutRNGstate();

  if (r_unif <= 0.5) {
    result = ((b-a)*(r_unif - 0.025) / 0.475) + a;
  }
  else {
    result = ((c-b)*(r_unif - 0.5) / 0.475) + b;
  }

  return result;
}

double rtriangle(double a, double b, double c){
  double F = (b - a) / (c - a);
  double r_unif;
  double result;

  GetRNGstate();
  r_unif = runif(0.0, 1.0);
  PutRNGstate();

  if ((0 < r_unif) & (r_unif < F)) {
    result = a + sqrt(r_unif*(c - a)*(b - a));
  }
  else {
    result = b - sqrt((1 - r_unif)*(c - a)*(c - b));
  }
  return result;
}

double vmvft(double* vec, double(*mat)[5], int column){
  double prod_tmp1, prod_tmp2=0;
  for(int a=0; a<5; a++){
    prod_tmp1=0;
    for(int b=0; b<5; b++){
      prod_tmp1 += vec[b] * mat[a][b];
    }
    prod_tmp2 += prod_tmp1 * vec[a];
  }

  return prod_tmp2;
}

// RegisteringDynamic Symbols
void R_init_LARisk(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
