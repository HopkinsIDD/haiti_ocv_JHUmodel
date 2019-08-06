## pomp model in C for the national vaccination campaign

## ---------------------------------------- ##
## Step function used in simulating process 

rSim.step <- Csnippet('
  // transition rates
  double Srate[12]; 
  double Erate[13];
  double Irate[12];
  double Arate[12];
  double Rrate[12];
  double S1rate[2];
  double E1rate[3];
  double I1rate[2];
  double A1rate[2];
  double R1rate[2];
  double S2rate[2];
  double E2rate[3];
  double I2rate[2];
  double A2rate[2];
  double R2rate[2];
  double S3rate[2];
  double E3rate[3];
  double I3rate[2];
  double A3rate[2];
  double R3rate[2];
  double S4rate[2];
  double E4rate[3];
  double I4rate[2];
  double A4rate[2];
  double R4rate[2];
  double S5rate[2];
  double E5rate[3];
  double I5rate[2];
  double A5rate[2];
  double R5rate[2];
  double S6rate[2];
  double E6rate[3];
  double I6rate[2];
  double A6rate[2];
  double R6rate[2];
  double S7rate[2];
  double E7rate[3];
  double I7rate[2];
  double A7rate[2];
  double R7rate[2];
  double S8rate[2];
  double E8rate[3];
  double I8rate[2];
  double A8rate[2];
  double R8rate[2];
  double S9rate[2];
  double E9rate[3];
  double I9rate[2];
  double A9rate[2];
  double R9rate[2];
  double S10rate[2];
  double E10rate[3];
  double I10rate[2];
  double A10rate[2];
  double R10rate[2];

  // transition numbers
  double Strans[12]; 
  double Etrans[13];
  double Itrans[12];
  double Atrans[12];
  double Rtrans[12];
  double S1trans[2];
  double E1trans[3];
  double I1trans[2];
  double A1trans[2];
  double R1trans[2];
  double S2trans[2];
  double E2trans[3];
  double I2trans[2];
  double A2trans[2];
  double R2trans[2];
  double S3trans[2];
  double E3trans[3];
  double I3trans[2];
  double A3trans[2];
  double R3trans[2];
  double S4trans[2];
  double E4trans[3];
  double I4trans[2];
  double A4trans[2];
  double R4trans[2];
  double S5trans[2];
  double E5trans[3];
  double I5trans[2];
  double A5trans[2];
  double R5trans[2];
  double S6trans[2];
  double E6trans[3];
  double I6trans[2];
  double A6trans[2];
  double R6trans[2];
  double S7trans[2];
  double E7trans[3];
  double I7trans[2];
  double A7trans[2];
  double R7trans[2];
  double S8trans[2];
  double E8trans[3];
  double I8trans[2];
  double A8trans[2];
  double R8trans[2];
  double S9trans[2];
  double E9trans[3];
  double I9trans[2];
  double A9trans[2];
  double R9trans[2];
  double S10trans[2];
  double E10trans[3];
  double I10trans[2];
  double A10trans[2];
  double R10trans[2];

  // vac numbers
  double vac1 = 0.0;
  double vac2 = 0.0;
  double vac3 = 0.0;
  double vac4 = 0.0;
  double vac5 = 0.0;
  double vac6 = 0.0;
  double vac7 = 0.0;
  double vac8 = 0.0;
  double vac9 = 0.0;
  double vac10 = 0.0;

  // some population demonitors
  int Nnv = S + E + I + A + R;
  int N1 = S1 + E1 + I1 + A1 + R1; // Centre
  int N2 = S2 + E2 + I2 + A2 + R2; // Artibonite
  int N3 = S3 + E3 + I3 + A3 + R3; // Ouest
  int N4 = S4 + E4 + I4 + A4 + R4; // Nord Ouest
  int N5 = S5 + E5 + I5 + A5 + R5; // Nord
  int N6 = S6 + E6 + I6 + A6 + R6; // Sud
  int N7 = S7 + E7 + I7 + A7 + R7; // Nippes
  int N8 = S8 + E8 + I8 + A8 + R8; // Nord Est
  int N9 = S9 + E9 + I9 + A9 + R9; // Sud Est
  int N10 = S10 + E10 + I10 + A10 + R10; // GrandAnse
  int N = Nnv + N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8 + N9 + N10;
  int births = rpois(mu*N*dt);

  // time checks for vaccination campaign (from covariate)
  if (vac_tcheck == 1) {
      vac1 = num_vacc/Nnv/dt;  
  } 
  if (vac_tcheck == 2){
      vac2 = num_vacc/Nnv/dt; 
  }
  if (vac_tcheck == 3){
      vac3 = num_vacc/Nnv/dt; 
  }
  if (vac_tcheck == 4) {
      vac4 = num_vacc/Nnv/dt;  
  } 
  if (vac_tcheck == 5){
      vac5 = num_vacc/Nnv/dt; 
  }
  if (vac_tcheck == 6){
      vac6 = num_vacc/Nnv/dt; 
  }
  if (vac_tcheck == 7) {
      vac7 = num_vacc/Nnv/dt;  
  } 
  if (vac_tcheck == 8){
      vac8 = num_vacc/Nnv/dt; 
  }
  if (vac_tcheck == 9){
      vac9 = num_vacc/Nnv/dt; 
  }
  if (vac_tcheck == 10) {
      vac10 = num_vacc/Nnv/dt;  
  } 

  // make seasonal beta term for current time
  double mybeta = beta1*seas1 + beta2*seas2 + beta3*seas3 + beta4*seas4 + beta5*seas5 + beta6*seas6;
  double foi = pow(I+I1+I2+I3+I4+I5+I6+I7+I8+I9+I10+(1-kappa)*(A+A1+A2+A3+A4+A5+A6+A7+A8+A9+A10), nu)*mybeta/N;  

  // make thetav (artificially from covariate dataframe)
  double thetav_d1 = ve_d1;
  double thetav_d2 = ve_d2;
  double thetav_d3 = ve_d3;
  double thetav_d4 = ve_d4;
  double thetav_d5 = ve_d5;
  double thetav_d6 = ve_d6;
  double thetav_d7 = ve_d7;
  double thetav_d8 = ve_d8;
  double thetav_d9 = ve_d9;
  double thetav_d10 = ve_d10;
 
  //compute the rates for all the transitions
  Srate[0]= foi;  //S -> E
  Srate[1]= delta;
  Srate[2]= vac1;
  Srate[3]= vac2;
  Srate[4]= vac3;
  Srate[5]= vac4;
  Srate[6]= vac5;
  Srate[7]= vac6;
  Srate[8]= vac7;
  Srate[9]= vac8;
  Srate[10]= vac9;
  Srate[11]= vac10;

  Erate[0]= sigma*(1-theta0); // E -> I
  Erate[1]= sigma*theta0; // E -> A
  Erate[2]= delta;
  Erate[3]= vac1;
  Erate[4]= vac2;
  Erate[5]= vac3;
  Erate[6]= vac4;
  Erate[7]= vac5;
  Erate[8]= vac6;
  Erate[9]= vac7;
  Erate[10]= vac8;
  Erate[11]= vac9;
  Erate[12]= vac10;

  Irate[0]= gamma; // I -> R
  Irate[1]= delta;
  Irate[2]= vac1;
  Irate[3]= vac2;
  Irate[4]= vac3;
  Irate[5]= vac4;
  Irate[6]= vac5;
  Irate[7]= vac6;
  Irate[8]= vac7;
  Irate[9]= vac8;
  Irate[10]= vac9;
  Irate[11]= vac10;

  Arate[0]= gamma; // A -> R
  Arate[1]= delta;
  Arate[2]= vac1;
  Arate[3]= vac2;
  Arate[4]= vac3;
  Arate[5]= vac4;
  Arate[6]= vac5;
  Arate[7]= vac6;
  Arate[8]= vac7;
  Arate[9]= vac8;
  Arate[10]= vac9;
  Arate[11]= vac10;

  Rrate[0]= alpha; // R -> S exponential waning immunity from natural infection
  Rrate[1]= delta;
  Rrate[2]= vac1;
  Rrate[3]= vac2;
  Rrate[4]= vac3;
  Rrate[5]= vac4;
  Rrate[6]= vac5;
  Rrate[7]= vac6;
  Rrate[8]= vac7;
  Rrate[9]= vac8;
  Rrate[10]= vac9;
  Rrate[11]= vac10;

  S1rate[0]= foi; 
  S1rate[1]= delta;
  E1rate[0]= sigma*(1-thetav_d1); 
  E1rate[1]= sigma*thetav_d1; 
  E1rate[2]= delta;
  I1rate[0]= gamma; 
  I1rate[1]= delta;
  A1rate[0]= gamma; 
  A1rate[1]= delta;
  R1rate[0]= alpha; 
  R1rate[1]= delta;

  S2rate[0]= foi; 
  S2rate[1]= delta;
  E2rate[0]= sigma*(1-thetav_d2); 
  E2rate[1]= sigma*thetav_d2; 
  E2rate[2]= delta;
  I2rate[0]= gamma; 
  I2rate[1]= delta;
  A2rate[0]= gamma; 
  A2rate[1]= delta;
  R2rate[0]= alpha; 
  R2rate[1]= delta;

  S3rate[0]= foi; 
  S3rate[1]= delta;
  E3rate[0]= sigma*(1-thetav_d3); 
  E3rate[1]= sigma*thetav_d3; 
  E3rate[2]= delta;
  I3rate[0]= gamma; 
  I3rate[1]= delta;
  A3rate[0]= gamma; 
  A3rate[1]= delta;
  R3rate[0]= alpha; 
  R3rate[1]= delta;

  S4rate[0]= foi; 
  S4rate[1]= delta;
  E4rate[0]= sigma*(1-thetav_d4); 
  E4rate[1]= sigma*thetav_d4; 
  E4rate[2]= delta;
  I4rate[0]= gamma; 
  I4rate[1]= delta;
  A4rate[0]= gamma; 
  A4rate[1]= delta;
  R4rate[0]= alpha; 
  R4rate[1]= delta;

  S5rate[0]= foi; 
  S5rate[1]= delta;
  E5rate[0]= sigma*(1-thetav_d5); 
  E5rate[1]= sigma*thetav_d5; 
  E5rate[2]= delta;
  I5rate[0]= gamma; 
  I5rate[1]= delta;
  A5rate[0]= gamma; 
  A5rate[1]= delta;
  R5rate[0]= alpha; 
  R5rate[1]= delta;

  S6rate[0]= foi; 
  S6rate[1]= delta;
  E6rate[0]= sigma*(1-thetav_d6); 
  E6rate[1]= sigma*thetav_d6; 
  E6rate[2]= delta;
  I6rate[0]= gamma; 
  I6rate[1]= delta;
  A6rate[0]= gamma; 
  A6rate[1]= delta;
  R6rate[0]= alpha; 
  R6rate[1]= delta;

  S7rate[0]= foi; 
  S7rate[1]= delta;
  E7rate[0]= sigma*(1-thetav_d7); 
  E7rate[1]= sigma*thetav_d7; 
  E7rate[2]= delta;
  I7rate[0]= gamma; 
  I7rate[1]= delta;
  A7rate[0]= gamma; 
  A7rate[1]= delta;
  R7rate[0]= alpha; 
  R7rate[1]= delta;

  S8rate[0]= foi; 
  S8rate[1]= delta;
  E8rate[0]= sigma*(1-thetav_d8); 
  E8rate[1]= sigma*thetav_d8; 
  E8rate[2]= delta;
  I8rate[0]= gamma; 
  I8rate[1]= delta;
  A8rate[0]= gamma; 
  A8rate[1]= delta;
  R8rate[0]= alpha; 
  R8rate[1]= delta;

  S9rate[0]= foi; 
  S9rate[1]= delta;
  E9rate[0]= sigma*(1-thetav_d9); 
  E9rate[1]= sigma*thetav_d9; 
  E9rate[2]= delta;
  I9rate[0]= gamma; 
  I9rate[1]= delta;
  A9rate[0]= gamma; 
  A9rate[1]= delta;
  R9rate[0]= alpha; 
  R9rate[1]= delta;

  S10rate[0]= foi; 
  S10rate[1]= delta;
  E10rate[0]= sigma*(1-thetav_d10); 
  E10rate[1]= sigma*thetav_d10; 
  E10rate[2]= delta;
  I10rate[0]= gamma; 
  I10rate[1]= delta;
  A10rate[0]= gamma; 
  A10rate[1]= delta;
  R10rate[0]= alpha; 
  R10rate[1]= delta;

  // compute the transition numbers
  reulermultinom(12,S,&Srate[0],dt,&Strans[0]);
  reulermultinom(13,E,&Erate[0],dt,&Etrans[0]);
  reulermultinom(12,I,&Irate[0],dt,&Itrans[0]);
  reulermultinom(12,A,&Arate[0],dt,&Atrans[0]);
  reulermultinom(12,R,&Rrate[0],dt,&Rtrans[0]);

  reulermultinom(2,S1,&S1rate[0],dt,&S1trans[0]);
  reulermultinom(3,E1,&E1rate[0],dt,&E1trans[0]);
  reulermultinom(2,I1,&I1rate[0],dt,&I1trans[0]);
  reulermultinom(2,A1,&A1rate[0],dt,&A1trans[0]);
  reulermultinom(2,R1,&R1rate[0],dt,&R1trans[0]);

  reulermultinom(2,S2,&S2rate[0],dt,&S2trans[0]);
  reulermultinom(3,E2,&E2rate[0],dt,&E2trans[0]);
  reulermultinom(2,I2,&I2rate[0],dt,&I2trans[0]);
  reulermultinom(2,A2,&A2rate[0],dt,&A2trans[0]);
  reulermultinom(2,R2,&R2rate[0],dt,&R2trans[0]);

  reulermultinom(2,S3,&S3rate[0],dt,&S3trans[0]);
  reulermultinom(3,E3,&E3rate[0],dt,&E3trans[0]);
  reulermultinom(2,I3,&I3rate[0],dt,&I3trans[0]);
  reulermultinom(2,A3,&A3rate[0],dt,&A3trans[0]);
  reulermultinom(2,R3,&R3rate[0],dt,&R3trans[0]);

  reulermultinom(2,S4,&S4rate[0],dt,&S4trans[0]);
  reulermultinom(3,E4,&E4rate[0],dt,&E4trans[0]);
  reulermultinom(2,I4,&I4rate[0],dt,&I4trans[0]);
  reulermultinom(2,A4,&A4rate[0],dt,&A4trans[0]);
  reulermultinom(2,R4,&R4rate[0],dt,&R4trans[0]);

  reulermultinom(2,S5,&S5rate[0],dt,&S5trans[0]);
  reulermultinom(3,E5,&E5rate[0],dt,&E5trans[0]);
  reulermultinom(2,I5,&I5rate[0],dt,&I5trans[0]);
  reulermultinom(2,A5,&A5rate[0],dt,&A5trans[0]);
  reulermultinom(2,R5,&R5rate[0],dt,&R5trans[0]);

  reulermultinom(2,S6,&S6rate[0],dt,&S6trans[0]);
  reulermultinom(3,E6,&E6rate[0],dt,&E6trans[0]);
  reulermultinom(2,I6,&I6rate[0],dt,&I6trans[0]);
  reulermultinom(2,A6,&A6rate[0],dt,&A6trans[0]);
  reulermultinom(2,R6,&R6rate[0],dt,&R6trans[0]);

  reulermultinom(2,S7,&S7rate[0],dt,&S7trans[0]);
  reulermultinom(3,E7,&E7rate[0],dt,&E7trans[0]);
  reulermultinom(2,I7,&I7rate[0],dt,&I7trans[0]);
  reulermultinom(2,A7,&A7rate[0],dt,&A7trans[0]);
  reulermultinom(2,R7,&R7rate[0],dt,&R7trans[0]);

  reulermultinom(2,S8,&S8rate[0],dt,&S8trans[0]);
  reulermultinom(3,E8,&E8rate[0],dt,&E8trans[0]);
  reulermultinom(2,I8,&I8rate[0],dt,&I8trans[0]);
  reulermultinom(2,A8,&A8rate[0],dt,&A8trans[0]);
  reulermultinom(2,R8,&R8rate[0],dt,&R8trans[0]);

  reulermultinom(2,S9,&S9rate[0],dt,&S9trans[0]);
  reulermultinom(3,E9,&E9rate[0],dt,&E9trans[0]);
  reulermultinom(2,I9,&I9rate[0],dt,&I9trans[0]);
  reulermultinom(2,A9,&A9rate[0],dt,&A9trans[0]);
  reulermultinom(2,R9,&R9rate[0],dt,&R9trans[0]);

  reulermultinom(2,S10,&S10rate[0],dt,&S10trans[0]);
  reulermultinom(3,E10,&E10rate[0],dt,&E10trans[0]);
  reulermultinom(2,I10,&I10rate[0],dt,&I10trans[0]);
  reulermultinom(2,A10,&A10rate[0],dt,&A10trans[0]);
  reulermultinom(2,R10,&R10rate[0],dt,&R10trans[0]);
  
  // balance the equations
  S += -Strans[0] - Strans[1] - Strans[2] - Strans[3] - Strans[4] - Strans[5] - Strans[6] - Strans[7] - Strans[8] - Strans[9] - Strans[10] - Strans[11] + Rtrans[0] + births;
  E += -Etrans[0] - Etrans[1] - Etrans[2] - Etrans[3] - Etrans[4] - Etrans[5] - Etrans[6] - Etrans[7] - Etrans[8] - Etrans[9] - Etrans[10] - Etrans[11] - Etrans[12] + Strans[0];
  I += -Itrans[0] - Itrans[1] - Itrans[2] - Itrans[3] - Itrans[4] - Itrans[5] - Itrans[6] - Itrans[7] - Itrans[8] - Itrans[9] - Itrans[10] - Itrans[11] + Etrans[0];
  A += -Atrans[0] - Atrans[1] - Atrans[2] - Atrans[3] - Atrans[4] - Atrans[5] - Atrans[6] - Atrans[7] - Atrans[8] - Atrans[9] - Atrans[10] - Atrans[11] + Etrans[1];
  R += -Rtrans[0] - Rtrans[1] - Rtrans[2] - Rtrans[3] - Rtrans[4] - Rtrans[5] - Rtrans[6] - Rtrans[7] - Rtrans[8] - Rtrans[9] - Rtrans[10] - Rtrans[11] + Itrans[0] + Atrans[0];

  S1 += -S1trans[0] - S1trans[1] + R1trans[0] + Strans[2];
  E1 += -E1trans[0] - E1trans[1] - E1trans[2] + S1trans[0] + Etrans[3];
  I1 += -I1trans[0] - I1trans[1] + E1trans[0] + Itrans[2];
  A1 += -A1trans[0] - A1trans[1] + E1trans[1] + Atrans[2];
  R1 += -R1trans[0] - R1trans[1] + I1trans[0] + A1trans[0] + Rtrans[2];

  S2 += -S2trans[0] - S2trans[1] + R2trans[0] + Strans[3];
  E2 += -E2trans[0] - E2trans[1] - E2trans[2] + S2trans[0] + Etrans[4];
  I2 += -I2trans[0] - I2trans[1] + E2trans[0] + Itrans[3];
  A2 += -A2trans[0] - A2trans[1] + E2trans[1] + Atrans[3];
  R2 += -R2trans[0] - R2trans[1] + I2trans[0] + A2trans[0] + Rtrans[3];

  S3 += -S3trans[0] - S3trans[1] + R3trans[0] + Strans[4];
  E3 += -E3trans[0] - E3trans[1] - E3trans[2] + S3trans[0] + Etrans[5];
  I3 += -I3trans[0] - I3trans[1] + E3trans[0] + Itrans[4];
  A3 += -A3trans[0] - A3trans[1] + E3trans[1] + Atrans[4];
  R3 += -R3trans[0] - R3trans[1] + I3trans[0] + A3trans[0] + Rtrans[4];

  S4 += -S4trans[0] - S4trans[1] + R4trans[0] + Strans[5];
  E4 += -E4trans[0] - E4trans[1] - E4trans[2] + S4trans[0] + Etrans[6];
  I4 += -I4trans[0] - I4trans[1] + E4trans[0] + Itrans[5];
  A4 += -A4trans[0] - A4trans[1] + E4trans[1] + Atrans[5];
  R4 += -R4trans[0] - R4trans[1] + I4trans[0] + A4trans[0] + Rtrans[5];

  S5 += -S5trans[0] - S5trans[1] + R5trans[0] + Strans[6];
  E5 += -E5trans[0] - E5trans[1] - E5trans[2] + S5trans[0] + Etrans[7];
  I5 += -I5trans[0] - I5trans[1] + E5trans[0] + Itrans[6];
  A5 += -A5trans[0] - A5trans[1] + E5trans[1] + Atrans[6];
  R5 += -R5trans[0] - R5trans[1] + I5trans[0] + A5trans[0] + Rtrans[6];

  S6 += -S6trans[0] - S6trans[1] + R6trans[0] + Strans[7];
  E6 += -E6trans[0] - E6trans[1] - E6trans[2] + S6trans[0] + Etrans[8];
  I6 += -I6trans[0] - I6trans[1] + E6trans[0] + Itrans[7];
  A6 += -A6trans[0] - A6trans[1] + E6trans[1] + Atrans[7];
  R6 += -R6trans[0] - R6trans[1] + I6trans[0] + A6trans[0] + Rtrans[7];

  S7 += -S7trans[0] - S7trans[1] + R7trans[0] + Strans[8];
  E7 += -E7trans[0] - E7trans[1] - E7trans[2] + S7trans[0] + Etrans[9];
  I7 += -I7trans[0] - I7trans[1] + E7trans[0] + Itrans[8];
  A7 += -A7trans[0] - A7trans[1] + E7trans[1] + Atrans[8];
  R7 += -R7trans[0] - R7trans[1] + I7trans[0] + A7trans[0] + Rtrans[8];

  S8 += -S8trans[0] - S8trans[1] + R8trans[0] + Strans[9];
  E8 += -E8trans[0] - E8trans[1] - E8trans[2] + S8trans[0] + Etrans[10];
  I8 += -I8trans[0] - I8trans[1] + E8trans[0] + Itrans[9];
  A8 += -A8trans[0] - A8trans[1] + E8trans[1] + Atrans[9];
  R8 += -R8trans[0] - R8trans[1] + I8trans[0] + A8trans[0] + Rtrans[9];

  S9 += -S9trans[0] - S9trans[1] + R9trans[0] + Strans[10];
  E9 += -E9trans[0] - E9trans[1] - E9trans[2] + S9trans[0] + Etrans[11];
  I9 += -I9trans[0] - I9trans[1] + E9trans[0] + Itrans[10];
  A9 += -A9trans[0] - A9trans[1] + E9trans[1] + Atrans[10];
  R9 += -R9trans[0] - R9trans[1] + I9trans[0] + A9trans[0] + Rtrans[10];

  S10 += -S10trans[0] - S10trans[1] + R10trans[0] + Strans[11];
  E10 += -E10trans[0] - E10trans[1] - E10trans[2] + S10trans[0] + Etrans[12];
  I10 += -I10trans[0] - I10trans[1] + E10trans[0] + Itrans[11];
  A10 += -A10trans[0] - A10trans[1] + E10trans[1] + Atrans[11];
  R10 += -R10trans[0] - R10trans[1] + I10trans[0] + A10trans[0] + Rtrans[11];

  incid += Etrans[0] + E1trans[0] + E2trans[0] + E3trans[0] + E4trans[0] + E5trans[0] + E6trans[0] + E7trans[0] + E8trans[0] + E9trans[0] + E10trans[0]; // incidence is cumulative entries into I state
  incidU += Etrans[0];
  incidV += E1trans[0] + E2trans[0] + E3trans[0] + E4trans[0] + E5trans[0] + E6trans[0] + E7trans[0] + E8trans[0] + E9trans[0] + E10trans[0];
  asymV += E1trans[1] + E2trans[1] + E3trans[1] + E4trans[1] + E5trans[1] + E6trans[1] + E7trans[1] + E8trans[1] + E9trans[1] + E10trans[1];
  newV += Strans[2] + Strans[3] + Strans[4] + Strans[5] + Strans[6] + Strans[7] + Strans[8] + Strans[9] + Strans[10] + Strans[11] + 
    Etrans[3] + Etrans[4] + Etrans[5] + Etrans[6] + Etrans[7] + Etrans[8] + Etrans[9] + Etrans[10] + Etrans[11] + Etrans[12] +
    Itrans[2] + Itrans[3] + Itrans[4] + Itrans[5] + Itrans[6] + Itrans[7] + Itrans[8] + Itrans[9] + Itrans[10] + Itrans[11] + 
    Atrans[2] + Atrans[3] + Atrans[4] + Atrans[5] + Atrans[6] + Atrans[7] + Atrans[8] + Atrans[9] + Atrans[10] + Atrans[11] + 
    Rtrans[2] + Rtrans[3] + Rtrans[4] + Rtrans[5] + Rtrans[6] + Rtrans[7] + Rtrans[8] + Rtrans[9] + Rtrans[10] + Rtrans[11];
  foival += foi;
  Str0 += Strans[0];
  Sout += Strans[0] + Strans[1] + Strans[2] + Strans[3] + Strans[4] + Strans[5] + Strans[6] + Strans[7] + Strans[8] + Strans[9] + Strans[10] + Strans[11];
  Sin += Rtrans[0] + births;

  //if (R_FINITE(incid)) 
  //  Rprintf(\" Srate %i  \\n\", -Strans[0] - Strans[1] + Rtrans[0] + births);
')
# pow((I+I1+I2+(1-kappa)*(A+A1+A2)), nu)*mybeta/N;

## ---------------------------------------- ##
## Simulation of the measurement process
rMeas <- Csnippet('
  cases = rnbinom_mu(theta, rho*incid);
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } else {
    cases = 0.0; // 0.0;
  }
')


## ---------------------------------------- ##
## Build pomp model

build.fc.mod <- function(pop, ##  = get.haiti.pop()
                          dat, ##  = get.mspp.agg.data()
                          my.times = "week", 
                          covar.times = "time",
                          my.t0 = 0,
                          covar){

  my.mod <- pomp(
    data = dat,
    times = my.times,
    t0 = my.t0,
    rmeasure = rMeas,
    rprocess = euler.sim(step.fun=rSim.step, delta.t=1/7),
    covar = covar,
    tcovar = covar.times, 
    statenames=c("S", sprintf("S%d", 1:10),
                 "E", sprintf("E%d", 1:10),
                 "I", sprintf("I%d", 1:10),
                 "A", sprintf("A%d", 1:10),
                 "R", sprintf("R%d", 1:10),
                 "incid","incidU","incidV","asymV","newV","foival","Str0","Sout","Sin"),
    paramnames=c(
        "rho","theta",
        "beta1","beta2","beta3",
        "beta4",
        "beta5",
        "beta6",
        "gamma","sigma","theta0","alpha","mu","delta","nu","kappa",
        "S.0","E.0","I.0","A.0","R.0",
        "S1.0","E1.0","I1.0","A1.0","R1.0",
        "S2.0","E2.0","I2.0","A2.0","R2.0",
        "S3.0","E3.0","I3.0","A3.0","R3.0",
        "S4.0","E4.0","I4.0","A4.0","R4.0",
        "S5.0","E5.0","I5.0","A5.0","R5.0",
        "S6.0","E6.0","I6.0","A6.0","R6.0",
        "S7.0","E7.0","I7.0","A7.0","R7.0",
        "S8.0","E8.0","I8.0","A8.0","R8.0",
        "S9.0","E9.0","I9.0","A9.0","R9.0",
        "S10.0","E10.0","I10.0","A10.0","R10.0",
        "incid.0","N0"),
    zeronames=c("incid","incidU","incidV","asymV","newV","foival","Str0","Sout","Sin"),
    initializer=function(params, t0, ...){
      all.state.names <- c("S", sprintf("S%d", 1:10),
                          "E", sprintf("E%d", 1:10),
                          "I", sprintf("I%d", 1:10),
                          "A", sprintf("A%d", 1:10),
                          "R", sprintf("R%d", 1:10),
                          "incid","incidU","incidV","asymV","newV","foival","Str0","Sout","Sin")
      x0 <- setNames(numeric(length(all.state.names)), all.state.names)
      fracs <- params[c("S.0","S1.0","S2.0","S3.0","S4.0","S5.0","S6.0","S7.0","S8.0","S9.0","S10.0",
                        "E.0","E1.0","E2.0","E3.0","E4.0","E5.0","E6.0","E7.0","E8.0","E9.0","E10.0",
                        "I.0","I1.0","I2.0","I3.0","I4.0","I5.0","I6.0","I7.0","I8.0","I9.0","I10.0",
                        "A.0","A1.0","A2.0","A3.0","A4.0","A5.0","A6.0","A7.0","A8.0","A9.0","A10.0",
                        "R.0","R1.0","R2.0","R3.0","R4.0","R5.0","R6.0","R7.0","R8.0","R9.0","R10.0")]
      initpop <- params[c("N0")]
      inc <- params[c("incid.0")]
      x0[c("S","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10",
           "E","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10",
           "I","I1","I2","I3","I4","I5","I6","I7","I8","I9","I10",
           "A","A1","A2","A3","A4","A5","A6","A7","A8","A9","A10",
           "R","R1","R2","R3","R4","R5","R6","R7","R8","R9","R10")] <- round(initpop*fracs/sum(fracs)) ## rm incid from statenames
      x0[c("incid")] <- round(inc)
      # print("******")
      # print(x0)
      # print(params)
      x0
    })

  return(my.mod)
}