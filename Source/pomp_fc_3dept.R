## pomp model in C for the 3-department model

## ---------------------------------------- ##
## Step function used in simulating process 

rSim.step <- Csnippet('
  // transition rates
  double Srate[5]; 
  double Erate[6];
  double Irate[5];
  double Arate[5];
  double Rrate[5];
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

  // transition numbers
  double Strans[5]; 
  double Etrans[6];
  double Itrans[5];
  double Atrans[5];
  double Rtrans[5];
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

  // vac numbers
  double vac1 = 0.0;
  double vac2 = 0.0;
  double vac3 = 0.0;

  // some population demonitors
  int Nnv = S + E + I + A + R;
  int N1 = S1 + E1 + I1 + A1 + R1; // Centre
  int N2 = S2 + E2 + I2 + A2 + R2; // Artibonite
  int N3 = S3 + E3 + I3 + A3 + R3; // Ouest
  int N = Nnv + N1 + N2 + N3;
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

  // make seasonal beta term for current time
  double mybeta = beta1*seas1 + beta2*seas2 + beta3*seas3 + beta4*seas4 + beta5*seas5 + beta6*seas6;
  double foi = pow(I+I1+I2+I3+(1-kappa)*(A+A1+A2+A3), nu)*mybeta/N;  

  // make thetav (artificially from covariate dataframe)
  double thetav_d1 = ve_d1;
  double thetav_d2 = ve_d2;
  double thetav_d3 = ve_d3;

  //compute the rates for all the transitions
  Srate[0]= foi;  //S -> E
  Srate[1]= delta;
  Srate[2]= vac1;
  Srate[3]= vac2;
  Srate[4]= vac3;

  Erate[0]= sigma*(1-theta0); // E -> I
  Erate[1]= sigma*theta0; // E -> A
  Erate[2]= delta;
  Erate[3]= vac1;
  Erate[4]= vac2;
  Erate[5]= vac3;

  Irate[0]= gamma; // I -> R
  Irate[1]= delta;
  Irate[2]= vac1;
  Irate[3]= vac2;
  Irate[4]= vac3;

  Arate[0]= gamma; // A -> R
  Arate[1]= delta;
  Arate[2]= vac1;
  Arate[3]= vac2;
  Arate[4]= vac3;

  Rrate[0]= alpha; // R -> S exponential waning immunity from natural infection
  Rrate[1]= delta;
  Rrate[2]= vac1;
  Rrate[3]= vac2;
  Rrate[4]= vac3;

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

  // compute the transition numbers
  reulermultinom(5,S,&Srate[0],dt,&Strans[0]);
  reulermultinom(6,E,&Erate[0],dt,&Etrans[0]);
  reulermultinom(5,I,&Irate[0],dt,&Itrans[0]);
  reulermultinom(5,A,&Arate[0],dt,&Atrans[0]);
  reulermultinom(5,R,&Rrate[0],dt,&Rtrans[0]);

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

  
  // balance the equations
  S += -Strans[0] - Strans[1] - Strans[2] - Strans[3] - Strans[4] + Rtrans[0] + births;
  E += -Etrans[0] - Etrans[1] - Etrans[2] - Etrans[3] - Etrans[4] - Etrans[5] + Strans[0];
  I += -Itrans[0] - Itrans[1] - Itrans[2] - Itrans[3] - Itrans[4] + Etrans[0];
  A += -Atrans[0] - Atrans[1] - Atrans[2] - Atrans[3] - Atrans[4] + Etrans[1];
  R += -Rtrans[0] - Rtrans[1] - Rtrans[2] - Rtrans[3] - Rtrans[4] + Itrans[0] + Atrans[0];

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

  incid += Etrans[0] + E1trans[0] + E2trans[0] + E3trans[0]; // incidence is cumulative entries into I state
  incidU += Etrans[0];
  incidV += E1trans[0] + E2trans[0] + E3trans[0];
  asymV += E1trans[1] + E2trans[1] + E3trans[1];
  newV += Strans[2] + Strans[3] + Strans[4] + Etrans[3] + Etrans[4] + Etrans[5] + Itrans[2] + Itrans[3] + Itrans[4] + Atrans[2] + Atrans[3] + Atrans[4] + Rtrans[2] + Rtrans[3] + Rtrans[4];
  foival += foi;
  Str0 += Strans[0];
  Sout += Strans[0] + Strans[1] + Strans[2] + Strans[3];
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
    statenames=c("S", sprintf("S%d", 1:3),
                 "E", sprintf("E%d", 1:3),
                 "I", sprintf("I%d", 1:3),
                 "A", sprintf("A%d", 1:3),
                 "R", sprintf("R%d", 1:3),
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
        "incid.0","N0"),
    zeronames=c("incid","incidU","incidV","asymV","newV","foival","Str0","Sout","Sin"),
    initializer=function(params, t0, ...){
      all.state.names <- c("S", sprintf("S%d", 1:3),
                          "E", sprintf("E%d", 1:3),
                          "I", sprintf("I%d", 1:3),
                          "A", sprintf("A%d", 1:3),
                          "R", sprintf("R%d", 1:3),
                          "incid","incidU","incidV","asymV","newV","foival","Str0","Sout","Sin")
      x0 <- setNames(numeric(length(all.state.names)), all.state.names)
      fracs <- params[c("S.0","S1.0","S2.0","S3.0",
                        "E.0","E1.0","E2.0","E3.0",
                        "I.0","I1.0","I2.0","I3.0",
                        "A.0","A1.0","A2.0","A3.0",
                        "R.0","R1.0","R2.0","R3.0")]
      initpop <- params[c("N0")]
      inc <- params[c("incid.0")]
      x0[c("S","S1","S2","S3",
           "E","E1","E2","E3",
           "I","I1","I2","I3",
           "A","A1","A2","A3",
           "R","R1","R2","R3")] <- round(initpop*fracs/sum(fracs)) ## rm incid from statenames
      x0[c("incid")] <- round(inc)
      # print("******")
      # print(x0)
      # print(params)
      x0
    })

  return(my.mod)
}