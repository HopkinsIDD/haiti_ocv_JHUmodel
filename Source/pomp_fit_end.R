## pomp model in C for the endemic fit

## ---------------------------------------- ##
## Step function used in simulating process 

rSim.step <- Csnippet('
  // transition rates
  double Srate[2]; 
  double Erate[3];
  double Irate[2];
  double Arate[2];
  double Rrate[2];

  // transition numbers
  double Strans[2]; 
  double Etrans[3];
  double Itrans[2];
  double Atrans[2];
  double Rtrans[2];

  // some population demonitors
  int N = S + E + I + A + R;
  int births = rpois(mu*N*dt);

  // make seasonal beta term for current time
  double mybeta = beta1*seas1 + beta2*seas2 + beta3*seas3 + beta4*seas4 + beta5*seas5 + beta6*seas6;
  double foi = pow(I, nu)*mybeta/N;  

  //compute the rates for all the transitions
  Srate[0]= foi;  //S -> E
  Srate[1]= delta;

  Erate[0]= sigma*(1-theta0); // E -> I
  Erate[1]= sigma*theta0; // E -> A
  Erate[2]= delta;

  Irate[0]= gamma; // I -> R
  Irate[1]= delta;

  Arate[0]= gamma; // A -> R
  Arate[1]= delta;

  Rrate[0]= alpha; // R -> S exponential waning immunity from natural infection
  Rrate[1]= delta;

  // compute the transition numbers
  reulermultinom(2,S,&Srate[0],dt,&Strans[0]);
  reulermultinom(3,E,&Erate[0],dt,&Etrans[0]);
  reulermultinom(2,I,&Irate[0],dt,&Itrans[0]);
  reulermultinom(2,A,&Arate[0],dt,&Atrans[0]);
  reulermultinom(2,R,&Rrate[0],dt,&Rtrans[0]);

  // balance the equations
  S += -Strans[0] - Strans[1] + Rtrans[0] + births;
  E += -Etrans[0] - Etrans[1] - Etrans[2] + Strans[0];
  I += -Itrans[0] - Itrans[1] + Etrans[0];
  A += -Atrans[0] - Atrans[1] + Etrans[1];
  R += -Rtrans[0] - Rtrans[1] + Itrans[0] + Atrans[0];

  incid += Etrans[0]; // incidence is cumulative entries into I state
  //if (R_FINITE(incid)) 
  //  Rprintf(\" Srate %i  \\n\", -Strans[0] - Strans[1] + Rtrans[0] + births);
')

## ---------------------------------------- ##
## Deterministic skeleton of process

sim.skel <- Csnippet('
  // transition rates
  double Srate[2]; 
  double Erate[3];
  double Irate[2];
  double Arate[2];
  double Rrate[2];

  // transition terms
  double Sterm[2]; 
  double Eterm[3];
  double Iterm[2];
  double Aterm[2];
  double Rterm[2];

  // some population demonitors
  int N = S + E + I + A + R;
  double births = mu*N;

  // make seasonal beta term for current time
  double mybeta = beta1*seas1 + beta2*seas2 + beta3*seas3 + beta4*seas4 + beta5*seas5 + beta6*seas6;
  double foi = pow(I, nu)*mybeta/N;  

  //compute the rates for all the transitions
  Srate[0]= foi;  //S -> E
  Srate[1]= delta;

  Erate[0]= sigma*(1-theta0); // E -> I
  Erate[1]= sigma*theta0; // E -> A
  Erate[2]= delta;

  Irate[0]= gamma; // I -> R
  Irate[1]= delta;

  Arate[0]= gamma; // A -> R
  Arate[1]= delta;

  Rrate[0]= alpha; // R -> S exponential waning immunity from natural infection
  Rrate[1]= delta;

  // compute the transition terms
  for(int i=0; i < 2; i++){
    double term = Srate[i]*S;
    Sterm[i] = term;}
  for(int i=0; i < 3; i++){
    double term = Erate[i]*E;
    Eterm[i] = term;}
  for(int i=0; i < 2; i++){
    double term = Irate[i]*I;
    Iterm[i] = term;}
  for(int i=0; i < 2; i++){
    double term = Arate[i]*A;
    Aterm[i] = term;}
  for(int i=0; i < 2; i++){
    double term = Rrate[i]*R;
    Rterm[i] = term;}


  // balance the equations
  DS = -Sterm[0] - Sterm[1] + Rterm[0] + births; 
  DE = -Eterm[0] - Eterm[1] - Eterm[2] + Sterm[0];
  DI = -Iterm[0] - Iterm[1] + Eterm[0];
  DA = -Aterm[0] - Aterm[1] + Eterm[1];
  DR = -Rterm[0] - Rterm[1] + Iterm[0] + Aterm[0];
  Dincid = Eterm[0]; // incidence is cumulative entries into I state

')

## ---------------------------------------- ##
## Simulation of the measurement process

rMeas <- Csnippet('
  cases = rnbinom_mu(theta, rho*incid+1e-20);
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } else {
    cases = 0.0;
  }
')

## ---------------------------------------- ##
## Measurement process density evaluation

dMeas <- Csnippet('
  lik = dnbinom_mu(cases, theta, rho*incid+1e-20, give_log); // 4th arg is give_log
  //lik = (give_log) ? lik : exp(lik);
  if (!R_FINITE(lik)) 
      Rprintf(\"%i %i %f %f %f %f %i\\n\", lik, cases, theta, rho, incid, rho*incid, give_log);
')

## ---------------------------------------- ##
## Parameter Transformations 

## from estimation scale to natural scale
par.trans.fromEstim <- Csnippet('
  Tbeta1=exp(beta1);
  Tbeta2=exp(beta2);
  Tbeta3=exp(beta3);
  Tbeta4=exp(beta4);
  Tbeta5=exp(beta5);
  Tbeta6=exp(beta6);
  Trho=expit(rho);
  Ttheta=exp(theta);
  Tnu=expit(nu);
  Tsigma=exp(sigma);
  Tgamma=exp(gamma);
  Tmu=exp(mu);
  Tdelta=exp(delta);
  Ttheta0=expit(theta0);
  Talpha=exp(alpha); 
  //from_log_barycentric(&TS_0, &S_0, 5);
')

## from natural scale to estimation scale
par.trans.toEstim <- Csnippet('
  Tbeta1=log(beta1);
  Tbeta2=log(beta2);
  Tbeta3=log(beta3);
  Tbeta4=log(beta4);
  Tbeta5=log(beta5);
  Tbeta6=log(beta6);
  Trho=logit(rho);
  Ttheta=log(theta);
  Tnu=logit(nu);
  Tsigma=log(sigma);
  Tgamma=log(gamma);
  Tmu=log(mu);
  Tdelta=log(delta);
  Ttheta0=logit(theta0);
  Talpha=log(alpha); 
  //to_log_barycentric(&TS_0, &S_0, 5);
')

## ---------------------------------------- ##
## Build pomp model

build.end.mod <- function(pop, ##  = get.haiti.pop()
                          dat, ##  = get.mspp.agg.data()
                          my.times = "week_end", 
                          covar.times = "time_end",
                          my.t0 = 0,
                          covar){

  my.mod <- pomp(
    data = dat,
    times = my.times,
    t0 = my.t0,
    dmeasure = dMeas,
    rmeasure = rMeas,
    rprocess = euler.sim(step.fun=rSim.step, delta.t=1/7),
    skeleton = vectorfield(sim.skel),
    covar = covar,
    tcovar = covar.times, 
    fromEstimationScale=par.trans.fromEstim,
    toEstimationScale=par.trans.toEstim,
    statenames=c(
        "S","E","I","A","R","incid"),
    paramnames=c(
        "rho","theta",
        "beta1","beta2","beta3",
        "beta4",
        "beta5",
        "beta6",
        "gamma","sigma","theta0","alpha","mu","delta","nu",
        "S.0","E.0","I.0","A.0","R.0","incid.0","N0"),
    zeronames=c("incid"),
    initializer=function(params, t0, ...){
      all.state.names <- c("S","E","I","A","R","incid")
      x0 <- setNames(numeric(length(all.state.names)), all.state.names)
      fracs <- params[c("S.0","E.0","I.0","A.0","R.0")]
      N <- params[c("N0")]
      inc <- params[c("incid.0")]
      x0[-length(x0)] <- round(N*fracs/sum(fracs)) ## rm incid from statenames
      x0[length(x0)] <- round(inc)
      print(x0)
      x0
    })

  return(my.mod)
}