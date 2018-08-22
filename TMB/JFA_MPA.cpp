#include <TMB.hpp>
// functions

template <class Type> Type square(Type x){return x*x;}

template <class Type> Type selec(Type l95, Type l50, Type size){
  return 1 / (1 + exp(( - log(19) * (size - l50)) / (l95 - l50)));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA SECTION
  DATA_INTEGER(nyears);
  DATA_INTEGER(nclass);
  DATA_VECTOR(lvec);
  DATA_VECTOR(weight);
  DATA_MATRIX(tmatrix);
  DATA_SCALAR(M);
  DATA_VECTOR(catch_obs);
  DATA_VECTOR(pcll); // Years with catch at size distribution
  DATA_MATRIX(cal_obs);
  DATA_VECTOR(pcpue); // Years with observed cpue
  DATA_MATRIX(cpue_obs);
  DATA_VECTOR(selec_f);
    // for projections
  DATA_INTEGER(nyproj);
  DATA_SCALAR(Fproj);

  // PARAMETER SECTION
  PARAMETER(dummy);
  PARAMETER(LogRbar);
  PARAMETER_VECTOR(LogNinit);
  PARAMETER_VECTOR(LOG_F);
  PARAMETER_VECTOR(Eps);
  PARAMETER_VECTOR(selec_t);
  // CREATING ARRAYS
  matrix<Type> N(nyears + nyproj + 1, nclass);
  matrix<Type> F(nyears + nyproj, nclass);
  matrix<Type> Z(nyears + nyproj, nclass);
  matrix<Type> CAL(nyears + nyproj, nclass);
  vector<Type> sel_f(nclass);
  vector<Type> sel_t(nclass);
  vector<Type> catch_w(nyears + nyproj);
  vector<Type> BioPred(nyears + nyproj);

  Type CALtot;
  Type Penal;
  Type LikeCatch;
  Type LikeBio;
  Type LikeCAL;
  Type obj_fun;

 // selectivities
  for( int y  = 0; y < nyears; y ++){
    for( int len = 0; len < nclass; len++){
      F(y, len)  =  exp(LOG_F(y)) * selec(selec_f(0), selec_f(1), lvec(len));
      sel_t(len)  =  selec(selec_t(0), selec_t(1), lvec(len));
    }
  }
    obj_fun = square(1);
  // // End of specifications section
  // // =============================

  // // First set F and Z by size-classs (note that Fproj applies after year nyears)
  // for (int Iyear=0; Iyear<nyears+nyproj; Iyear++)
  //  for (int Iclass=0;Iclass<nclass;Iclass++)
  //   {
  //    if (Iyear < nyears)
  // 	  F(Iyear,Iclass) = exp(LogFullF(Iyear))*S(Iclass);
  // 	 else
  // 	  F(Iyear,Iclass) = Fproj*S(Iclass);
  // 	 Z(Iyear,Iclass) = M + F(Iyear,Iclass);

  //   }

  // // Now set the N matrix
  // for (int Iclass=0;Iclass<nclass;Iclass++) N(0,Iclass) = exp(LogNinit(Iclass));
  // for (int Iyear=0;Iyear<nyears+nyproj;Iyear++)
  //  {
  //   // Catch-at-length
  //   CALtot = 0; catch_w(Iyear) = 0;
  //   for (int Iclass=0;Iclass<nclass;Iclass++)
  //    {
  //     CAL(Iyear,Iclass) = F(Iyear,Iclass)/Z(Iyear,Iclass)*N(Iyear,Iclass)*(1.0-exp(-Z(Iyear,Iclass)));
  //     CALtot += CAL(Iyear,Iclass);
  //     catch_w(Iyear) += Weight(Iclass)*CAL(Iyear,Iclass);
  //    }
  //   for (int Iclass=0;Iclass<nclass;Iclass++) CAL(Iyear,Iclass) /= CALtot;

  //   // Numbers-at-age
  //   for (int Iclass=0;Iclass<nclass;Iclass++)
  //    {
  //     N(Iyear+1,Iclass) = 0;
  //     for (int Jclass=0;Jclass<nclass;Jclass++)
  //      N(Iyear+1,Iclass) += N(Iyear,Jclass)*exp(-Z(Iyear,Jclass))*X(Jclass,Iclass);
  //    }

  //   // Recruitment (watch for the index for Eps - and N)
  //   N(Iyear+1,0) += exp(LogRbar)*exp(Eps[Iyear]);
  //  }

  // // Catch Likelihood
  // Type SS = 0;
  // for (int Iyear=0; Iyear<nyears; Iyear++)
  //  SS += square(log(catch_wObs(Iyear)) - log(catch_w(Iyear)));
  // LikeCatch = SS /(2.0*0.05*0.05);

  // // Biomass predictions
  // for (int Iyear=0; Iyear<(nyears+nyproj); Iyear++)
  //  {
  //   BioPred(Iyear) = 0;
  //   for (int Iclass=0;Iclass<nclass;Iclass++) BioPred(Iyear) += N(Iyear,Iclass)*SurveyS(Iclass)*Weight(Iclass);
  //  }
  // Type Top = 0; Type Bot = 0; Type q;
  // for (int Iyear=0; Iyear<nyears; Iyear++)
  //  { Top += log(BioIndex(Iyear)/BioPred(Iyear)); Bot += 1.0; }
  // q = exp(Top/Bot);

  // // Likelihood
  // SS = 0;
  // for (int Iyear=0; Iyear<nyears; Iyear++)
  //  SS += square(log(BioIndex(Iyear))-log(q*BioPred(Iyear)));
  // LikeBio = SS/(2*BioSig*BioSig);

  // // CAL Likelihood
  // LikeCAL = 0;
  // for (int Iyear=0; Iyear<nyears; Iyear++)
  //  for (int Iclass=0;Iclass<nclass;Iclass++)
  //   if (CALObs(Iyear,Iclass) > 0)
  //    LikeCAL -= Neff*CALObs(Iyear,Iclass)*log(CAL(Iyear,Iclass)/CALObs(Iyear,Iclass));

  // // Recruitment penalty (include years after nyears)
  // Penal = 0;
  // for (int Iyear=0; Iyear<(nyears+nyproj); Iyear++)
  //  Penal += Eps(Iyear)*Eps(Iyear);
  // Penal = Penal / (2.0*0.6*0.6);

  // obj_fun = dummy*dummy + LikeCatch+LikeBio+LikeCAL+Penal;

  // // Stuff to report
  // REPORT(N);
  // REPORT(LikeCatch);
  // REPORT(LikeBio);
  // REPORT(LikeCAL);
  // REPORT(Penal);
  // REPORT(BioPred);
  // REPORT(obj_fun);

  return(obj_fun);
}
