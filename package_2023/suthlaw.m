function [rho, nu, mu] = suthlaw(patm,T)

%%%%%%%% 	[rho, nu, mu] = suthlaw(patm(Pa),T(C))
%%%%%%%%
%%%%%%%% Calculates air properties using Sutherland's relation
%%%%%%%%
%%%%%%%%                        patm (Pa)
%%%%%%%%    rho = -------------------------------------
%%%%%%%%            R (universal gas constant) x T (K)
%%%%%%%%
%%%%%%%% Valid for temperatures 0<T<500C and pressures < 34.5 bar.
%%%%%%%% 
%%%%%%%%             1.4578*10^-6 x T(K) 
%%%%%%%%     mu = ---------------------------
%%%%%%%%                 T(K) + 110.4
%%%%%%%%
%%%%%%%%     nu = mu/rho 


TK = T+273.15;
rho = patm./(287.*TK);
%mu = 1.4578e-6.*TK.^1.5./(TK+110.4);
mu = 1.7894e-5.*(TK/273.11).^1.5.*(273.11+110.56)./(TK+110.56);

nu = mu./rho;