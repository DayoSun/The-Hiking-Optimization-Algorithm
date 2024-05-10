function L_xmin=Interpolation(xi,xj,xk,fi,fj,fk,l,u)
%%%%%%%%%%%% Quadratic interpolation %%%%%%%%%%%
  %Eq.(5)
  a=(xj^2-xk^2)*fi+(xk^2-xi^2)*fj+(xi^2-xj^2)*fk;
  b=2*((xj-xk)*fi+(xk-xi)*fj+(xi-xj)*fk);

  L_xmin=a/(b+eps); 
  if isnan(L_xmin) || isinf(L_xmin) || L_xmin>u || L_xmin<l
      L_xmin=(rand*(u-l)+l);
  end
   
 
