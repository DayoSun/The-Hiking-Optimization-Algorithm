function [lb,ub,dim,fobj]=hybrid(F,nVar)
global initial_flag
initial_flag=0;
switch F
    case 'F34'
        fobj = @hybrid_func1;
        lb=-5;
        ub=5;
        dim=nVar;
        
    case 'F35'
        fobj = @hybrid_rot_func1;
        lb=-5;
        ub=5;
        dim=nVar;
        
    case 'F36'
        fobj = @hybrid_rot_func1_noise;
        lb=-5;
        ub=5;
        dim=nVar;
        
    case 'F37'
        fobj = @hybrid_rot_func2;
        lb=-5;
        ub=5;
        dim=nVar;
        
    case 'F38'
        fobj = @hybrid_rot_func2_narrow;
        lb=-5;
        ub=5;
        dim=nVar;
        
    case 'F39'
        fobj = @hybrid_rot_func2_onbound;
        lb=-5;
        ub=5;
        dim=nVar;
        
    case 'F40'
        fobj = @hybrid_rot_func3;
        lb=-5;
        ub=5;
        dim=nVar;
        
    case 'F41'
        fobj = @hybrid_rot_func3_highcond;
        lb=-5;
        ub=5;
        dim=nVar;
        
    case 'F42'
        fobj = @hybrid_rot_func3_noncont;
        lb=-5;
        ub=5;
        dim=nVar;
        
    case 'F43'
        fobj = @hybrid_rot_func4;
        lb=-5;
        ub=5;
        dim=nVar;
        
    case 'F44'
        fobj = @hybrid_rot_func3;
        lb=-5;
        ub=5;
        dim=nVar;
end
%---------------------------------------------------
%   15.Hybrid Composition Function 1
function fit=hybrid_func1(x)
global initial_flag
persistent  fun_num func o sigma lamda bias M
if initial_flag==0
    [ps,D]=size(x);
    initial_flag=1;
    fun_num=10;

    o=-5+10*rand(fun_num,D);
    
    func.f1=@(x) frastrigin(x);
    func.f2=@(x) frastrigin(x);
    func.f3=@(x) fweierstrass(x);
    func.f4=@(x) fweierstrass(x);
    func.f5=@(x) fgriewank(x);
    func.f6=@(x) fgriewank(x);
    func.f7=@(x) fackley(x);
    func.f8=@(x) fackley(x);
    func.f9=@(x) fsphere(x);
    func.f10=@(x) fsphere(x);

    bias=((1:fun_num)-1).*100;
    sigma=ones(1,fun_num);
    lamda=[1; 1; 10; 10; 5/60; 5/60; 5/32; 5/32; 5/100; 5/100];
    lamda=repmat(lamda,1,D);
    for i=1:fun_num
        eval(['M.M' int2str(i) '=diag(ones(1,D));']);
    end
end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);

%---------------------------------------------------------------------
%   16.Rotated Hybrid Composition Function 1	
function fit=hybrid_rot_func1(x)
global initial_flag
persistent  fun_num func o sigma lamda bias M
if initial_flag==0
    [ps,D]=size(x);
    initial_flag=1;
    fun_num=10;
    %load hybrid_func1_data % saved the predefined optima
%     if length(o(1,:))>=D
%          o=o(:,1:D);
%     else
         o=-5+10*rand(fun_num,D);
%     end
    func.f1=@(x) frastrigin(x);
    func.f2=@(x) frastrigin(x);
    func.f3=@(x) fweierstrass(x);
    func.f4=@(x) fweierstrass(x);
    func.f5=@(x) fgriewank(x);
    func.f6=@(x) fgriewank(x);
    func.f7=@(x) fackley(x);
    func.f8=@(x) fackley(x);
    func.f9=@(x) fsphere(x);
    func.f10=@(x) fsphere(x);
    bias=((1:fun_num)-1).*100;
    sigma=ones(1,fun_num); 
    lamda=[1; 1; 10; 10; 5/60; 5/60; 5/32; 5/32; 5/100; 5/100];
    lamda=repmat(lamda,1,D);
    c=[2,2,2,2,2,2,2,2,2,2,2];
    if D==2 %,load hybrid_func1_M_D2,
    elseif D==10 %,load hybrid_func1_M_D10,
    elseif D==30 %,load hybrid_func1_M_D30,
    elseif D==50 %,load hybrid_func1_M_D50,
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
%----------------------------------------------------------------
%   17.	Rotated Hybrid Composition Function 1 with Noise in Fitness	
function fit=hybrid_rot_func1_noise(x)
[ps,D]=size(x);
fit=hybrid_rot_func1(x).*(1+0.2.*abs(normrnd(0,1,ps,1)));
%----------------------------------------------------------------
%   18.	Rotated Hybrid Composition Function 2
function fit=hybrid_rot_func2(x)
global initial_flag
persistent  fun_num func o sigma lamda bias M
if initial_flag==0
    [ps,D]=size(x);
    initial_flag=1;
    fun_num=10;
%     load hybrid_func2_data % saved the predefined optima
%     if length(o(1,:))>=D
%          o=o(:,1:D);
%     else
         o=-5+10*rand(fun_num,D);
%     end
    o(10,:)=0;
    func.f1=@(x) fackley(x);
    func.f2=@(x) fackley(x);
    func.f3=@(x) frastrigin(x);
    func.f4=@(x) frastrigin(x);
    func.f5=@(x) fsphere(x);
    func.f6=@(x) fsphere(x);
    func.f7=@(x) fweierstrass(x);
    func.f8=@(x) fweierstrass(x);
    func.f9=@(x) fgriewank(x);
    func.f10=@(x) fgriewank(x);
    bias=((1:fun_num)-1).*100;
    sigma=[1 2 1.5 1.5 1 1 1.5 1.5 2 2];
    lamda=[2*5/32; 5/32; 2*1; 1; 2*5/100; 5/100; 2*10; 10; 2*5/60; 5/60];
    lamda=repmat(lamda,1,D);
    c=[2 3 2 3 2 3 20 30 200 300];
    if D==2%,load hybrid_func2_M_D2,
    elseif D==10%,load hybrid_func2_M_D10,
    elseif D==30%,load hybrid_func2_M_D30,
    elseif D==50%,load hybrid_func2_M_D50,
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
%----------------------------------------------------------------
%   19.	Rotated Hybrid Composition Function 2 with a Narrow Basin for the Global Optimum
function fit=hybrid_rot_func2_narrow(x)
global initial_flag
persistent  fun_num func o sigma lamda bias M
if initial_flag==0
    [ps,D]=size(x);
    initial_flag=1;
    fun_num=10;
%     load hybrid_func2_data % saved the predefined optima
%     if length(o(1,:))>=D
%          o=o(:,1:D);
%     else
         o=-5+10*rand(fun_num,D);
%     end
    o(10,:)=0;
    func.f1=@(x) fackley(x);
    func.f2=@(x) fackley(x);
    func.f3=@(x) frastrigin(x);
    func.f4=@(x) frastrigin(x);
    func.f5=@(x) fsphere(x);
    func.f6=@(x) fsphere(x);
    func.f7=@(x) fweierstrass(x);
    func.f8=@(x) fweierstrass(x);
    func.f9=@(x) fgriewank(x);
    func.f10=@(x) fgriewank(x);
    bias=((1:fun_num)-1).*100;
    sigma=[0.1 2 1.5 1.5 1 1 1.5 1.5 2 2];
    lamda=[0.1*5/32; 5/32; 2*1; 1; 2*5/100; 5/100; 2*10; 10; 2*5/60; 5/60];
    lamda=repmat(lamda,1,D);
    c=[2 3 2 3 2 3 20 30 200 300];
    if D==2%,load hybrid_func2_M_D2,
    elseif D==10%,load hybrid_func2_M_D10,
    elseif D==30%,load hybrid_func2_M_D30,
    elseif D==50%,load hybrid_func2_M_D50,
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
%----------------------------------------------------------------
%   20.	Rotated Hybrid Composition Function 2 with the Global Optimum on the Bounds	
function fit=hybrid_rot_func2_onbound(x)
global initial_flag
persistent  fun_num func o sigma lamda bias M
if initial_flag==0
    [ps,D]=size(x);
    initial_flag=1;
    fun_num=10;
%     load hybrid_func2_data % saved the predefined optima,
%     if length(o(1,:))>=D
%          o=o(:,1:D);
%     else
         o=-5+10*rand(fun_num,D);
%     end
    o(10,:)=0;
    o(1,2.*[1:floor(D/2)])=5;
    func.f1=@(x) fackley(x);
    func.f2=@(x) fackley(x);
    func.f3=@(x) frastrigin(x);
    func.f4=@(x) frastrigin(x);
    func.f5=@(x) fsphere(x);
    func.f6=@(x) fsphere(x);
    func.f7=@(x) fweierstrass(x);
    func.f8=@(x) fweierstrass(x);
    func.f9=@(x) fgriewank(x);
    func.f10=@(x) fgriewank(x);
    bias=((1:fun_num)-1).*100;
    sigma=[1 2 1.5 1.5 1 1 1.5 1.5 2 2];
    lamda=[2*5/32; 5/32; 2*1; 1; 2*5/100; 5/100; 2*10; 10; 2*5/60; 5/60];
    lamda=repmat(lamda,1,D);
    c=[2 3 2 3 2 3 20 30 200 300];
    if D==2%,load hybrid_func2_M_D2,
    elseif D==10%,load hybrid_func2_M_D10,
    elseif D==30%,load hybrid_func2_M_D30,
    elseif D==50%,load hybrid_func2_M_D50,
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
%-------------------------------------------------
%    21.Rotated Hybrid Composition Function 3		
function fit=hybrid_rot_func3(x)
global initial_flag
persistent  fun_num func o sigma lamda bias M
if initial_flag==0
    [ps,D]=size(x);
    initial_flag=1;
    fun_num=10;
%     load hybrid_func3_data % saved the predefined optima, a 10*1000 matrix;
%     if length(o(1,:))>=D
%          o=o(:,1:D);
%     else
         o=-5+10*rand(fun_num,D);
%     end
    func.f1=@(x) fE_ScafferF6(x);
    func.f2=@(x) fE_ScafferF6(x);
    func.f3=@(x) frastrigin(x);
    func.f4=@(x) frastrigin(x);
    func.f5=@(x) fEF8F2(x);
    func.f6=@(x) fEF8F2(x);
    func.f7=@(x) fweierstrass(x);
    func.f8=@(x) fweierstrass(x);
    func.f9=@(x) fgriewank(x);
    func.f10=@(x) fgriewank(x);
    bias=((1:fun_num)-1).*100;
    sigma=[1,1,1,1,1,2,2,2,2,2];
    lamda=[5*5/100; 5/100; 5*1; 1; 5*1; 1; 5*10; 10; 5*5/200; 5/200];
    lamda=repmat(lamda,1,D);
    c=ones(1,D);
    if D==2%,load hybrid_func3_M_D2,
    elseif D==10%,load hybrid_func3_M_D10,
    elseif D==30%,load hybrid_func3_M_D30,
    elseif D==50%,load hybrid_func3_M_D50,
    else 
        for i=1:fun_num
            A=normrnd(0,1,D,D);
            eval(['M.M' int2str(i) '=cGram_Schmidt(A));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
%-----------------------------------------
%   22.	Rotated Hybrid Composition Function 3 with High Condition Number Matrix
function fit=hybrid_rot_func3_highcond(x)
global initial_flag
persistent  fun_num func o sigma lamda bias M
if initial_flag==0
    [ps,D]=size(x);
    initial_flag=1;
    fun_num=10;
%     load hybrid_func3_data % saved the predefined optima, a 10*1000 matrix;
%     if length(o(1,:))>=D
%          o=o(:,1:D);
%     else
         o=-5+10*rand(fun_num,D);
%     end
    func.f1=@(x) fE_ScafferF6(x);
    func.f2=@(x) fE_ScafferF6(x);
    func.f3=@(x) frastrigin(x);
    func.f4=@(x) frastrigin(x);
    func.f5=@(x) fEF8F2(x);
    func.f6=@(x) fEF8F2(x);
    func.f7=@(x) fweierstrass(x);
    func.f8=@(x) fweierstrass(x);
    func.f9=@(x) fgriewank(x);
    func.f10=@(x) fgriewank(x);
    bias=((1:fun_num)-1).*100;
    sigma=[1,1,1,1,1,2,2,2,2,2];
    lamda=[5*5/100; 5/100; 5*1; 1; 5*1; 1; 5*10; 10; 5*5/200; 5/200];
    lamda=repmat(lamda,1,D);
    c=[10 20 50 100 200 1000 2000 3000 4000 5000];
    if D==2%,load hybrid_func3_HM_D2,
    elseif D==10%,load hybrid_func3_HM_D10,
    elseif D==30%,load hybrid_func3_HM_D30,
    elseif D==50%,load hybrid_func3_HM_D50,
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
%-----------------------------------------
%   23.	Non-Continuous Rotated Hybrid Composition Function 3
function fit=hybrid_rot_func3_noncont(x)
global initial_flag
persistent  o 
[ps,D]=size(x);
if initial_flag==0
%     load hybrid_func3_data % saved the predefined optima, a 10*1000 matrix;
%     o=o(1,1:D);
end
o=repmat(o,ps,1);
x=(abs(x-o)<0.5).*x+(abs(x-o)>=0.5).*(round(x.*2)./2);
fit=hybrid_rot_func3(x);
%-----------------------------------------
%   24.	Rotated Hybrid Composition Function 4	
function fit=hybrid_rot_func4(x)
global initial_flag
persistent  fun_num func o sigma lamda bias M
if initial_flag==0
    [ps,D]=size(x);
    initial_flag=1;
    fun_num=10;
%     load hybrid_func4_data % saved the predefined optima, a 10*1000 matrix;
%     if length(o(1,:))>=D
%          o=o(:,1:D);
%     else
         o=-5+10*rand(fun_num,D);
%     end
    func.f1=@(x) fweierstrass(x);
    func.f2=@(x) fE_ScafferF6(x);
    func.f3=@(x) fEF8F2(x);
    func.f4=@(x) fackley(x);
    func.f5=@(x) frastrigin(x);
    func.f6=@(x) fgriewank(x);
    func.f7=@(x) fE_ScafferF6_noncont(x);
    func.f8=@(x) frastrigin_noncont(x);
    func.f9=@(x) felliptic(x);
    func.f10=@(x) fsphere_noise(x);
    bias=((1:fun_num)-1).*100;
    sigma=[2,2,2,2,2,2,2,2,2,2];
    lamda=[10; 5/20; 1; 5/32; 1; 5/100 ; 5/50; 1; 5/100; 5/100; ];
    lamda=repmat(lamda,1,D);
    c=[100 50 30 10 5 5 4 3 2 2];
    if D==2%,load hybrid_func4_M_D2,
    elseif D==10%,load hybrid_func4_M_D10,
    elseif D==30%,load hybrid_func4_M_D30,
    elseif D==50%,load hybrid_func4_M_D50,
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M);
%----------------------------------
function fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M)
[ps,D]=size(x);
weight=[];
for i=1:fun_num
    oo=repmat(o(i,:),ps,1);
    weight(:,i)=exp(-sum((x-oo).^2,2)./2./(D*sigma(i)^2)); %#ok<*AGROW>
    %[-sum((x-oo).^2,2), 2*D*sigma(i)^2, weight(:,i) ]
end

[tmp,tmpid]=sort(weight,2);
for i=1:ps
    weight(i,:)=(weight(i,:)==tmp(i,fun_num)).*weight(i,:)+(weight(i,:)~=tmp(i,fun_num)).*(weight(i,:).*(1-tmp(i,fun_num).^10));
end
weight=weight./repmat(sum(weight,2),1,fun_num);

fit=0;
for i=1:fun_num
    oo=repmat(o(i,:),ps,1);
    eval(['f=feval(func.f' int2str(i) ',((x-oo)./repmat(lamda(i,:),ps,1)).*' int2str(i) ');']);
    x1=5*ones(1,D);
    eval(['f1=feval(func.f' int2str(i) ',(x1./lamda(i,:)).*' int2str(i) ');']);
    %eval(['(x1./lamda(i,:))*M.M' int2str(i)])
    fit1=2000.*f./f1;
    %[fit1 f1 weight(:,i)]
    %weight(:,i).*(fit1+bias(i))
    fit(isinf(fit))=0;
    fit=fit+weight(:,i).*(fit1+bias(i));
end
%-------------------------------------------------
%basic functions

function f=fsphere(x)
%Please notice there is no use to rotate a sphere function, with rotation
%here just for a similar structure as other functions and easy programming
[ps,D]=size(x);
f=sum(x.^2,2);
%--------------------------------
function f=fsphere_noise(x)
[ps,D]=size(x);
f=sum(x.^2,2).*(1+0.1.*normrnd(0,1,ps,1));
%--------------------------------
function f=fgriewank(x)
[ps,D]=size(x);
f=1;
for i=1:D
    f=f.*cos(x(:,i)./sqrt(i));
end
f=sum(x.^2,2)./4000-f+1;
%--------------------------------
function f=fackley(x)
[ps,D]=size(x);
f=sum(x.^2,2);
f=20-20.*exp(-0.2.*sqrt(f./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);
%--------------------------------
function f=frastrigin(x)
[ps,D]=size(x);
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
%--------------------------------
function f=frastrigin_noncont(x)
[ps,D]=size(x);
x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
%--------------------------------
function [f]=fweierstrass(x)
[ps,D]=size(x);
x=x+0.5;
a = 0.5;
b = 3;
kmax = 20;
c1(1:kmax+1) = a.^(0:kmax);
c2(1:kmax+1) = 2*pi*b.^(0:kmax);
f=0;
c=-w(0.5,c1,c2);
for i=1:D
f=f+w(x(:,i)',c1,c2);
end
f=f+c*D;

function y = w(x,c1,c2)
y = zeros(length(x),1);
for k = 1:length(x)
	y(k) = sum(c1 .* cos(c2.*x(:,k)));
end
%--------------------------------
function f=fE_ScafferF6(x)
fhd=@(x) ScafferF6(x);
[ps,D]=size(x);

f=0;
for i=1:(D-1)
    f=f+feval(fhd,(x(:,i:i+1)));
end
    f=f+feval(fhd,x(:,[D,1]));
%--------------------------------    
function f=fE_ScafferF6_noncont(x)
fhd=@(x) ScafferF6(x);
[ps,D]=size(x);
x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
f=0;
for i=1:(D-1)
    f=f+feval(fhd,(x(:,i:i+1)));
end
    f=f+feval(fhd,x(:,[D,1]));
%------------------------------
function f=fEF8F2(x)
[ps,D]=size(x);
f=0;
for i=1:(D-1)
    f=f+F8F2(x(:,[i,i+1]));
end
    f=f+F8F2(x(:,[D,1]));
%--------------------------------
function f=fschwefel_102(x)
[ps,D]=size(x);
f=0;
for i=1:D
    f=f+sum(x(:,1:i),2).^2;
end
%--------------------------------
function f=felliptic(x)
[ps,D]=size(x);
a=1e+6;
f=0;
for i=1:D
f=f+a.^((i-1)/(D-1)).*x(:,i).^2;
end
%--------------------------------
% classical Gram Schmid 
 function [q,r] = cGram_Schmidt (A)
% computes the QR factorization of $A$ via
% classical Gram Schmid 
% 
 [n,m] = size(A); 
 q = A;    
 for j=1:m
     for i=1:j-1 
         r(i,j) = q(:,j)'*q(:,i);
     end
     for i=1:j-1   
       q(:,j) = q(:,j) -  r(i,j)*q(:,i);
     end
     t =  norm(q(:,j),2 ) ;
     q(:,j) = q(:,j) / t ;
     r(j,j) = t  ;
 end
  
function M=rot_matrix(D,c)
  A=normrnd(0,1,D,D);
  P=cGram_Schmidt(A);
  A=normrnd(0,1,D,D);
  Q=cGram_Schmidt(A);
  u=rand(1,D);
  D=c.^((u-min(u))./(max(u)-min(u)));
  D=diag(D);
  M=P*D*Q;
