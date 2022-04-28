%GAMMA MIX - 4 COMPONENTSto obtain Table 3
%--------------------------------------------------------
function [out] = gamma4_mixture(k_P)
 %this function generates the posterior draws for mixture of gamma using MCMC method
 load('data_dominance_paper.mat');% load the dataset
 x = income_1999_perm(:,1);
 x = x(1:round(10/10*length(x)));
 X1=x(1:3/10*length(x));
 X2=x(4/10*length(x):6/10*length(x));
 X3=x(6/10*length(x):8/10*length(x));
 X4=x(8/10*length(x):length(x));
 n=size(x,1);
 %Preliminaries
 k=4;		%number of components
 
%Number of replications
 b=1000;	%number of burn-in replications
 nit=10000;	%number of retained replications / sample size
 s=b+nit;	%total

%Store draws after burn-in in the following matrices 
 W=zeros(nit,4);
 M=zeros(nit,k);
 V=zeros(nit,k);
 totalz=zeros(n,1);
 %Including burn-in draws
 Wtot=zeros(s,4);
 Mtot=zeros(s,k);
 Vtot=zeros(s,k);
  
%Priors
 phi1=1;
 phi2=1;
 phi3=1;
 phi4=1;
%  theta1=(mean(x)^2)/var(x);
%  theta2=(mean(x)^2)/var(x);
%  theta3=(mean(x)^2)/var(x);
%  theta4=(mean(x)^2)/var(x);
 theta1=0.01;
 theta2=0.01;
 theta3=0.01;
 theta4=0.01;
 alpha1=2.2;
 alpha2=2.2;
 alpha3=2.2;
 alpha4=2.2;
 beta1=60;
 beta2=100;
 beta3=150;
 beta4=200;
% %Step 1 -> [t=0, set initial values w0, m0, v0]
  %Initial value for w (weight):
   w1=0.05;
   w2=0.42;
   w3=0.42;
   w4=1-w1-w2-w3;
  %Initial value for m (mean):
  m1 = 180;
  m2 = 395;
  m3 = 445;
  m4 = 4334.8;
  %Initial value for v (shape parameter):
   v1=37;
   v2=0.96;
   v3=5.348;
   v4=186.88;

   %defining logs for logPz
   logx=log(x);

   accept1 = 0;
   accept2 = 0;
   accept3 = 0;
   accept4=0;
for i=1:s
    i
%Step 2 - Draw z from p(z|x,k,wt,mt,vt)
 %z=w*G(v1,m1,x); %(proportional to)
 Z1=w1*gampdf(x,v1,m1/v1);
 Z2=w2*gampdf(x,v2,m2/v2);
 Z3=w3*gampdf(x,v3,m3/v3);
 Z4=w4*gampdf(x,v4,m4/v4);
 PZ1=Z1./(Z1+Z2+Z3+Z4);
 PZ2=Z2./(Z1+Z2+Z3+Z4);
 PZ3=Z3./(Z1+Z2+Z3+Z4);
 PZ4=Z4./(Z1+Z2+Z3+Z4);
 Z=rand(n,1);
 
 z1 = (Z <=PZ1);
 z2 = (Z > PZ1) & (Z <= PZ1+PZ2);
 z3 = (Z > PZ1+PZ2) & (Z <= PZ1 + PZ2 + PZ3);
 logx=log(x);

 N1=sum(z1);
 N2=sum(z2);
 N3=sum(z3);
 N4=sum((1-z1-z2-z3));
 S1=sum(x.*(z1));
 S2=sum(x.*(z2));
 S3=sum(x.*(z3));
 S4=sum(x.*((1-z1-z2-z3)));
 logP1=sum(logx.*(z1));
 logP2=sum(logx.*(z2));
 logP3=sum(logx.*(z3));
 logP4=sum(logx.*((1-z1-z2-z3)));
        
%Step 3 - Draw w from p(w|x,k,zt+1,mt,vt)
 weight = dirich_rnd([phi1+N1 ; phi2+N2 ; phi3+N3 ; phi4+N4]);
 w1 = weight(1,1);
 w2 = weight(2,1);
 w3 = weight(3,1);
 w4 = weight(4,1);
%Step 4 - Draw m from p(m|x,k,zt+1,wt+1,vt) for z=1,2
  m1=1/(gamrnd(alpha1+N1*v1,1/(beta1+S1*v1)));
  m2=1/(gamrnd(alpha2+N2*v2,1/(beta2+S2*v2)));
  m3=1/(gamrnd(alpha3+N3*v3,1/(beta3+S3*v3)));
  m4=1/(gamrnd(alpha4+N4*v4,1/(beta4+S4*v4)));
  
%Step 5 - Draw v from p(v|x,k,zt+1,wt+1,mt+1) for z=1,2
	%Metropolis step
	 r1=100;
     r2=100;
     r3=100;
     r4=100;
     A1=rand(1);
     A2=rand(1);
     A3=rand(1);
     A4=rand(1);
	%Draw candidate V from G(r,r/v) for z=1,2
	 V1=gamrnd(r1,v1/r1);
     V2=gamrnd(r2,v2/r2);
     V3=gamrnd(r3,v3/r3);
     V4=gamrnd(r4,v4/r4);
	%Calculate posterior to accept/reject
        logvv1=v1*N1*log(v1)-N1*gammaln(v1)-v1*(theta1+(S1/m1)+N1*log(m1)-logP1);
        logvc1=V1*N1*log(V1)-N1*gammaln(V1)-V1*(theta1+(S1/m1)+N1*log(m1)-logP1);
        logV_v1=log(gampdf(v1,r1,V1/r1));
        logv_V1=log(gampdf(V1,r1,v1/r1));
        MH1=exp(logvc1+logV_v1-logvv1-logv_V1);
        C1=min(1,MH1);
            if A1<=C1
            v1=V1;
            accept1 = accept1+1;
            end
        logvv2=v2*N2*log(v2)-N2*gammaln(v2)-v2*(theta2+(S2/m2)+N2*log(m2)-logP2);
        logvc2=V2*N2*log(V2)-N2*gammaln(V2)-V2*(theta2+(S2/m2)+N2*log(m2)-logP2);
        logV_v2=log(gampdf(v2,r2,V2/r2));
        logv_V2=log(gampdf(V2,r2,v2/r2));
        MH2=exp(logvc2+logV_v2-logvv2-logv_V2);
        C2=min(1,MH2);
            if A2<=C2
            v2=V2;
            accept2 = accept2+1;
            end
        logvv3=v3*N3*log(v3)-N3*gammaln(v3)-v3*(theta3+(S3/m3)+N3*log(m3)-logP3);
        logvc3=V3*N3*log(V3)-N3*gammaln(V3)-V3*(theta3+(S3/m3)+N3*log(m3)-logP3);
        logV_v3=log(gampdf(v3,r3,V3/r3));
        logv_V3=log(gampdf(V3,r3,v3/r3));
        MH3=exp(logvc3+logV_v3-logvv3-logv_V3);
        C3=min(1,MH3);
            if A3<=C3
            v3=V3;
            accept3 = accept3+1;
            end
        logvv4=v4*N4*log(v4)-N4*gammaln(v4)-v4*(theta4+(S4/m4)+N4*log(m4)-logP4);
        logvc4=V4*N4*log(V4)-N4*gammaln(V4)-V4*(theta4+(S4/m4)+N4*log(m4)-logP4);
        logV_v4=log(gampdf(v4,r4,V4/r4));
        logv_V4=log(gampdf(V4,r4,v4/r4));
        MH4=exp(logvc4+logV_v4-logvv4-logv_V4);
        C4=min(1,MH4);
            if A4<=C4
            v4=V4;
            accept4 = accept4+1;
            end
                           
%Order resriction (optional)
if m2>m3
qq=m3;
m3=m2;
m2=qq;
qq=v3;
v3=v2;
v2=qq;
qq=w3;
w3=w2;
w2=qq;
end

if m3>m4
qq=m4;
m4=m3;
m3=qq;
qq=v4;
v4=v3;
v3=qq;
qq=w4;
w4=w3;
w3=qq;
end


if m1>m2
qq=m2;
m2=m1;
m1=qq;
qq=v2;
v2=v1;
v1=qq;
qq=w2;
w2=w1;
w1=qq;
end

%Storing draws
    if i>b
      %Store all draws after discarding burn-in
      W(i-b,1)=w1;
      W(i-b,2)=w2;
      W(i-b,3)=w3;
      W(i-b,4)=w4;
      M(i-b,1)=m1;
      M(i-b,2)=m2;
      M(i-b,3)=m3;
      M(i-b,4)=m4;
      V(i-b,1)=v1;
      V(i-b,2)=v2;
      V(i-b,3)=v3;
      V(i-b,4)=v4;
      
      %totalz=totalz+z;
    end
      %Store all draws including burn-in
      Wtot(i,1)=w1;
      Wtot(i,2)=w2;
      Wtot(i,3)=w3;
      Wtot(i,4)=w4;
      Mtot(i,1)=m1;
      Mtot(i,2)=m2;
      Mtot(i,3)=m3;
      Mtot(i,4)=m4;
      Vtot(i,1)=v1;
      Vtot(i,2)=v2;   
      Vtot(i,3)=v3;
      Vtot(i,4)=v4;









end

%Analysis
prop1 = accept1/s;
prop2 = accept2/s;
prop3 = accept3/s;
prop4 = accept4/s;

mW2 = mean(W(:,2))
stdW2 = std(W(:,2));
minW2 = min(W(:,2));
maxW2 = max(W(:,2));


mW1=mean(W(:,1))
stdW = std(W(:,1));
minW = min(W(:,1));
maxW = max(W(:,1));


mW3 = mean(W(:,3))
mM1=mean(M(:,1))
stdM1 = std(M(:,1));
minM1 = min(M(:,1));
maxM1 = max(M(:,1));

mM2=mean(M(:,2))
stdM2 = std(M(:,2));
minM2 = min(M(:,2));
maxM2 = max(M(:,2));

mM3 = mean(M(:,3));
mM4 = mean(M(:,4));
mV1=mean(V(:,1))
stdV1 = std(V(:,1));
minV1 = min(V(:,1));
maxV1 = max(V(:,1));
 
mV2=mean(V(:,2))
stdV2 = std(V(:,2));
minV2 = min(V(:,2));
maxV2 = max(V(:,2));

mV3 = mean(V(:,3))
mV4 = mean(V(:,4))
 subplot(3,5,1); plot(W(:,1));title('W1 draws');
 subplot(3,5,2); plot(W(:,2));title('W2 draws');
 subplot(3,5,3); plot(W(:,3));title('W3 draws');
 subplot(3,5,4); plot(M(:,1));title('M1 draws');
 subplot(3,5,5); plot(M(:,2));title('M2 draws');
 subplot(3,5,6); plot(M(:,3));title('M3 draws');
 subplot(3,5,7); plot(V(:,1));title('V1 draws');
 subplot(3,5,8); plot(V(:,2));title('V2 draws');
 subplot(3,5,9); plot(V(:,3));title('V3 draws');
W1_1999 = W(:,1);
W2_1999 = W(:,2);
W3_1999 = W(:,3);
W4_1999 = W(:,4);
M1_1999 = M(:,1);
M2_1999 = M(:,2);
M3_1999 = M(:,3);
M4_1999 = M(:,4);
V1_1999 = V(:,1);
V2_1999 = V(:,2);
V3_1999 = V(:,3);
V4_1999 = V(:,4);

nm = ['output_4gamma_1999',k_P,'.mat']
save(nm,'W1_1999','W2_1999','W3_1999','W4_1999','M1_1999','M2_1999','M3_1999','M4_1999','V1_1999','V2_1999','V3_1999','V4_1999')

end
