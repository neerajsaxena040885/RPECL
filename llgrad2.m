% Calculate logit probability for chosen alternatives for each person
%    with multiple choice situations for each person and multiple people,
%    using globals for all inputs except coefficients.
% Written by Kenneth Train, first version on July 14, 2006..
%    lastest edits on Sept 24, 2006.
%
% Simulated Mixed Logit Probability and gradient
% Logit probability is Prod_t { exp(V_*t)/SUM_j[exp(V_jt)] } 
%             = Prod_t { 1 / (1+ Sum_j~=* [exp(V_jt-V-*t)] }
%    where * denotes chosen alternative, j is alternative, t is choice situation.
% Using differences from chosen alternative reduces computation time (with
%    one less alternative), eliminates need to retain and use the dependent
%    variable, and avoids numerical problems when exp(V) is very large or
%    small, since prob is 1/(1+k) which can be evaluated for any k, including
%    infinite k and infinitesimally small k. In contrast. e(V)/sum e(V) can
%    result in numerical "divide by zero" if denominator is sufficiently small
%    and NaN if numerator and denominator are both numerically zero.
% 
% Input f contains the fixed coefficients, and has dimension NFX1.
% Input c contains the random coefficients for each person, and has dimensions NV x NP.
% Either input can be an empty matrix. 
% Output p contains the logit probabilities, which is a row vector of dimension 1xNP.
% Output g contains the gradients of log(p), which is a matrix (NF+NV+NV) x NP
% Code assumes that all GLOBALS already exist.


function [p, g]=llgrad2(f,b,w,ec); 

global NV NF NP NDRAWS NALTMAX NCSMAX NMEM NEC
global X S XF DR 
global MDR ECVECTOR

p=zeros(1,NP);
g=zeros(NF+NV+NV+NEC,NP);

c=trans(b,w,DR);   %Transforms draws into random coefficients c is NV x NP x NMEM
v=zeros(NMEM,NALTMAX-1,NCSMAX,NP);
if NF > 0
   ff=reshape(f,1,1,NF,1);
   ff=repmat(ff,[NALTMAX-1,NCSMAX,1,NP]);
   vf=reshape(sum(XF.*ff,3),NALTMAX-1,NCSMAX,NP);  %vf is (NALTMAX-1) x NCSMAX x NP
else
   vf=zeros(NALTMAX-1,NCSMAX,NP);
end
vf=repmat(vf,[1,1,1,NMEM]);

if NV >0
   
   cc=reshape(c,1,1,NV,NP,NMEM);
   cc=repmat(cc,[NALTMAX-1,NCSMAX,1,1,1]);
   v=(repmat(X,[1,1,1,1,NMEM]).*cc); %v is (NALTMAX-1) x NCSMAX x NV x NP x NMEM
   v=reshape(sum(v,3),NALTMAX-1,NCSMAX,NP,NMEM);             %v is (NALTMAX-1) x NCSMAX x NP x NMEM
   v=v+vf; 
else
   v=vf;
end

% 16/08/2016: Code change to accommodate Error components (ECVECTOR) in the Utility
% equations
if NEC > 0
   v=v+ec.*ECVECTOR; 
end

v=exp(v);
v(isinf(v))=10.^20;  %As precaution when exp(v) is too large for machine
v=v.*repmat(S,[1,1,1,NMEM]);
pp=1./(1+sum(v,1)); %pp is 1 x NCSMAX x NP x NMEM

%Calculate gradient -- Using Simulated Gradient as of now...


%Back to prob
pp=reshape(pp,NCSMAX,NP,NMEM);     %pp is now NCSMAX x NP x NMEM
pp=prod(pp,1);      %pp is 1xNPxNMEM
pp=sum(pp,3);       %pp is 1xNP

p=p+pp;

p=p./NDRAWS;
p(1,isnan(p))=1; %Change missing values to 1, as a precaution. 
%Gradient
   g=g./NDRAWS;
   g=g./repmat(p,NF+NV+NV+NEC,1);