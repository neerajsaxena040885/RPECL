% Calculates log-likelihood function value for mixed logit model
% Written by Kenneth Train, July 27, 2006, revised July 31, 2006
%
% This code is input to Matlab's funcmin command
%
% Input param is a column vector of parameters, dimension (NF+NV+NV)x1
%     containing the fixed coefficients, the first parameters of the random
%     coefficients, and then the second parameters of the random coefficients
% Output ll is the scalar value of the negative of the simulated log-likelihood 
%     at the input parameters

function [ll, g] =loglik(param)

global NV NF IDV NEC

if NF>0
  f=param(1:NF,1);
else
  f=[];
end

if NV>0
  if sum(IDV(:,2) == 5) >0;
     b=zeros(NV,1);
     b(IDV(:,2) ~= 5,1)=param(NF+1:NF+sum(IDV(:,2) ~= 5),1);
     w=param(NF+sum(IDV(:,2) ~= 5)+1:end,1);
  else;
     b=param(NF+1:NF+NV,1);
     w=param(NF+NV+1:NF+NV+NV,1);
  end;
else
  b=[];
  w=[];
end

if NEC>0
  ec=param(NF+NV+NV+NEC,1);
else
  ec=[];
end

[p g]=llgrad2(f,b,w,ec); 


ll=-sum(log(p),2);
g=-sum(g,2);



