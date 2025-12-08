function index = refineCodePhase(correlation)
%Functions refines the code phase estimation in the acquisition
%
%index = refineCodePhase(correlation)
%
%   Inputs: 
%       correlation  - The correlation results of the acquired signal
%
%   Outputs:
%       index        - The refined index of the peak

leng=length(correlation);
[prompt,indP]=max(correlation);

indE=indP-1;indVE=indP-2;
indL=indP+1;indVL=indP+2;

% Early Correlator
temp=rem(indE+leng,leng);
if temp==0
    temp=leng;
end
E=correlation(temp);
% Very Early Correlator
temp=rem(indVE+leng,leng);
if temp==0
    temp=leng;
end
VE=correlation(temp);
% Late Correlator
temp=rem(indL,leng);
L=correlation(temp);
% Very Early Correlator
temp=rem(indVL,leng);
VL=correlation(temp);

if E>=L
    a1=indE;b1=E;a2=indP;b2=prompt;
    a3=indVE;b3=VE;a4=indL;b4=L;
else
    a1=indP;b1=prompt;a2=indL;b2=L; 
    a3=indE;b3=E;a4=indVL;b4=VL;
end
hL=(b1-b3)/(a1-a3);
hR=(b2-b4)/(a2-a4);
index=(hL*a1-hR*a2+b2-b1)/(hL-hR);
index=index-1;
