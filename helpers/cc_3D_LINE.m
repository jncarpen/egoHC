function [ CCM ] = cc_3D_LINE( A,B,method )
%%CC_3D 3D Matrix Cross Correlation
% Calclates the 3D Cross Correlation for matrices A and B
%%Inputs
%{
A: 1st Matrix for Cross Correlation
B: 2nd Cross Correlation
method: Method for the cross correlation
   (1): 'raw': Nomal Crosscorrelation   
%}
%%Outputs
%{
CCM: Cross Correlation Matrix
%}
%%Check and Patch the matrices
as=size(A);
bs=size(B);
if bs>=as  %Matrix B larger than A
    A = padarray(A,(bs-as)/2); % Pad the arrays to the same size
elseif as>=bs
    B = padarray(B,(as-bs)/2); % Pad the arrays to the same size
end
% if the matrices have the same size, nothing has to be done
as=size(A);
CCM=zeros(2*as(1)-1,2*as(2)-1,2*as(3)-1);
%%Cross-Correlation Calculation  
% Data Loops for all the values to be calculated, given 3D matrix with
% size x:size(A,1)*2-1; y:size(A,2)*2-1; z:size(A,3)*2-1 [cc].
% 
if strcmp(method,'raw');
for ii=1:((as(1)*2-1)*(as(2)*2-1)*(as(3)*2-1));
      [i,j,k]=ind2sub(size(CCM),ii);
      %%Search for the summation boundary values
  %%%%%%%%%%%%%%%% x Values
      if i<=as(1)
          axi=1;          axt=i;
          bxi=as(1)-i+1;    bxt=as(1);
      elseif i>as(1)
          axi=i-as(1)+1;    axt=as(1);
          bxi=1;          bxt=2*as(1)-i;
      end
  %%%%%%%%%%%%%%%%  y Values
      if j<=as(2)
          ayi=1;            ayt=j;
          byi=as(2)-j+1;    byt=as(2);
      elseif j>as(2)
          ayi=j-as(2)+1;    ayt=as(2);
          byi=1;            byt=2*as(2)-j;
      end                
  %%%%%%%%%%%%%%%%  z Values
      if k<=as(3)
          azi=1;            azt=k;
          bzi=as(3)-k+1;    bzt=as(3);
      elseif k>as(3)
          azi=k-as(3)+1;    azt=as(3);
          bzi=1;            bzt=2*as(3)-k;
      end        
      %%Cross Correlation Operation  
      % First, Calculate the vector product of the two matrices
      CDummy=...
      A(axi:axt,ayi:ayt,azi:azt).*B(bxi:bxt,byi:byt,bzi:bzt);
      CCM(i,j,k)=sum(CDummy(:));
end
end
end