%                 Not_uniform_Set_Shaping_Theory
%                        Contact info
%                Christian.Schmidt55u@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The program performs the following operations:
% 1) generates a random sequence s with not uniform distribution
% 2) calculate the frequencies of the symbols present in the sequence 
% 3) use this information to calculate the N*H0(s), with H0(s) the emprical
%    entropy of the sequence and N lenght of the sequence
% 4) apply the transform f(s) the new sequence have lenght N+1
% 6) compares the empircal entropy multiplied for lenght N*H0(s) of the 
%    generated sequence s witht he empircal entropy multiplied for 
%    lenght (N+1)*H0(f(s)) of transfromated sequence f(s) of lenght N+1
% 7) repeats all these steps a number of times defined by the parameter history
% 8) display the average values obtained
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Important
% The non-uniform source has this characteristic: it generates N symbols 
% from a alphabet of size ns. Among these symbols there is a most probable 
% symbol defined by the parameter mustp, the other symbols have 
% uniform frequency (1-pmax)/(ns-1).
% 
% the parameters that define the source have these limitations:
%
% N lenght of the sequence >= 40
% ns size ofmthe alphabet >= 20
% pmax probaibility of the most frequnet symbol <= 0.8 and > (1-pmax)/(ns-1)
%
% If you want to use the algorithm for sources with different distributions, 
% please contact us. We are able to create versions that work for any type 
% of non-uniform function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Results for different settings
%
%  Ps=probability with which (N+1)*H0(f(s))<N*H0(s)
%
%  ns= number of symbols
%  N=lenght of the sequence
%  pmax=probaibility of the most frequnet symbol
%
%  ns    N    pmax        ps   Average((N+1)*H0(f(s))<N*H0(s))
%  30   400   0.5         88%                8.0 bit
%  40   400   0.5         89%               10.4 bit
%  50   400   0.5         92%               13.0 bit
%  60   400   0.5         95%               15.2 bit
%
%  As you can see, increasing the symbol number increases the probability Ps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

clear all;

% parameters of the source

ns=40;
len=400;
pmax=0.5;
history=1000;

cs=0;
cs2=0;
totcodel=0;
totinfc=0;
tottinfc=0;
tdife=0;
itnent=0;
infc=0;
tinfc=0;
itinfc=0;
dife=0;
con=0;

for i=1:history
   
 index=0;

 % Genaration of the sequence with a not uniform distribution

 ns2=ns-1;

 symbols=1:(ns);
 prob(1,1:ns2)=(1-pmax)/ns2;
 prob(ns2+1)=pmax;

 seq=randsrc(1,len,[symbols; prob]);

 % the empircal entropy multiplied for lenght N*H0(s) of the generated sequences

 infc=0;

 for i=1:len

  sy=seq(1,i);
  fs=nnz(seq==sy)/len;
  infc=infc-log2(fs);

 end

 % Start trasformation

 mcodel=10000;
 
 nseq=fSST(seq);

 % The new sequence is long nlen=len+1, because in the last position we have added the parameter of the function
 
 nlen=len+1;

 % the empircal entropy multiplied for lenght (N+1)*H0(f(s)) of the transformed sequence f(s) of length nlen=len+1

  tinfc=0;

  for i=1:nlen

   sy=nseq(1,i);
   fs=nnz(nseq==sy)/nlen;
   tinfc=tinfc-log2(fs);

  end 

 if tinfc < infc

   cs2=cs2+1;

 end

 dife=infc-tinfc;

 % We apply the inverse transform and we obtain the initial sequence

 iseq=invfSST(nseq);

 % we check that the obtained sequence is equal to the initial sequence

 flag=isequal(seq,iseq);
 con=con+1;

 if flag == false

   fprintf('Error, sequence not equal to the initial sequence\n');
   con

 end

 totinfc=totinfc+infc;
 tottinfc=tottinfc+tinfc;
 tdife=tdife+dife;

 dife=0;

end

% We calculate the average of the empircal entropy multiplied for lenght N*H0(s) of the generated sequences and The average of the  empircal entropy multiplied for lenght (N+1)*H0(f(s)) of the transformed sequences
  
 medinfc=totinfc/history;
 medcodel=totcodel/history;
 medtinfc=tottinfc/history;
 mdife=tdife/history
 cs2

% We calculate the percentage of sequences where (N+1)*H0(f(s)) < N*H0(s)

 pcs=(cs2/history)*100;

% We display the average values obtained

 fprintf('The average of the empircal entropy multiplied for lenght N*H0(s) of the generated sequences\n');
 medinfc

 fprintf('The average of the  empircal entropy multiplied for lenght (N+1)*H0(f(s)) of the transformed sequences\n');
 medtinfc

 fprintf('Number of sequences where (N+1)*H0(f(s)) < N*H0(s) \n');
 cs2

 fprintf('There is a percentage of %2.0f%% that (N+1)*H0(f(s)) < N*H0(s) \n',pcs);

