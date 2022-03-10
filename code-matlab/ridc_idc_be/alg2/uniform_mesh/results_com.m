  % % p=4,      % % example 1 % %

 u_my=2.718281827021971;  % n=100
 u_ex=2.718281828459046; % exact
 u_sf=2.718281826174;     % ridc soft for n=100

 u_my =2.718281828374721;  % n=200 
 u_ex= 2.718281828459046; % exact
 u_sf =2.718281828314;    % ridc soft for n=200

 u_my =2.718281828453943;  % n=400 
 u_ex= 2.718281828459046; % exact
 u_sf =2.718281828450;  % ridc soft for n=400
 
 %%  example 2 %%   p=4  % 
 
 u_my=1.000000189400767;  % n=100
 u_ex=1;                  % exact
 u_sf=1.000000227105;     % ridc soft for n=100

 u_my =1.000000011776812;  % n=200 
 u_ex= 1;                 % exact
 u_sf =1.000000013987;    % ridc soft for n=200
 
 u_my =1.000000000734351;  % n=400 
 u_ex= 1;                  % exact
 u_sf =1.000000000867;  % ridc soft for n=400
 
 % order test of ridc soft for p=4, example 1

e1=exp(1)-2.718281826174;
e2=exp(1)-2.718281828314;
e3=exp(1)-2.718281828450;

p1=log2(e1/e2)
p2=log2(e2/e3)

% order test of ridc soft for p=4, example 2

er1=1.000000227105-1;
er2=1.000000013987-1;
er3=1.000000000867-1;

p3=log2(er1/er2)
p4=log2(er2/er3)
