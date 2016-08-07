m=6;
n=61;
S0=100;
S1=130;
K=115;
r=0.05;
T=1.0;
sigma=0.3;

f = inline('(-x*x*0.5)/sqrt(2*pi)', 'x');


time = transpose(linspace(T,0,m));
S=linspace(S0,S1,n);

d1 =(log(S/K) + ( (r + sigma ^2/2)*(T-time)))./( sigma*
sqrt(T-time));
d2 = (log(S/K) + ( (r - sigma ^2/2)*(T-time)))./( sigma*
sqrt(T-time));

part1 = bsxfun(@times, normcdf(d1),S);
part2 = bsxfun( @times , K*exp(-r*(T-time)), normcdf(d2));
VC= part1-part2;
Delta = normcdf(d1);

g1 = 1./bsxfun ( @times , sqrt (2*pi)*sigma*sqrt(T-time
), S);
g2 = exp(-d1 .^2/2);
Gamma = g1 .* g2;



subplot(1,3,1);
plot(S,VC);
title("Asset-Option prices, sigma = 0.3")
xlabel('Asset price');
ylabel('Option Value');


subplot(1,3,2);
plot(S,Delta);
title("Delta")
xlabel('Asset price');
ylabel('Delta');

subplot(1,3,3);
plot(S,Gamma);
title("Gamma");
xlabel('Asset price');
ylabel('Gamma');
 
