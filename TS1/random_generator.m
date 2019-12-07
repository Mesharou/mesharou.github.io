function x=random_generator(n,pdfname,opts)
% number of samples
%n = 1000; 

% uniform distribution for y=CDF(x)
% x is then obtained by inverting the CDF
y = rand(n,1);

switch lower(pdfname)
  case 'cauchy'
  % Cauchy's distribution
  x = tan( pi * (y-.5) );
  %semilogy(abs(x))

  case 'normal'
  % Normal distribution
  x = sqrt(2)*erfinv( 2*y-1 );
  %  plot(x)

  case 'exponential'
  % Exponential
  x = - log( 1-y );
  %  plot(x)

  case {'chi-squared','chi2'}
  % Chi-squared
  k = opts;
  x = 2 * gammaincinv(y,k/2);
  %plot(x)
  
  case 'pareto'
  % Pareto
  alpha = opts;
  x = ( 1 - y ).^(-1/alpha);
  %plot(x)
 otherwise
  % uniform
  x = y;
end
  
return

db=0.1;
b=0:db:1e2;

h=hist(x,b);
loglog(b,h)
