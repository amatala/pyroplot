function Y = filtNs(X,N)
% FILTNS        Symmetric filter for noisy data
%
% Y = filtNs(x,N)
%
% Symmetric N-point averaging. N must be odd integer.
%
% See also: filtN
%
if mod(N,2) == 0,
   error('N must be an odd integer.')
end
M = length(X);
Y(1,:) = X(1,:);
for i = 2:(N-1)/2
   Y(i,:) = sum(X(1:i+i-1,:))/(2*i-1);
end
for i = 1+(N-1)/2 : M-(N-1)/2 
   Y(i,:) = sum(X(i-(N-1)/2:i+(N-1)/2,:))/N;
end
for i = M-(N-1)/2+1:M-1
   Y(i,:) = sum(X(i-(M-i):M,:))/(1+2*(M-i));
end
Y(M,:) = X(M,:);
