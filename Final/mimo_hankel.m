function [H0,H1] = mimo_hankel(MP,s)
%MIMO_HANKEL Summary of this function goes here
%   Detailed explanation goes here
[n_output,n_input,l] = size(MP);

if (~exist('s', 'var'))
        s = floor(l/3);
end


H0 = zeros(n_output*s,n_input*s);
H1 = zeros(n_output*s,n_input*s);
for i=1:s
    for j=1:l-s-1
        H0(n_output*(i-1)+1:n_output*(i),n_input*(j-1)+1:n_input*j) = MP(:,:,i+j);
        H1(n_output*(i-1)+1:n_output*(i),n_input*(j-1)+1:n_input*j) = MP(:,:,i+j+1);
    end
end

end

