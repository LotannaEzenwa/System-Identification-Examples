function [RMP] = recover_SYSMP(MP,n,m,r)
% RECOVER_SYSMP 
%   Given Markov Parameters, n, and the inputs and outputs, this function
%   calculates n System Markov Parameteres

[a,k] = size(MP);
D = MP(:,1:r);
MP = reshape(MP(:,r+1:end),m,r+m,[]);
[qq,qq,p] = size(MP);


if k == 1
    error('Chain Not Long Enough')
    return
end
% if mod(k,2) == 0
%     error('Include all Markov Parameters')
%     return
% end

RMP = zeros(m,r,n);
RMP(:,:,1) = D;
for k=2:n
    if k<=p+1
        S = 0;
        for i=1:k-1
            
            if k==i
                S = S + -MP(:,r+1:end,i)*D;
            else
                S = S + -MP(:,r+1:end,i)*RMP(:,:,k-i);
            end
        end
        RMP(:,:,k) = MP(:,1:r,k-1) - S;
    else
        S = 0;
        for i=1:p
                S = S + -MP(:,r+1:end,i)*RMP(:,:,k-i);
         
        end
        RMP(:,:,k) = -S;
    end
end
     

 
end

