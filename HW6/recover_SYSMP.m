function [RMP] = recover_SYSMP(MP,n,r,m)
% RECOVER_SYSMP 
%   Given Markov Parameters, n, and the inputs and outputs, this function
%   calculates n System Markov Parameteres

[a,k] = size(MP);
p = (k-1)/(r+m);
D = MP(:,1);
if n<2
    RMP = [D];
end

if k == 1
    error('Chain Not Long Enough')
    return
end
% if mod(k,2) == 0
%     error('Include all Markov Parameters')
%     return
% end

RMP = zeros(r,n);
RMP(:,1) = D;

    function [idx] = observer_MP_index(mp_k, expo,inn,oot)
        idx = 2;
        strt = (inn+oot)*(mp_k-1);
        idx = idx + strt;
        if expo==1
            idx = idx:idx+inn-1;
        elseif expo==2
            idx = idx+inn:idx+inn+oot-1;
        else
            idx = idx:idx+inn+oot-1;
        end
        
        
    end


for kk=1:n-1
    
    S = 0;
    if kk<p+1
        for i=1:kk
            
            idx = observer_MP_index(i,2,m,r);
            S = S + -MP(:,idx)*RMP(:,kk-i+1);
        end
        id1 = observer_MP_index(kk,1,m,r);
        RMP(:,kk+1) = MP(:,id1) - S;
        
    else
        
        for i=1:p
        
            idx = observer_MP_index(i,2,m,r);
            S = S + -MP(:,idx)*RMP(:,kk-i+1);
            
        end
        
        RMP(:,kk+1) = -S;
    end
    
end
end

