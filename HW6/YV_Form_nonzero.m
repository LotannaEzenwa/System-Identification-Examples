function [ Y, V ] = YV_Form_nonzero( input_data, output_data, p )
% YV_Form
% Calculates the Y and V Matrix for ARX models

[r, k1] = size(input_data);
[m, k2] = size(output_data);
k = k2;
l = k;
s = (m+r)*p+r;

y = zeros(m,l);
V_p = zeros(s,l);
v_p = zeros(s,1);


if k1~=k2
    error('Input and Output Size Mismatch')
end

for i=1:k
    for j=1:p
        
        if i-j < 1
            v_p(r+j*(r+m):end) = 0;
            break
        end
        
        
        ind1 = i-j;
        b = length(1+r+j*(r+m):1+r+(j+1)*(r+m)-1);
        
        v_p(r+(j-1)*b+1:r+j*b) = [input_data(:,ind1);output_data(:,ind1)];
        %v_p(2*(i-j-p+1):2*(i-j-p+1)+1) = [output_data(i-j);input_data(i-j)];
        
    end
   
    V_p(:,i) = v_p;
    
    V_p(1:r,i) = input_data(:,i);
    
    v_p(:) = 0;
    
    
end
Y = output_data(:,p+1:end)';
V = V_p(:,p+1:end);


end

