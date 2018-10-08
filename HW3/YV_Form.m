function [ Y, V ] = YV_Form( input_data, output_data, p )
% YV_Form
% Calculates the Y and V Matrix for ARX models

[r, k1] = size(input_data);
[m, k2] = size(output_data);
k = k2-1;
s = (m+r)*p+r;

V_p = zeros(s,k+1);
v_p = zeros(p*(r+m),1);


if k1~=k2
    error('Input and Output Size Mismatch')
end

for i=1:k1
    for j=1:p
        
        if i-j < 1
            v_p(2*j:end) = 0;
            break
        end
        
        
        ind1 = i-j;
        
        v_p(2*j-1:2*j) = [output_data(:,ind1);input_data(:,ind1)];
        %v_p(2*(i-j-p+1):2*(i-j-p+1)+1) = [output_data(i-j);input_data(i-j)];
        
    end

    V_p(1,i) = input_data(:,i);
    V_p(2:end,i) = v_p;
    v_p(:) = 0;
    
    
end
Y = output_data(:)';
V = V_p;


end

