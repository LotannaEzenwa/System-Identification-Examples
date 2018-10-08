function [ Y, V ] = YV_Form_zero_test( input_data, output_data, p)
% YV_Form
% Calculates the Y and V Matrix for ARX models


[k1, r] = size(input_data);
[k2, m] = size(output_data);

Y = zeros(m,k1);
V = 0;
V_lower = zeros(r+m,p,k1);



k = k1;
V_query =[input_data(:,:),output_data(:,:)].';
for i=0:k-1
    if i<p
        idx = i;
        if idx<1
            continue
        end
        vv = flip(V_query(:,1:i),2);
        V_lower(:,1:i,i+1) = reshape(vv(:),r+m,[]);
    
    else
        vv = flip(V_query(:,i-p+1:i),2);
        V_lower(:,:,i+1) = reshape(vv(:),r+m,[]);
        
    end
    
    
    
end
V_ll = reshape(V_lower,p*(r+m),[]);
V = [input_data(:,:)';V_ll];
V = V(:,p+1:end);
Y = output_data(:,p+1:end)';


end

