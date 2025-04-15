a = find_alpha(3);
disp(a);

function alpha = find_alpha(k)
    n = 1000000;
    cos_array = zeros(1,k);
    sin_array = zeros(1,k);
    l = 50;
    r = 5;
    connect = 0;
    not_connect = 0;
    
    for i = 1:n
        random_array = rand(1, k)*2*pi;
    
        SUM = 0;
        cos_sum = 0;
        sin_sum = 0;
        
        for j = 1:k
            cos_array(j) = cos(random_array(j));
            sin_array(j) = sin(random_array(j));  
        end
        
        for o = 1:k
            cos_sum = cos_sum + cos_array(o);
            sin_sum = sin_sum + sin_array(o);
            
        end
        %disp(cos_sum)
    
        SUM = l^2 * (cos_sum^2 + sin_sum^2);
        
        if SUM <= 4*r^2
            connect = connect + 1;
            %disp(SUM)
            %disp(connect)
        else
            not_connect = not_connect + 1;
            %disp(SUM)
        end
    end
    alpha = 2*connect/(connect+not_connect);
end

