k = 6;

n = 1000000;
cos_array = zeros(1,k);
sin_array = zeros(1,k);
l = 50;
all_distance = zeros(1,n);

for i = 1:n
    random_array = rand(1, k)*2*pi;

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
    
    all_distance(1,i) = sqrt(l^2 * (cos_sum^2 + sin_sum^2));
    %disp(all_distance)
end

% Define the bin size (width)
binSize = 1; % Adjust this to your desired bin size

% Create a histogram plot
histogram(all_distance, 'BinWidth', binSize);

% Customize the plot (optional)
title('Histogram with Bin Size 2');
xlabel('X-axis Label');
ylabel('Frequency');
