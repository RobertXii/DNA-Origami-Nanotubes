initialTemp = 326;
ramp = 1;
finalTemp = 293;

dt = 1;
total_t = 3600*(initialTemp-finalTemp)/ramp+3600*0;
disp(total_t/3600)
mu = 0.06;
H = -114.5*4184.0 / 8.3;
S = -313.5*4.184 / 8.3;

kinds = 12;
dpm = ramp/60;
n = round(total_t/dt);
npm = 60/dt;
o = 0;
t = 1;

mol_m = 15 * 10^-9;%mol of monomers
m_nm3 = mol_m * 0.1;%1nM = 10^-6 * 6*10^23/10^27 = 6*10^-10 nm^-3
m = mol_m * 6.02 * 10^23 * 10^3;

%disp(origamis)
origamis = zeros(n, kinds);
const = zeros(kinds,kinds);
form_loop = zeros(1,kinds);
origamis(1,1) = m;
loop = zeros(n,kinds);
Additional_num = round((initialTemp-finalTemp)*60/ramp);
summation = zeros(n+Additional_num,1);
summation(1,1) = m;

%define alpha
for p = 1:kinds
    if mod(p, 3) == 0
        form_loop(1,p) = find_alpha(p);
    else
        form_loop(1,p) = 0;
    end
end


for r = 1:n+Additional_num
    if mod(o, npm) == 0 && initialTemp - finalTemp >= 0.1 %change the kernal every 1min
        initialTemp = initialTemp - dpm;
        M = (initialTemp + 273) * 2 * 1.38 * 10^-23 * melt(initialTemp, H, S, m_nm3)^2 / (3 * mu);
        for j = 1:kinds
            for k = 1:kinds
                const(j,k) = M*(j+k)^2 / (3*k*j);
            end
        end
        o = o+1;
        %disp(initialTemp);
        %disp(melt(initialTemp, H, S, m_nm3));

    else
        dori = zeros(1,kinds);
        %fraction = melt(initialTemp, H, S, m_nm3)^2;
        for i = 1:kinds 
            for j = 1:kinds %rate of loosing
                if j == i 
                    dori(1,i) = dori(1,i) - 2*origamis(t,i)*origamis(t,j)*const(i,j)  ;
                else
                    dori(1,i) = dori(1,i) - origamis(t,i)*origamis(t,j)*const(i,j);
                end
            end
            dori(1,i) = dori(1,i) - form_loop(1,i)*origamis(t,i); %from string to loop

            add = 0;
            for k = 1:i-1 %rate of getting
                if i == 1
                    add = add + 0;
                elseif i-k == k
                    add = add + origamis(t,k)*origamis(t,i-k)*const(k,i-k);
                else 
                    add = add + 0.5*origamis(t,k)*origamis(t,i-k)*const(k,i-k);
                end
            end
            %disp(melt(initialTemp, H, S, m_nm3))

            loop(t+1,i) = loop(t,i) + form_loop(1,i)*origamis(t,i)*dt; %loop incraesed
            dori(1,i) = dori(1,i) + add;
            origamis(t+1,i) = origamis(t,i) + dt*dori(1,i);
        end

        for c = 1:kinds
            summation(t+1,1) = summation(t+1,1) + origamis(t+1,c)*c+loop(t+1,c)*c;
        end
        t = t+1;
        o = o+1;
    end
end

l = length(transpose(origamis(:, 1)));

%create a new matrix include all
all = zeros(t,kinds+round(kinds/3));
e = 0;
for i = 1:kinds+round(kinds/3)
    if mod(i+1,4) == 0
        all(:,i) = loop(:,i-e)*(i-e);
        e = e+1;
    else 
        all(:,i) = origamis(:,i-e)*(i-e);
    end
end

disp(all(t,3)/all(t,7))


time = 0:dt:(l-1)*dt;
plot(time, transpose(all(:, 1)), 'LineWidth', 2, 'Color', '#ff1900');  % Line for y1 will be blue
hold on;  % This command allows you to add additional lines to the same plot
plot(time, transpose(all(:, 2)), 'LineWidth', 2, 'Color', '#ff8c00');   % Line for y2 will be red
plot(time, transpose(all(:, 3)), 'LineWidth', 2, 'Color', '#ffb820'); % Line for y3 will be green
plot(time, transpose(all(:, 4)), 'LineWidth', 2, 'Color', '#1aff00');
plot(time, transpose(all(:, 5)), 'LineWidth', 2, 'Color', '#00ffe5');
plot(time, transpose(all(:, 6)), 'LineWidth', 2, 'Color', '#0084ff');
plot(time, transpose(all(:, 7)), 'LineWidth', 2, 'Color', '#0400ff');
plot(time, transpose(all(:, 8)), 'LineWidth', 2, 'Color', '#8800ff');
plot(time, transpose(all(:, 9)), 'LineWidth', 2, 'Color', '#ff00e6');
plot(time, transpose(all(:, 10)), 'LineWidth', 2, 'Color', '#0008ff');
plot(time, transpose(all(:, 11)), 'LineWidth', 2, 'Color', '#007bff');
plot(time, transpose(summation(1:t, 1)), 'LineWidth', 2, 'Color', 'black');
%plot(time, transpose(origamis(:, 11)), '-', 'Color', '#4c00ff');
%plot(time, transpose(origamis(:, 12)), '-', 'Color', 'black');
xlabel('time/s');
ylabel('Portion');
title('distribution');
grid off;
legend('1mers', '2mers', '3Loop','3mers','4mers', '5mers', '6Loop','6mers', '7mers', '8mers','9Loop', 'all');
hold off;
%{

number = [all(l, :)];
bar(number);
xlabel('-mers');
ylabel('intensity');
title('length distribution');


%}

%disp(transpose(all(end,:)))



%find alpha value
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
    alpha = connect/(connect+not_connect);
end


%melting curve
function f = melt(temp, H, S, C)
    x = C * exp(-(H - temp*S)/temp);
    f = 1-(-1+sqrt(4*x+1))/(x*2);
    %disp(f)
end