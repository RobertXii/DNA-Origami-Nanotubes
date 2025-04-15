initialTemp = 326;
ramp =1;
finalTemp = 298;
eff_alpha = 0.67;
N = 4.436 * 10^17 * 9 * 10^-16 / eff_alpha; %times r^2
dt = 1;
total_t = 3600*(initialTemp-finalTemp)/ramp+3600*10;
disp(total_t/3600)
mu = 0.0006;
H = -58.3*4184.0 / 8.3;
S = -154.3*4.184 / 8.3;


kinds = 45;
dpm = ramp/60;
n = total_t/dt;
npm = 60/dt;
o = 0;
t = 1;

mol_triL = 2 * 10^-9;
triL = mol_triL * 6.02 * 10^23 * 10^3;
m_nm3 = mol_triL * 0.1;%1nM = 10^-6 * 6*10^23/10^27 = 6*10^-10 nm^-3

%disp(origamis)
origamis = zeros(n, kinds+1);
const = zeros(kinds,kinds);
origamis(1,1) = triL;

for r = 1:n+(initialTemp-finalTemp)*60/ramp
    if mod(o, npm) == 0 && initialTemp - finalTemp >= 0.1 %change the kernal every 1min
        initialTemp = initialTemp - dpm;
        M = (initialTemp) * 2 * 1.38 * 10^-23 * melt(initialTemp, H, S, m_nm3)^2 / (3 * mu);
        for j = 1:kinds
            for k = 1:kinds
                const(j,k) = 2 * M *(j+k)^2 / (k*j + N *(k+j)*(k*j)^3/(k^3+j^3));
            end
        end

        o = o+1;
        %disp(initialTemp);
        %disp(const);

    else
        dori = zeros(1,kinds);
        for i = 1:kinds 
            for j = 1:kinds %rate of loosing
                if j == i 
                    dori(1,i) = dori(1,i) - 2*origamis(t,i)*origamis(t,j)*const(i,j);
                else
                    dori(1,i) = dori(1,i) - origamis(t,i)*origamis(t,j)*const(i,j);
                end
            end
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
            dori(1,i) = dori(1,i) + add;
            origamis(t+1,i) = origamis(t,i) + dt*dori(1,i);
        end
        t = t+1;
        o = o+1;
    end
end


time = 0:dt:t-1;
plot(time, transpose(origamis(:, 1)), 'LineWidth', 2, 'Color', '#ff1900');  % Line for y1 will be blue
hold on;  % This command allows you to add additional lines to the same plot
plot(time, transpose(origamis(:, 2)*2), 'LineWidth', 2, 'Color', '#ff8c00');   % Line for y2 will be red
plot(time, transpose(origamis(:, 3)*3), 'LineWidth', 2, 'Color', '#fff200'); % Line for y3 will be green
plot(time, transpose(origamis(:, 4)*4), 'LineWidth', 2, 'Color', '#1aff00');
plot(time, transpose(origamis(:, 5)*5), 'LineWidth', 2, 'Color', '#00ffe5');
plot(time, transpose(origamis(:, 6)*6), 'LineWidth', 2, 'Color', '#0084ff');
plot(time, transpose(origamis(:, 7)*7), 'LineWidth', 2, 'Color', '#0400ff');
plot(time, transpose(origamis(:, 8)*8), 'LineWidth', 2, 'Color', '#8800ff');
plot(time, transpose(origamis(:, 9)*9), 'LineWidth', 2, 'Color', '#ff00e6');
plot(time, transpose(origamis(:, 10)*10), 'LineWidth', 2, 'Color', '#1EFF00');
plot(time, transpose(origamis(:, 11)*11), 'LineWidth', 2, 'Color', '#4c00ff');
plot(time, transpose(origamis(:, 12)*12), 'LineWidth', 2, 'Color', 'black');
xlabel('time/s');
ylabel('Count');
% title('distribution');
set(gca, 'FontSize', 20);
set(gca, 'FontName', 'Times New Roman'); % Change font for axes labels
grid off;
legend('1', '2', '3','4', '5', '6', '7','8','9','10','11','12');
hold off;
%{


number = [origamis(t, :)];

bar(number);

xlabel('-mers');
ylabel('intensity');
title('length distribution');

disp(transpose(origamis(end,:)))
%}


function f = melt(temp, H, S, C)
    x = C * exp(-(H - temp*S)/temp);
    f = 1-(-1+sqrt(4*x+1))/(x*2);
    %disp(f)
end
