mol_m = 5;  %6 * 10^-9;
m = mol_m * 6*10^-10;%1nM = 10^-6 * 6*10^23/10^27 = 6*10^-10 nm^-3
H = -83.5*4184.0 / 8.3;%79.8 * 10^3
S = -218.9*4.184 / 8.3;%-222.6
i=283;
j=350;
%{
%vertical se single
H = -83.5*4184.0 / 8.3;%79.8 * 10^3
S = -218.9*4.184 / 8.3;%-222.6

%Horizental se single
H = -107.5*2*4184.0 / 8.3;
S = -313.5*2*4.184 / 8.3;%looks more likely 

%Horizental se 18
H4 = -58.3*4184.0 / 8.3;%79.8 * 10^3
S4 = -154.3*4.184 / 8.3;%-222.6
%}

frac = zeros(4,j-i);
Tem = zeros(1,j-i);


for temp = i:j
    x = m * exp(-(H - temp*S)/temp);
    f = 1-(-1+sqrt(4*x+1))/(x*2);
    Tem(1,temp-i+1)=temp-273;
    frac(1,temp-i+1)=f;
%{
    x2 = m * exp(-(H2 - temp*S2)/temp);
    x3 = m * exp(-(H3 - temp*S3)/temp);
    x4 = m * exp(-(H4 - temp*S4)/temp);
    f2 = 1-(-1+sqrt(4*x2+1))/(x2*2);
    f3 = 1-(-1+sqrt(4*x3+1))/(x3*2);
    f4 = 1-(-1+sqrt(4*x4+1))/(x4*2);
    frac(1,temp-i+1)=f;
    frac(2,temp-i+1)=f2;
    frac(3,temp-i+1)=f3;
    frac(4,temp-i+1)=f4;
%}
end

plot(Tem,frac(1,:),'LineWidth', 2, 'Color', '#ff1900');
hold on;
plot(Tem,frac(3,:),'LineWidth', 2, 'Color', '#0062ff');
plot(Tem,frac(4,:),'LineWidth', 2, 'Color', '#ff6a00');
plot(Tem,frac(2,:),'LineWidth', 2, 'Color', '#7700ff');


%legend('6 Horizental SE on Origami','Single Horizental SE', '18 Vertical SE', 'Singel Vertical SE');
xlabel('Temperature(Celsius)');
ylabel('fraction of monomers');
title('Approximate Melting Curve');
hold off;
