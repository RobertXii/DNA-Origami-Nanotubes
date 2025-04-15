mol_m = 6;
m = mol_m * 6.02 * 10^23;
d = 0;
tri = 0;
triL = 0;
te = 0;
p = 0;
h = 0;
hL = 0;
hep = 0;

md = 0;
dd = 0;
trid = 0;
triLd = 0;
ted = 0;
pd = 0;
trimers = 0;
hd = 0;
hLd = 0;
hexmers = 0;
hepd = 0;

Temperature = 40;
M = (Temperature + 273) * 1.53 * 10^-15;
N = 3.48 * 10^16;
alpha = 0.02;
beta = 0.01;
dt = 1;
total_t =1000000;
n = total_t/dt;

listm = zeros(1, n);
listd =  zeros(1, n);
listtri =  zeros(1, n);
listtriL =  zeros(1, n);
listte =  zeros(1, n);
listp =  zeros(1, n);
listtrimers =  zeros(1,n);
listh =  zeros(1,n);
listhL =  zeros(1,n);
listhexmers =  zeros(1,n);
listhep =  zeros(1,n);
time = zeros(1,n);
origami = zeros(1,n);

const = zeros(7,7);

for j = 1:7
    for k = 1:7
        const(j,k) = M*(j+k)^2 / (k*j+N*(k+j)*(k*j)^3/(k^3+j^3));
        if mod(k, 3) == 1 && mod(j, 3) == 1
            const(j,k) = const(j,k);
        elseif mod(k, 3) == 0 && mod(j, 3) == 0
            const(j,k) = const(j,k);
        else 
            const(j,k) = const(j,k)*2;
        end
    end
end

disp(const);

for i = 1:n
    
    md = -2*m*m*const(1,1) - m*d*const(1,2) - m*tri*const(1,3) - m*te*const(1,4) - m*p*const(1,5) - m*h*const(1,6) - m*hep*const(1,7);
    dd = m*m*const(1,1) - d*m*const(2,1) - 2*d*d*const(2,2) - d*tri*const(2,3) - d*te*const(2,4) - d*p*const(2,5) - d*h*const(2,6) - d*hep*const(2,7);
    trid = m*d*const(2,1) - m*tri*const(3,1)- tri*d*const(3,2) - 2*tri*tri*const(3,3) - tri*te*const(3,4) - tri*p*const(3,5) - tri*h*const(3,6) - tri*hep*const(3,7) - tri*alpha;
    triLd = alpha * tri;
    ted = d*d*const(2,2) + m*tri*const(3,1) - m*te*const(4,1)- te*d*const(4,2) - te*tri*const(4,3) - 2*te*te*const(4,4) - te*p*const(4,5) - te*h*const(4,6) - te*hep*const(4,7);
    pd = tri*d*const(3,2) + m*te*const(4,1) - m*p*const(5,1)- p*d*const(5,2) - p*tri*const(5,3) - p*te*const(5,4) - 2*p*p*const(5,5) - p*h*const(5,6) - p*hep*const(5,7);
    hd = tri*tri*const(3,3) + d*te*const(4,2) + m*p*const(5,1) - m*h*const(6,1)- h*d*const(6,2) - h*tri*const(6,3) - h*te*const(6,4) - h*p*const(6,5) - 2*h*h*const(6,6) - h*hep*const(6,7) - beta*h;
    hLd =  beta*h;
    hepd = tri*te*const(4,3) + d*p*const(5,2) + m*h*const(6,1) - m*hep*const(7,1)- hep*d*const(7,2) - hep*tri*const(7,3) - hep*te*const(7,4) - hep*p*const(7,5) - hep*h*const(7,6) - hep*hep*const(7,7);

    m = m + dt * md;
    d = d + dt * dd;
    tri = tri + dt * trid;
    triL = triL + dt * triLd;
    te = te + dt * ted;
    p = p + dt * pd;
    h = h + dt * hd;
    hL = hL + dt * hLd;
    hep = hep + dt * hepd;
    trimers = tri + triL;
    hexmers = h + hL;
    
    total_origami = m+2*d+3*trimers+4*te+5*p+6*hexmers+7*hep;
     %{
    listm(i) = m;
    listd(i) = d;
    listtri(i) = tri;
    listtriL(i) = triL;
    listte(i) = te;
    listp(i) = p;
    listh(i) = h;
    listhL(i) = hL;
    listhep(i) = hep;
    listtrimers(i) = trimers;
    listhexmers(i) = hexmers;
    time(i) = i*dt;
    origami(i) = total_origami;
 %} 

    listm(i) = m;
    listd(i) = d*2;
    listtri(i) = tri*3;
    listtriL(i) = triL*3;
    listte(i) = te*4;
    listp(i) = p*5;
    listh(i) = h*6;
    listhL(i) = hL*6;
    listhep(i) = hep*7;
    listtrimers(i) = trimers*3;
    listhexmers(i) = hexmers*6;
    time(i) = i*dt;
    origami(i) = total_origami;
    %count


end

plot(time, listm, '-', 'Color', '#ff1900');  % Line for y1 will be blue
hold on;  % This command allows you to add additional lines to the same plot
plot(time, listd, '-', 'Color', '#ff8c00');   % Line for y2 will be red
%plot(time, listtrimers, '-', 'Color', '#fff200'); % Line for y3 will be green
%plot(time, listtri, '-', 'Color', '#1aff00');
plot(time, listtriL, '-', 'Color', '#00ffe5');
plot(time, listte, '-', 'Color', '#0084ff');
plot(time, listp, '-', 'Color', '#0400ff');
%plot(time, listh, '-', 'Color', '#8800ff');
plot(time, listhL, '-', 'Color', '#ff00e6');
%plot(time, listhexmers, '-', 'Color', '#1EFF00');
plot(time, listhep, '-', 'Color', '#4c00ff');
plot(time, origami, '-', 'Color', 'black');
xlabel('time');
ylabel('numbers');
title('time distribution (intensity)');
grid on;
legend('1', '2', '3L','4', '5', '6L', '7','ori');
hold off;
%{
disp(listm(end));
disp(listd(end));
disp(listtriL(end));
disp(listtri(end));
disp(listte(end));
disp(listp(end));
disp(listhL(end));
disp(listh(end));
disp(listhep(end));

l = length(listm);
value = {'m','d','triL','tri','4','5','6L','6','7'};
number = [listm(l), listd(l), listtriL(l), listtri(l), listte(l), listp(l), listhL(l), listh(l), listhep(l)];

bar(number);
xticks(1:length(value));
xticklabels(value);

xlabel('-mers');
ylabel('intensity');
title('Histogram');
%}