m = 6;
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

k = 0.01;
alpha = 0.1;
beta = 0.1;
dt = 0.01;
total_t =500;
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

for i = 1:n
    
    %K=(1/a+1/b)^1/2
    md = -k*6*(sqrt(2)*m*m+sqrt(3/2)*m*d+sqrt(4/3)*m*tri+sqrt(5/4)*m*te+sqrt(6/5)*m*p+sqrt(7/6)*m*h+0.5*sqrt(8/7)*m*hep);
    dd = k*3*sqrt(2)*m*m-k*6*(sqrt(3/2)*d*m+2*d*d+sqrt(5/6)*d*tri+sqrt(3/4)*d*te+sqrt(7/10)*d*p+sqrt(2/3)*d*h+sqrt(9/14)*d*hep);
    trid = k*6*sqrt(3/2)*m*d-k*6*(sqrt(4/3)*tri*m+sqrt(5/6)*tri*d+sqrt(2/3)*tri*tri+sqrt(7/12)*tri*te+sqrt(8/15)*tri*p+0.5*sqrt(1/2)*tri*h+sqrt(10/21)*tri*hep)-alpha*tri;
    triLd = alpha*tri;
    ted = k*6*(d*d+sqrt(4/3)*tri*m)-k*6*(0.5*sqrt(5/4)*te*m+sqrt(3/4)*te*d+sqrt(7/12)*te*tri+sqrt(1/2)*te*te+sqrt(9/20)*te*p+sqrt(5/12)*te*h+0.5*sqrt(11/28)*te*hep);
    pd = k*6*(sqrt(5/6)*d*tri+0.5*sqrt(5/4)*h*m)-k*6*(sqrt(6/5)*p*m+sqrt(7/10)*p*d+sqrt(8/15)*p*tri+sqrt(9/20)*p*te+2*sqrt(2/5)*p*p+sqrt(11/30)*p*h+sqrt(12/35)*p*hep);
    hd = k*6*(sqrt(6/5)*m*p+sqrt(3/4)*d*h+0.5*sqrt(2/3)*tri*tri)-k*6*(sqrt(7/6)*h*m+sqrt(2/3)*h*d+sqrt(1/2)*h*tri+sqrt(5/12)*h*te+sqrt(11/30)*h*p+sqrt(1/3)*h*h+sqrt(13/42)*h*hep)-beta*h;
    hLd = beta*h;
    hepd = k*6*(sqrt(7/6)*m*h+sqrt(7/10)*d*p+sqrt(7/12)*tri*te)-k*6*(0.5*sqrt(8/7)*hep*m+sqrt(9/14)*hep*d+sqrt(10/21)*hep*tri+0.5*sqrt(11/28)*hep*te+sqrt(12/35)*hep*p+sqrt(13/42)*hep*h+sqrt(14/49)*hep*hep);
%{   
    %K=(1/a+b)^1/2
    md = -k*6*(sqrt(1/2)*m*m+sqrt(1/3)*m*d+sqrt(1/4)*m*tri+sqrt(1/5)*m*te+sqrt(1/6)*m*p+sqrt(1/7)*m*h+0.5*sqrt(1/8)*m*hep);
    dd = k*3*sqrt(1/2)*m*m-k*6*(sqrt(1/3)*d*m+2*sqrt(1/4)*d*d+sqrt(1/5)*d*tri+sqrt(1/6)*d*te+sqrt(1/7)*d*p+sqrt(1/8)*d*h+sqrt(1/9)*d*hep);
    trid = k*6*sqrt(1/3)*m*d-k*6*(sqrt(1/4)*tri*m+sqrt(1/5)*tri*d+sqrt(1/6)*tri*tri+sqrt(1/7)*tri*te+sqrt(1/8)*tri*p+0.5*sqrt(1/9)*tri*h+sqrt(1/10)*tri*hep)-alpha*tri;
    triLd = alpha*tri;
    ted = k*6*(sqrt(1/4)*d*d+sqrt(1/4)*tri*m)-k*6*(0.5*sqrt(1/5)*te*m+sqrt(1/6)*te*d+sqrt(1/7)*te*tri+sqrt(1/8)*te*te+sqrt(1/9)*te*p+sqrt(1/10)*te*h+0.5*sqrt(1/11)*te*hep);
    pd = k*6*(sqrt(1/5)*d*tri+0.5*sqrt(1/5)*h*m)-k*6*(sqrt(1/6)*p*m+sqrt(1/7)*p*d+sqrt(1/8)*p*tri+sqrt(1/9)*p*te+2*sqrt(1/10)*p*p+sqrt(1/11)*p*h+sqrt(1/12)*p*hep);
    hd = k*6*(sqrt(1/6)*m*p+sqrt(1/6)*d*h+0.5*sqrt(1/6)*tri*tri)-k*6*(sqrt(1/7)*h*m+sqrt(1/8)*h*d+sqrt(1/9)*h*tri+sqrt(1/10)*h*te+sqrt(1/11)*h*p+sqrt(1/12)*h*h+sqrt(1/13)*h*hep)-beta*h;
    hLd = beta*h;
    hepd = k*6*(sqrt(1/7)*m*h+sqrt(1/7)*d*p+sqrt(1/7)*tri*te)-k*6*(0.5*sqrt(1/8)*hep*m+sqrt(1/9)*hep*d+sqrt(1/10)*hep*tri+0.5*sqrt(1/11)*hep*te+sqrt(1/12)*hep*p+sqrt(1/13)*hep*h+sqrt(1/14)*hep*hep);
 %}    

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
plot(time, listtrimers, '-', 'Color', '#fff200'); % Line for y3 will be green
%plot(time, listtri, '-', 'Color', '#1aff00');
plot(time, listtriL, '-', 'Color', '#00ffe5');
plot(time, listte, '-', 'Color', '#0084ff');
plot(time, listp, '-', 'Color', '#0400ff');
plot(time, listh, '-', 'Color', '#8800ff');
plot(time, listhL, '-', 'Color', '#ff00e6');
plot(time, listhexmers, '-', 'Color', '#1EFF00');
plot(time, listhep, '-', 'Color', '#4c00ff');
plot(time, origami, '-', 'Color', 'black');
xlabel('time');
ylabel('numbers');
title('time distribution (intensity)');
grid on;
legend('1', '2','3', '3L','4', '5', '6', '6L', '6all', '7','ori');
hold off;



