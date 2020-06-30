%%
r1=1.4;
r0=0.4;
r2=0.9;

gridparam = 8;
shift = [0,1/3];%[1/3,1/3];


th=linspace(0,2*pi,100);

r=linspace(r0,r1,2^gridparam); 
r3=r(end-40);

plot(r2*cos(th)+shift(1),r2*sin(th)+shift(2)); 
hold on; 
plot(r0*cos(th),r0*sin(th));
plot(r1*cos(th),r1*sin(th)); 
plot(r3*cos(th),r3*sin(th)); 

hold off;
axis([-r1 r1 -r1 r1])

dr=r(2)-r(1)

R = sqrt(shift(1)^2 + shift(2)^2);
d=r1-(R+r2); d/dr
d=(R-r2)+r0; d/dr
d=r1-r3;d/dr


