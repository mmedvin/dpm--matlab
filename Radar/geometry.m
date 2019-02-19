%%
r1=1.4;
r0=0.4;
r2=0.9;

r=linspace(r0,r1,2^11); 
r3=r(end-20);

plot(r2*cos(th)+1/3,r2*sin(th)+1/3); 
hold on; 
plot(r0*cos(th),r0*sin(th));
plot(r1*cos(th),r1*sin(th)); 
plot(r3*cos(th),r3*sin(th)); 

hold off;
axis([-r1 r1 -r1 r1])

dr=r(2)-r(1)

d=r1-(sqrt(2*(1/3)^2)+0.9); d/dr
d=(sqrt(2*(1/3)^2)-0.9)+r0; d/dr
d=r1-r3;d/dr