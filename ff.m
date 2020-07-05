h=1.5;
l=4.9;
count=0; r=l; R=[];L=[];alf=1;
while alf>0.01
    r=r+0.02;
    R=[R r];
    alf=atan(h/r);
    %Alf=[Alf alf];
    b=h;
    a=-b/r;
    L=[L b+a*l];
end

L=L/h/2*1;

figure; plot(R,L);


legend('L_T(R)')

R1=[5 6 7 8 9 10 11 12 13];
L1=[5.54 26.97 43.38 55.81 65.60 73.01 79.10 84.39 88.89];
L1=L1/264.54;

hold on;
plot (R1,L1,'or')

xlabel('R,расстояние до точки, м')
ylabel('L,отношение положения точки в кадре к высоте кадра')
grid on
legend('L_T(R)','L_Э(R)')