% % % % 
FLAG=0;
if FLAG==1
    clear all
    clc
    load Is
    load Xcentr
    load Ycentr
    A=Is;
    N=length(Is);
    X=fix(length(Is)/2);
    Y=X;
    C=0; %есть цель
    %Xc=1;Yc=1;%%%!!!
else
  clear all
  clc
    N=200;
    n=15;
     
   A=zeros(N);
   A(randi(N^2,[n 1]))=1;
    [x0,y0]=find(A==1);
    x0(x0>(N-20))=N-20;
    A=zeros(N);
    A(x0+( y0-1 )*N)=1;
    
    X=N;Y=round(N/2);
    
    C=0; %наличие цели
end

stupid=1308; % в 308 они не умеют в имшоу!!!


Z=1;%0.75;%коэффициент забывания
    
figure(1)
if stupid==308
    image(A*255)
else
    imshow(A)
end


%% камера в 200;100

X0=X;Y0=Y;
Xend=X;Yend=Y;
[x0,y0]=find(A==1);
i=1:10:length(x0);
x0=x0(i);
y0=y0(i);

%% зациклим
flag=1;
Pfi=pi/2; % направление зрения

N1=100;
I2=ones(N);
Ind=zeros(N+2*N1);

[XX,YY]=meshgrid(1:N);
con=0;cof=0;co=1;
while flag==1
x=x0;y=y0;
%% область зрения:

fi=pi/6;%pi/6; %60  pi/3; % если не видно 30 градусов

[hfi,r]=cart2pol(XX-X,YY-Y);
hfi2=hfi;
hfi2(hfi>0)=hfi(hfi>0)-pi/2;
hfi2(hfi<=0)=hfi(hfi<=0)+3*pi/2;
hfi=hfi2;
i=find(hfi<(Pfi-fi) );
i=unique([ i ; find(hfi>(Pfi+fi) ) ]);
%% карта глубины с точки зрения робота:
I=A;

I=zeros(N);
I(i)=1; 
I01=I';
I=I';
r=[];ri=[];
%I(x+( y-1 )*N)=1;
[x y]=find( A.*(1-I)==1 );
% i=find(A==1);
% I(i)=1;

%figure(18); imshow(I)
In_s=zeros(size(I));
r=[];

tic
for i=1:length(x)
    r(i)=sqrt((X-x(i))^2+(Y-y(i))^2)/sqrt(N^2+(N/2)^2); %0..1
    
    if C==1
        xa=find(x(i)==Xc);
        ya=find(y(i)==Yc);
        if length(find(xa==ya'))
            % если это цель
            Mu=-1;
            r(i)=2*r(i);
            Shab=ShabRound(r(i),Mu);
            [n m]=size(Shab);
            n=fix(n/2);
            m=fix(m/2);
            ri(i)=n;
            if ri(i)==0
                ri(i)=1;
            end
            In_=zeros(N+2*N1);
            In_( x(i)+N1-n:x(i)+N1+n , y(i)+N1-m:y(i)+N1+m )=Shab;

            In_=In_(1+N1:N+N1,1+N1:N+N1);
            
            In_s=min(In_s,In_);
            
        else
            % если просто препятствия
            Mu=1;
            Shab=ShabRound(r(i));
            [n m]=size(Shab);
            n=fix(n/2);
            m=fix(m/2);
            ri(i)=n;
            if ri(i)==0
                ri(i)=1;
            end
            In=zeros(N+2*N1);
            In( x(i)+N1-n:x(i)+N1+n , y(i)+N1-m:y(i)+N1+m )=Shab;

            In=In(1+N1:N+N1,1+N1:N+N1);

            I=max(I,In);%*(1-r(i)/2); %%!!!
        end
        
    else
        % если просто препятствия как обычно
        Shab=ShabRound(r(i));
        [n m]=size(Shab);
        n=fix(n/2);
        m=fix(m/2);
        ri(i)=n;
        if ri(i)==0
            ri(i)=1;
        end
        In=zeros(N+2*N1);
        In( x(i)+N1-n:x(i)+N1+n , y(i)+N1-m:y(i)+N1+m )=Shab;

        In=In(1+N1:N+N1,1+N1:N+N1);

        I=max(I,In);%*(1-r(i)/2); 
    end

end
I=I/max(max(I));
I=max(I,I01);
In_s=In_s/max(max(abs(In_s)));
if C==1
    I=I+2*In_s;
end
%figure(45);imshow(I)
t_shab=toc

tic
if C==1
    I1=Trass(I,X,Y,x,y,ri,Xc,Yc);
else
    I1=Trass(I,X,Y,x,y,ri);
end
t_trass=toc
%figure(2); imshow(I1)
I=I/max(max(I));

I1(I1>1)=1;
I2=(I1.*(Z*I2+(1-Z)/2)); %забывание
figure(3)
hold off;
if stupid==308
    image((A+I2)*255)
else
    imshow(A+I2)
end
hold on
plot(Yend,Xend,'ro')

% figure(31)
% [hfi,r]=cart2pol(XX-X,YY-Y);
% hfi2=hfi;
% hfi2(hfi>0)=hfi(hfi>0)-pi/2;
% hfi2(hfi<=0)=hfi(hfi<=0)+3*pi/2;
% i=find(hfi2<(Pfi-fi) );
% i=unique([ i ; find(hfi2>(Pfi+fi) ) ]);
% I3=I2';I3(i)=0;
% surfl(I3)

%
cof=0;
while cof<3
    tic
[X1 Y1] = GoodWay_Roy( I2,X,Y );

%[X1 Y1] = GoodWay( I2,X,Y );

%[X1 Y1] = GoodWay_pole( I2,X,Y, C )
t_wayroy=toc
% Shab=ShabRound(0.1);
% [n m]=size(Shab);
% n=fix(n/2);
% m=fix(m/2);
% 
% Ind( X+N1-n:X+N1+n , Y+N1-m:Y+N1+m )=Shab;   
% Ind( X1+N1-n:X1+N1+n , Y1+N1-m:Y1+N1+m )=Shab;   
% Indi=Ind(1+N1:N+N1,1+N1:N+N1);  
% figure(4)
% if stupid==308
%     image((A+Indi)*255)
% else
%     imshow(A+Indi)
% end
con=con+1;


Y1(X1>=(N))=[];
X1(X1>=(N))=[];
X1(Y1>=(N))=[];
Y1(Y1>=(N))=[];

Y1(X1<1)=[];
X1(X1<1)=[];
X1(Y1<1)=[];
Y1(Y1<1)=[];
nu=5;
Y1( abs(X-X1)<nu )=[];
X1( abs(X-X1)<nu )=[];
X1( abs(Y-Y1)<nu )=[];
Y1( abs(Y-Y1)<nu )=[];
% nu=0;
% X2=[200 125 80 32];
% Y2=[100 85 85 102 100];
% co=co+1; cof=3;
% X1=X2(co);
% Y1=Y2(co);

if length(find( X1<nu) ) 
    flag=0;
    Y1=Y1(find( X1<nu));
    Y1=Y1(1);
    X1=X1(find( X1<nu));
    X1=X1(1);
    Xend=[Xend;X1];Yend=[Yend;Y1];
    cof=3;
elseif C==1
    flag=0;
    Xend=[Xend;X1(1)];Yend=[Yend;Y1(1)];
    cof=3;

else
    
if length(X1)>=1
    %cof=0;
    flag=1;
%     figure;
%     plot([X X1],[Y Y1],X,Y,'*')
    conti=1;
    nu=10;
    xend=round(Xend/nu)+1;
    yend=round(Yend/nu)+1;
    xend=xend+( yend-1 )*N;
    cofi=0;
    while conti>0 && cofi<=length(X1)
        [a b]=find( repmat(xend,[1 length(xend)])==xend');
        if length(a)>length(xend)
            flag=0;
            conti=-2;
        else
            x=round(X1(conti)/nu)+1;
            y=round(Y1(conti)/nu)+1;
            x=x+( y-1 )*N;
            if length(find(xend==x))>=1
                %conti=conti+1;
                cofi=cofi+1;
            else
                X1=X1(conti);
                Y1=Y1(conti);
                conti=0;
            end
        end
    end
    if conti==0    

    Pfi=atan( abs((Y-Y1) / (X-X1)) );
    if X1<=X && Y1>=Y
        Pfi=atan( abs((X-X1) / (Y-Y1)) );
    elseif X1>=X && Y1>=Y
        Pfi=Pfi-pi/2;
    elseif X1>=X && Y1<=Y
        Pfi=atan( abs((X-X1) / (Y-Y1)) )+pi;
    else
        Pfi=Pfi+pi/2;
        
    end
    X=X1;
    Y=Y1;
    Xend=[Xend;X];Yend=[Yend;Y];
    cof=3;
%     else
%         flag=0;
    end
    if con>=30
        flag=0;
    end
    %Pfi=pi/2;
else
    cof=cof+1;
    if cof==3
        flag=0;
    end
end

end


end
% figure(4)
% hold off;
% if stupid==308
%     image((A+I2)*255)
% else
%     imshow(A+I2)
% end
% hold on
% plot(Yend,Xend,'r')


% figure(5)
% hold off
% if stupid==308
%     image((A+I1)*255)%+Indi
% else
%     imshow(A+I1)%+Indi
% end
% hold on
% plot(Yend,Xend,'r')

if C==1
    flag=0;
end

end

figure(1)
hold on
plot(Yend,Xend,'g')

con
%% шоб записать себе картинку
%imwrite(I,'testmap.jpg')

%% используемые функции
% function S=ShabRound(r) %r=0..1
% %% функция для шаблона круга с градиентом
% % чем меньше r тем меньше круг, r - соотношение удаленности и максимальной удаленности
% end

% function S=Trass(I,X,Y,x,y,ri)
% %% функция для зарисовки в невидимую область всего за объектом (функция определения зон с вероятностью объекта в них)
% % поступает карта с кругами, положение камеры, положение центра кругов и их радиус
% end

%%

% % % %Монте-Карло для n и m
% % % N=200;
% % % for i=1:1000
% % %     n(i)=rand(1)*N;
% % %     m(i)=rand(1)*N;
% % %     [con(i),rmin(i)]=test_way(n(i),m(i));
% % % end
% % % An=[n' m' con' rmin' con'./rmin'];
% % % An=sortrows(An,5);
% % % n=An(1,1)
% % % m=An(1,2)
% % % 
% % % %Q-learning для n и m
% % % 
% % % n=rand(1)*N;
% % % m=rand(1)*N;
% % % flag=1;
% % % while flag==1
% % % [con,rmin]=test_way(n,m);
% % % r1=con/rmin;
% % % if abs(r1-r)<2
% % %     flag=0;
% % % else
% % %     r=r1;
% % %     n=n+rand(1)*r;
% % %     m=m+rand(1)*r;
% % % end
% % % end
% % % n
% % % m


    


    
    
