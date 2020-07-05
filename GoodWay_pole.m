function [X1 Y1] = GoodWay_pole( I,X,Y, C)

if C==0
    %вверх
    [X1 Y1] = GoodWay_Roy( I,X,Y );
else
    %в цель
    
%xc=sum(Xc)/length(Xc);
%yc=sum(Yc)/length(Yc); %среднее, где лежит объект.

[x y]=meshgrid(1:length(I));
x=x-X;
y=y-Y;
[ Fit Rt]=cart2pol(y,x); %радиус до каждой точки


Ia=Rt.*I; %красивая штука, нужен минимум по прямой. Используем рой.

Ia(Rt<10)=0;
%figure; surfl(Ia)

%% рой
%% по поиску минимума суммы движения

N=length(Ia);
R=sqrt(N^2+(N/2)^2);
I0=Ia;
% Rii=[];
% for fi=80*pi/100:pi/100:130*pi/100
%     Ri=0:0.33:(N/2-1);
%     
%     [x y]=pol2cart(fi,Ri);
%     x1=round(x)+X;y1=round(y)+Y;
%     cx=x1+( y1-1 )*N;
%     Rii=[Rii; Ia(cx)];
% end
% figure(118); surfl(Rii)
% figure(120); plot(Rii(35,:))
% figure(120); hold on; plot(cumsum(Rii(35,:)));    
n=15;%число разведчиков
m=15;%число рабочих

fi=rand([1 n])*2*pi;
dfi=pi/12;
dr=20;

M=sum(sum(Ia>0))+1;
flag=1;count=0; minR_=0; co=0;
while flag==1 && count<100
count=count+1;
    %разведчики
for i=1:length(fi)
[x1 y1]=pol2cart(fi(i),0:0.33:R);
x1=round(x1)+X;y1=round(y1)+Y;
x1( y1<=0 ) =[];
y1( y1<=0 ) =[];
y1( x1<=0 ) =[];
x1( x1<=0 ) =[];
x1( y1>N ) =[];
y1( y1>N ) =[];
y1( x1>N ) =[];
x1( x1>N ) =[];

R1(i)=M;
if x1
    
    cx=x1+( y1-1 )*N;
    Rr=cumsum(Ia(cx));
    k=find(Rr==min(Rr));
    k=k(end);
    
    if i==1
        min(Rr)
    end
    x=x1(k)-X;y=y1(k)-Y;
    
    R1(i)=sqrt(x^2+y^2);
    K1(i)=sign(min(Rr));
    
  
end
    %рабочие
    fi2=rand([1 m])*dfi/2-dfi+fi(i);
    
    R2=[];
    for j=1:length(fi2)
        [x1 y1]=pol2cart(fi2(j),0:0.33:R);
        x1=round(x1)+X;y1=round(y1)+Y;
        x1( y1<=0 ) =[];
        y1( y1<=0 ) =[];
        y1( x1<=0 ) =[];
        x1( x1<=0 ) =[];
        x1( y1>N ) =[];
        y1( y1>N ) =[];
        y1( x1>N ) =[];
        x1( x1>N ) =[];
        
        R2(j)=M;
        if x1
            
            cx=x1+( y1-1 )*N;
            
            Rr=cumsum(Ia(cx));
            
            k=find(Rr==min(Rr));
            k=k(end);
               
            x=x1(k)-X;y=y1(k)-Y;

            R2(j)=sqrt(x^2+y^2);
            K2(j)=sign(min(Rr));
            
       end     
            
        if R2(j)<M
            if R2(j)>R1(i) && K2(j)<=K1(i)
                
                fi(i)=fi2(j);
                R1(i)=R2(j);
                K1(i)=K2(j)
            end
        end
        
    end
       
end

if abs(minR_-min(R1.*K1))<2
   co=co+1;
   if co==3
        flag=0;  
   end
elseif minR_>min(R1.*K1)
    minR_=min(R1.*K1); 
end

end

An=[fi' R1' K1'];
An=sortrows(An,[3 2],{'ascend','descend'});
pfi=An(:,1);
Rpfi=An(:,2);
pfi(Rpfi<1)=[];
Rpfi(Rpfi<1)=[];
[x1 y1]=pol2cart(pfi,Rpfi*7/9);
X1=X+round(x1);Y1=Y+round(y1);

end
end

