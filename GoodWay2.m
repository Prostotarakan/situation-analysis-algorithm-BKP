function [X1 Y1] = GoodWay2( I,X,Y,PFI )
%% функция выбора направления движения по картинке с вероятностями. 
% возвращает следующую точку для осмотра при движении по прямой?
N=length(I);
R=sqrt(max(X,abs(N-X))^2+max(Y,abs(N-Y))^2); % максимально удаленная точка
I0=I;
gran=0.1;
%% развертка 
count=0;
dfi=atan(1/R); % максимальный угол для различия двух ближайших клеток на краю
r=0:0.33:R;
% fi=pi/3;
% PFI=PFI+pi/2;
% fi=PFI-fi:dfi:PFI+fi;
fi=0:dfi:2*pi; % на весь круг

%[RR FIFI]=meshgrid(r,fi);
count=1;i=0;h=[];
for hfi=fi
    i=i+1;
[x1 y1]=pol2cart(hfi,r);
x1=round(x1)+X;y1=round(y1)+Y;
Rr=r;
Fi=hfi;
x1=reshape(x1,[],1); % переразбивка координат и вероятностей в вектор
y1=reshape(y1,[],1);
Rr=reshape(Rr,[],1);
Rr( y1<=0 ) =[];
x1( y1<=0 ) =[];
y1( y1<=0 ) =[];
Rr( x1<=0 ) =[];
y1( x1<=0 ) =[];
x1( x1<=0 ) =[];
Rr( y1>N ) =[];
x1( y1>N ) =[];
y1( y1>N ) =[];
Rr( x1>N ) =[];
y1( x1>N ) =[];
x1( x1>N ) =[];

XY=[x1 y1]; % выбивание уникальных номеров
[C, c]=unique(XY,'rows');
if c
        Rr=Rr(c);
        x1=x1(c);
        y1=y1(c);
        j=C(:,1)+( C(:,2)-1 )*N;
        C=I(j);
        li=find(C>gran);
        if li
            r2(i)=Rr(li(1));
        else
            r2(i)=length(C);
        end
else
        r2(i)=0;
end
if i>1 && abs( r2(i)-r2(i-1) )>=50
    r2(count:i-1)=min(r2(count:i-1));   
    h=[h; i-count r2(i-1)];
    count=i;
end    
          
end
r2(count:i)=min(r2(count:i));        

if count>1
    if h(1,1)==1
                h(2,1)=h(2,1)+1;
                h(1,:)=[];
        end
    if find(h(:,1)==1)
        h(find(h(:,1)==1)-1,1)=h(find(h(:,1)==1)-1,1)+1; %против выбросов
        h(find(h(:,1)==1),:)=[];
        if h(1,1)==2
                h(2,1)=h(2,1)+2;
                h(1,:)=[];
        end
        if find(h(:,1)==2)
            h(find(h(:,1)==2)-1,1)=h(find(h(:,1)==2)-1,1)+2; %против выбросов
            h(find(h(:,1)==2),:)=[];
        end
    end
    h=[h h(:,2).*h(:,1)];
    n=find(h(:,3)==max(h(:,3)));
    count=sum(h(1:n-1,1))+round(h(n,1)/2);
    phi=fi(count);
    Rphi=h(n,2)*2/3;
    [x1 y1]=pol2cart(phi,Rphi);
    X1=X+round(x1);Y1=Y+round(y1);
else
    X1=X-round(3/4*min(r2));Y1=Y;
    %X1=N+1;Y1=N+1;
end






end

