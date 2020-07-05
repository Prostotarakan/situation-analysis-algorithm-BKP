function [X1 Y1] = GoodWay_Roy( I,X,Y)
%роевой алгоритм пчелинной колонии для поиска следующей точки на карте I

N=length(I);
R=sqrt(N^2+(N/2)^2);
I0=I;

n=15;%число разведчиков
m=15;%число рабочих

fi=rand([1 n])*2*pi;
dfi=pi/6;
r=0:0.33:R;
maxR_o=0;
flag=1;count=0;
while flag==1 && count<100
count=count+1;
    %разведчики
for i=1:length(fi)
[x1 y1]=pol2cart(fi(i),r);
x1=round(x1)+X;y1=round(y1)+Y;
x1( y1<=0 ) =[];
y1( y1<=0 ) =[];
y1( x1<=0 ) =[];
x1( x1<=0 ) =[];
x1( y1>N ) =[];
y1( y1>N ) =[];
y1( x1>N ) =[];
x1( x1>N ) =[];
XY=[x1' y1'];
[C c]=unique(XY,'rows');
if c
    c=C(:,1)+( C(:,2)-1 )*N;
    R1=I(c)';
    R1=flip(R1);
    R1(1)=0;
    R1=find(R1>0);
    
    if R1
        r1(i)=r(R1(1));
    else
        r1(i)=0;
    end
end
    %рабочие
    fi2=rand([1 m])*dfi-dfi/2+fi(i);
    for j=1:length(fi2)
        [x1 y1]=pol2cart(fi2(j),r);
        x1=round(x1)+X;y1=round(y1)+Y;
        x1( y1<=0 ) =[];
        y1( y1<=0 ) =[];
        y1( x1<=0 ) =[];
        x1( x1<=0 ) =[];
        x1( y1>N ) =[];
        y1( y1>N ) =[];
        y1( x1>N ) =[];
        x1( x1>N ) =[];
        XY=[x1' y1'];
        [C c]=unique(XY,'rows');
        if c
            c=C(:,1)+( C(:,2)-1 )*N;
            R1=I(c)';
            R1=flip(R1);
            R1=find(R1>0);
            if R1
                if r(R1(1))*(1+sin(fi2(j)))>r1(i)*(1+sin(fi(i)))
                    r1(i)=r(R1(1));
                    fi(i)=fi2(j);
                end
            end
        end
    end
    
    
end
%dfi=dfi*99/100;% чуть уменьшим следующую область

maxR_n=max(r1.*(1+sin(fi)));
if abs(maxR_n-maxR_o)<2
    flag=0;
    RR=max(r1);
    FI=fi(find(r1==max(r1)));
else
    maxR_o=maxR_n;
end
end

An=[fi' r1' (r1.*(1+sin(fi)))'];
An=sortrows(An,3,'descend');
pfi=An(:,1);
Rpfi=An(:,2);
pfi(Rpfi<1)=[];
Rpfi(Rpfi<1)=[];
% [Rpfi d]=sortrows(r1',1,'descend');
% phi=fi(d);

[x1 y1]=pol2cart(pfi,Rpfi*8/9);
X1=X+round(x1);Y1=Y+round(y1);


end

