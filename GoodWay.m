function [X1 Y1] = GoodWay( I,X,Y )
%% функция выбора направления движения по картинке с вероятностями. 
% возвращает следующую точку для осмотра при движении по прямой?
N=length(I);
R=sqrt(N^2+(N/2)^2);
I0=I;

%% развертка 
dfi=atan(1/R); % максимальный угол для различия двух ближайших клеток на краю
r=0:0.33:R;
count=0;

% fi=pi/3;
% PFI=PFI+pi/2;
% fi=PFI-fi:dfi:PFI+fi;
fi=0:dfi:2*pi;
cofi=0.01+sin((0:0.01*dfi:2*pi)/2);%[ 0.01:0.01/round(length(fi)/2-1):1 1:(-0.01/round(length(fi)/2-1)):0.01];
R1=ones( [5*round(R) length(fi)]);
Nn=round((N-X)*cofi);
for hfi=fi
    count=count+1;
    [x1 y1]=pol2cart(hfi,r);
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
        R1(5*round(R)-length(c)+1:5*round(R),count)=I(c)';
        %R1(1:length(c))=sort(I(c)');
    end
end

X1=0;Y1=0;
%RR=R1(1:2:2*round(R),1:2:length(fi));
%figure(40); surfl(R1); shading interp

%% выбор оптимального пути. по максимальному столбику R1 и ширине этого столбика
% чем дальше от нуля соседние столбики, тем уже может быть наш.
R2=R1;
R2(R2>0)=1;
count=1;
h=[];
for i=1:length(fi)
    a=flip(R2(:,i));
    a(1)=0;
    a=find(a>0);
    if a
        r2(i)=a(1);
    else
        r2(i)=1;
    end
%     if i>1 && mod(i,20)==0
%         r2(i-19:i)=min(r2(i-19:i));
%     end
    if i>1 && abs( r2(i)-r2(i-1) )>=50
        r2(count:i-1)=min(r2(count:i-1));
        
        h=[h; i-count r2(i-1)+Nn(count+round((i-count)/2))];
        count=i;
    end
    if (i-count)>=100
        r2(count:i-1)=min(r2(count:i-1));
        
        h=[h; i-count r2(i-1)+Nn(count+round((i-count)/2))];
        count=i;
    end
        
end
r2(count:i)=min(r2(count:i));
h=[h; i-count r2(i-1)+Nn(count+round((i-count)/2))];
%figure(50); plot(r2);

if count>1
    if h(1,1)==1
                h(2,1)=h(2,1)+1;
                h(1,:)=[];
    end
    if find(h(:,1)==1)
        k=find(h(:,1)==1);
        for i=length((k)):-1:1
            h(k(i)-1,1)=h(k(i)-1,1)+h(k(i),1); %против выбросов
        end
        h(k,:)=[];
        if h(1,1)==2
                h(2,1)=h(2,1)+2;
                h(1,:)=[];
        end
        if find(h(:,1)==2)
            k=find(h(:,1)==2);
            for i=length((k)):-1:1
                h(k(i)-1,1)=h(k(i)-1,1)+h(k(i),1); %против выбросов
            end
            h(k,:)=[];
            %h(find(h(:,1)==2)-1,1)=h(find(h(:,1)==2)-1,1)+2; %против выбросов
            %h(find(h(:,1)==2),:)=[];
        end
    end
    h=[h h(:,2).*h(:,1) h*0 h*0 h*0];
    h(1,4)=fi(round(h(1,1)/2));
    for i=2:length(h(:,1))
        h(i,4)=fi(sum(h(1:i-1,1))+round(h(i,1)/2));
        h(i,5)=cofi(sum(h(1:i-1,1))+round(h(i,1)/2));
        h(i,6)=sum(h(1:i-1,1))+round(h(i,1)/2);
    end
    %figure(50); hold off; plot(r2); hold on; plot(h(:,6),h(:,2),'o');
    %figure(51); hold off; plot(h(:,6),h(:,3),'o');
    h(:,3)=h(:,3).*h(:,5);
    h( h(:,2)<5,3)=1;
    %figure(51); hold on; plot(h(:,6),h(:,3),'*');
    
    n=find(h(:,3)==max(h(:,3)));
    
    h=sortrows(h,3,'descend');
    
    n=1:size(h,1);%(fix(size(h,1)/2)+1);
    phi=h(n,4);
    Rphi=h(n,2)*2/3;
    [x1 y1]=pol2cart(phi,Rphi);
    X1=X+round(x1);Y1=Y+round(y1);
else
    X1=X-round(3/4*min(r2));Y1=Y;
    %X1=N+1;Y1=N+1;
end
% % [X1 Y1 RR] = GoodWay( I,X,Y );
% % R=RR(1:2:1118,1:2:1757);
% % figure(4); surfl(R); shading interp
