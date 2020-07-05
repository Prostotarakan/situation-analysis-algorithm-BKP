% % % % ��� ����� � �������� ����, ��� ����� ����������� � ������. - ����
% % % % 
% % % % ����������� � ����� ������!!!
% % % % 
% % % % ���������� �����, ���� ����� ������ �� �������� ���� ������.
% % % % ���� �� ���.




    clear all
    clc
    N=200;
    n=15;
    
    load X
    load Y
    
    Y=Y+100;
    
    x0=round(X);
    y0=round(Y);
    N=200;
    y0(x0<=0)=[];
    x0(x0<=0)=[];
    x0(y0<=0)=[];
    y0(y0<=0)=[];
    y0(x0>N)=[];
    x0(x0>N)=[];
    x0(y0>N)=[];
    y0(y0>N)=[];
    N=200;
    A=zeros(N);
    A(x0+( y0-1 )*N)=1;
    
    load Xk
    load Yk
    N=200;
    Xk=round(Xk);
    Yk=round(Yk);
    X=Xk(1);Y=Yk(1);
    
    C=0; %������� ����


stupid=1308; % � 308 ��� �� ����� � �����!!!


Z=1;%0.75;%����������� ���������
    
figure(1)
if stupid==308
    image(A*255)
else
    imshow(A)
end

hold on; plot(X,Y,'or')

%% ������ � 200;100

X0=X;Y0=Y;
Xend=X;Yend=Y;
[x0,y0]=find(A==1);

%% ��������
flag=1;
Pfi=pi/2; % ����������� ������

N1=100;
I2=ones(N);
Ind=zeros(N+2*N1);

[XX,YY]=meshgrid(1:N);
con=2;
while flag==1
x=x0;y=y0;
%% ������� ������:

fi=pi/6; %60  pi/3; % ���� �� ����� 30 ��������

[hfi,r]=cart2pol(XX-X,YY-Y);
hfi2=hfi;
hfi2(hfi>0)=hfi(hfi>0)-pi/2;
hfi2(hfi<=0)=hfi(hfi<=0)+3*pi/2;
hfi=hfi2;
i=find(hfi<(Pfi-fi) );
i=unique([ i ; find(hfi>(Pfi+fi) ) ]);
%% ����� ������� � ����� ������ ������:

I=zeros(N);
I(i)=1; 
I=I';
r=[];ri=[];
%I(x+( y-1 )*N)=1;
Ia=I+A;
[x y]=find( A.*(1-I)==1 );
I=Ia;
%figure(18); imshow(I)
In_s=zeros(size(I));

for i=1:length(x)
    r(i)=sqrt((X-x(i))^2+(Y-y(i))^2)/sqrt(N^2+(N/2)^2); %0..1
    
    if C==1
        xa=find(x(i)==Xc);
        ya=find(y(i)==Yc);
        if length(find(xa==ya'))
            % ���� ��� ����
            Mu=-1;
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
            % ���� ������ �����������
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
        % ���� ������ ����������� ��� ������
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

I=I+In_s;
%figure(45);imshow(I)
if C==1
    I1=Trass(I,X,Y,x,y,ri,Xc,Yc);
else
    I1=Trass(I,X,Y,x,y,ri);
end
%figure(2); imshow(I1)
%I=I/max(max(I));

I1(I1>1)=1;
I2=(I1.*(Z*I2+(1-Z)/2)); %���������
figure(3)
hold off;
if stupid==308
    image((A+I2)*255)
else
    imshow(A+I2)
end
hold on
plot(Yend,Xend,'r')
plot(Yend(end),Xend(end),'or')   
X1=Xk(con);
Y1=Yk(con);
con=con+1;

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

if con>length(Xk)
    flag=0;

end



end

figure(1)
hold on
plot(Yend,Xend,'o-g')

con
%% ��� �������� ���� ��������
%imwrite(I,'testmap.jpg')

%% ������������ �������
% function S=ShabRound(r) %r=0..1
% %% ������� ��� ������� ����� � ����������
% % ��� ������ r ��� ������ ����, r - ����������� ����������� � ������������ �����������
% end

% function S=Trass(I,X,Y,x,y,ri)
% %% ������� ��� ��������� � ��������� ������� ����� �� �������� (������� ����������� ��� � ������������ ������� � ���)
% % ��������� ����� � �������, ��������� ������, ��������� ������ ������ � �� ������
% end

%%



    


    
    
