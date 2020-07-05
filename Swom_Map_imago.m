clear all
clc

tic
%угол обзора камеры примерно 50х30 градусов 
%видео 3, 29 секунда (доехал до крайней точки чтоб развернуться)

I=imread('D:\BMSTU\I.jpg'); %сам кадр
I_gl=imread('D:\BMSTU\I_gl.jpg'); %глубина
I_r=imread('D:\BMSTU\I_r.jpg'); %кадр с сегментацией
%кадр с рамкой.
I_r(1:10,:,:)=[];
I_r(:,1:10,:)=[];
I_r(484:493,:,:)=[];
I_r(:,647:656,:)=[];

%image(I_r);

a=[0.0452   0    0.6597    0.1337];

[X Y ~]=size(I);

[X_gl Y_gl ~]=size(I_gl);
if X_gl~=X
    I=imresize(I,[X_gl Y_gl]);
    [X Y ~]=size(I);
    I_r=imresize(I_r,[X_gl Y_gl]);
    I_r=rgb2ind(I_r,3);
end

A=a.*[X Y X Y]; %вершины области с объектом-целью для I

cv=[];co=[];
R=[];
N=0;
for i=1:Y
    Ix=I_r(:,i);
    r=length(find(Ix==1));
    R(i)=round(0.5*0.1/(0.1-r/(X_gl)*2*0.1)*100);
    
    if length(find(Ix==2))
        cv(i)=1;
    end
    if length(find(Ix==0))
        co(i)=1;
    end
    
    
    if length(find(Ix==0))==0
        N=max(N,R(i));
    end
        
end
R(R>N)=N;
R(abs(R-N)<N/20)=N;
co(R>=N)=0;
cv(R>=N)=1;

%figure(2); plot(R)
H=round(N)+1; % подвязать это к карте глубины?!!!
R=H-abs(R-N); %или минус макс R
%figure(3); plot(R)

Is=zeros(2*H);
X0=fix(length(Is)/2);
Y0=X0;

Fi=pi/3;
dfi=Fi/(Y-1);
fi=[0:dfi:Fi]-Fi/2;
Xc=[];Yc=[];
for i=1:Y
if co(i)==1   
    if R(i)~=N
    [xa ya]=pol2cart(fi(i),R(i));
    x1=round(xa)+X0;y1=round(ya)+Y0;
    x1=length(Is)-x1;
    Is(x1,y1)=1;

    if i<Y
        if R(i+1)~=R(i)
            if R(i)>R(i+1)
                r1=R(i+1):0.6:R(i);
            else
                r1=R(i):0.6:R(i+1);
            end
            [xa ya]=pol2cart(fi(i),r1);
            x1=round(xa)+X0;
            y1=round(ya)+Y0;
            x1=length(Is)-x1;
            c=x1+( y1-1 )*length(Is);
            Is(c)=1;
        end
    end
    end
end

if i<A(4) %если искомый объект
    Xc=[Xc x1];
    Yc=[Yc y1];
end
  
% figure(5);
% imshow(Is)
if cv(i)==1
    [xa ya]=pol2cart(fi(i),R(i));
    y1=round(ya)+Y0;
    x1=round(H-N); %разница между H и размером картины
    Is(x1,y1)=1;
end

  

end

toc

figure(5);
imshow(Is)
hold on; plot(Yc,Xc,'g')

% [xa ya]=pol2cart(fi(23),r(23));
% x1=round(xa)+X0;y1=round(ya)+Y0;
% x1=length(Is)-x1;
% hold on; plot(y1,x1,'or')




% условный результат сегментации. подвязать к Диме!
% x=[0    55    55   165   165   181   181   235   235   418   418   500];
% y=[180   180   167   167   199   199   174   174   194   194   176   176];

% figure(1);
% imshow(I)
% hold on; plot(x,y,'r')

%Im=zeros(size(I));
% L=[];
% x0=0;y0=1;
% 
% for i=1:2:12
%     x0=x(i);
%     if x0==0
%         x0=1;
%     end
%     x1=x(i+1);
%     y1=y(i);
%     
%     Im(y0:y1,x0:x1)=1;
%     
%     L=[L; sum(sum(I_gl((y1-10):y1 ,x0:x1 ,1))) sum(sum(I_gl((y1-10):y1 ,x0:x1 ,2))) sum(sum(I_gl((y1-10):y1 ,x0:x1 ,3)))];
% end

% figure(2);
% image(Im)
% 
% figure(3);
% imshow(I_gl)
% hold on; plot(x,y,'r')

% L
% sum(L')

% H=2; % считаю, что до стены 2 метра - подвязать это к карте глубины!!!
% R=100*H-(y-y(3));

% 
% Is=zeros(3*H*100);
% X0=fix(length(Is)/2);
% Y0=X0;
% 
% r=zeros([1 Y]);
% for i=1:Y
%     x0=x(i);
%     if x0==0
%         x0=1;
%     end
%     x1=x(i+1);
%     
%     r(x0:x1)=R(i);
%     
% end
% 
% Fi=pi/3;
% dfi=Fi/(Y-1);
% fi=[0:dfi:Fi]-Fi/2;
% for i=1:Y
%     
% [xa ya]=pol2cart(fi(i),r(i));
% x1=round(xa)+X0;y1=round(ya)+Y0;
% x1=length(Is)-x1;
% Is(x1,y1)=1;
% 
% if i<Y
%     if r(i+1)~=r(i)
%         if r(i)>r(i+1)
%             r1=r(i+1):0.6:r(i);
%         else
%             r1=r(i):0.6:r(i+1);
%         end
%         [xa ya]=pol2cart(fi(i),r1);
%         x1=round(xa)+X0;
%         y1=round(ya)+Y0;
%         x1=length(Is)-x1;
%         c=x1+( y1-1 )*length(Is);
%         Is(c)=1;
%     end
% end
% 
% end
% 
% figure(5);
% imshow(Is)
% 
% [xa ya]=pol2cart(fi(23),r(23));
% x1=round(xa)+X0;y1=round(ya)+Y0;
% x1=length(Is)-x1;
% hold on; plot(y1,x1,'or')
% 
