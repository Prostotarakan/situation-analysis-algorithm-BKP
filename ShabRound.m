function S=ShabRound(r,Mu) %r=0..1
%% функци€ дл€ шаблона круга с градиентом 
% чем меньше r тем меньше круг, r - соотношение удаленности и максимальной удаленности
n=30; % 15 n - коэффициент дл€ размера получившегос€ круга
N=fix(n*r)*2+1;

S=zeros(N);
R=fix(N/2);
[x y]=meshgrid(-fix(n*r):fix(n*r));
[fi,r] = cart2pol(x,y);
d=r/R; % функци€ €ркости круга - 1-r/R
d=exp(r/R); %по гауссу эта и следующа€ строчки
%d=flip(d);
S(r<=R)=3-d(r<=R);
if N<=1
    S=ones(3);
end

if nargin==2
    S=S*Mu;
end
end

