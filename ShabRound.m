function S=ShabRound(r,Mu) %r=0..1
%% ������� ��� ������� ����� � ���������� 
% ��� ������ r ��� ������ ����, r - ����������� ����������� � ������������ �����������
n=30; % 15 n - ����������� ��� ������� ������������� �����
N=fix(n*r)*2+1;

S=zeros(N);
R=fix(N/2);
[x y]=meshgrid(-fix(n*r):fix(n*r));
[fi,r] = cart2pol(x,y);
d=r/R; % ������� ������� ����� - 1-r/R
d=exp(r/R); %�� ������ ��� � ��������� �������
%d=flip(d);
S(r<=R)=3-d(r<=R);
if N<=1
    S=ones(3);
end

if nargin==2
    S=S*Mu;
end
end

