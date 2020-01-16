function [xx,ww]=retgauss_precision(a,b,n,m,xd,d);
% function [xx,ww]=retgauss(a,b,n,m,xd,d);
%
% integration methode de gauss degre n de a jusque à    b (éventuellement séparée en m morceaux)
%  xx est le tableau des points ou calculer la fonction a integrer (vecteur ligne)
%  ww est le tableau des poids associes (vecteur ligne)
% si n>200 methode des rectangles a n*abs(m) points
% si n<0 methode des trapezes a abs(n) points
% si m<0 on rajoute les bornes a et b avec des poids nuls
%
%  si xd est precise: xd points de discontinuite,(eventuellement definis modulo une periode d )
%   on repartit au moins m*n points entre les discontinuites
%
%  exemples
%%%%%%%%%%%%%%%%
%  [xx,ww]=retgauss(0,1,3)   xx=[0.1127,0.5, 0.8873]    ww=[0.27778,0.44444,0.27778]
%  [xx,ww]=retgauss(0,1,3,-1);xx=[0 , 0.1127, 0.5, 0.8873,1]    ww=[0, 0.27778 ,0.44444, 0.27778,0]  
%  [xx,ww]=retgauss(0,1,-3)  ou  [xx,ww]=retgauss(0,1,-3,-2); xx=[0 , 0.5 , 1]    ww=[0.25, 0.5 , 0.25]  

%n=331;retplot;[x,w]=retgauss(-1,1,n);for ii=1:1.3*n;a=2/(1-(2*(ii-1))^2);retplot(ii,log10(abs   ((sum(w.*cos(acos(x)*2*(ii-1)))-a)  /a)));end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1;%exist('retgauss_precision.mat','file')~=2; % calcul des coefficients et ecriture sur fichier
if exist('retgauss_precision_prv.mat','file')==2; % calcul des coefficients et ecriture sur fichier
load retgauss_precision_prv;
nn0=nn+1;
else;
nmax=501;
x=cell(1,nmax);w=cell(1,nmax);
x{1}=[];x{2}=[0];
nn0=3;
end;
figure;
for nn=nn0:nmax;mm=length(x{nn-1});tic;
ipair=(mod(nn,2)==0);
t=([1:mm-1]/mm);t=.5*(1-.97*t);
x{nn}=x{nn-1};x{nn}(1:mm-1)=x{nn}(1:mm-1).*(1-t)+x{nn}(2:mm).*t;% approximation de depart
y=[0,x{nn},1];
xx=retprecis(zeros(ipair+length(x{nn}),1));
for ii=1:length(x{nn});
% x{nn}(ii)=abs(real(retcadilhac(@legendre,[],x{nn}(ii),nn)));
% [zp,itermax,erz,erfonc,test]=retcadilhac(@legendre,struct('tol',1.e-30,'tolf',1.e-30'),retprecis(x{nn}(ii)),nn);

[zp,itermax,erz,erfonc,test]=retcadilhac(@legendre,struct('tol',1.e-8,'tolf',1.e-8,'nitermin',4,'precis',i,'tol_precis',20),x{nn}(ii)+[0;y(ii)-x{nn}(ii);y(ii+2)-x{nn}(ii)]/2,nn),
disp(rettexte(test,nn,ii,itermax));
zp=abs(real(zp));xx(ii+ipair)=zp;%  xx contient les x en precision etendue
x{nn}(ii)=double(zp);
end; 
if mod(nn,2)==0;x{nn}=[0,x{nn}];end;

tx=toc;tic;
% calcul des poids
mm=length(x{nn});
%aaa=legendre(xx,[1:2:2*mm-1]);if mod(nn,2)==0;aaa(:,1)=aaa(:,1)/2;end;bbb=zeros(mm,1);bbb(1)=1;bbb=retprecis(bbb);
%aaa=legendre(double(xx),[1:2:2*mm-1]);if mod(nn,2)==0;aaa(:,1)=aaa(:,1)/2;end;bbb=zeros(mm,1);bbb(1)=1;%bbb=retprecis(bbb);


aaa=zeros(0,mm);bbb=[];for ii=1:mm;aaa=[aaa;cos(acos(x{nn})*2*(ii-1))];bbb=[bbb;1/(1-(2*(ii-1))^2)];end;if mod(nn,2)==0;aaa(:,1)=aaa(:,1)/2;end;
%for ii=1:mm-1;aaa=[aaa;legendre(x{nn},2*ii+1)];bbb=[bbb;0];end;if mod(nn,2)==0;aaa(:,1)=aaa(:,1)/2;end;









w{nn}=(aaa\bbb)';

%w{nn}=double(w{nn});
subplot(2,1,1); plot(x{nn},'.');title(rettexte(nn));subplot(2,1,2); plot(w{nn},'.');drawnow;
tw=toc;disp(rettexte(nn,tx,tw));
save retgauss_precision_prv nmax x w nn;
x=x(1:nn);nmax=nn;for nnn=2:2:nn;x{nnn}=x{nnn}(2:end);end;save retgauss_precision nmax x w;
load retgauss_precision_prv;





end;
%for nn=2:2:nmax;x{nn}=x{nn}(2:end);end;
stop;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4;m=1;end; 
if m<0;bornes=1;else;bornes=-1;end;m=abs(m);% pour rajouter les bornes

if nargin<5; % pas de discontinuitees

if n<0;  % methode des trapezes d'ordre abs(n);
n=abs(n);
if n==1;xx=(a+b)/2;ww=(b-a);
else;xx=linspace(a,b,n);ww=((b-a)/(n-1))*ones(1,n);ww([1,n])=.5*ww([1,n]);end;
else;    % methode de gauss ou des rectangles
    
    
if n*m==0;xx=[];ww=[];return;end;

load retgauss_precision nmax;
if (n<nmax)&(n>1)% methode de gauss
load retgauss_precision;x=x{n+1};w=w{n+1};
if mod(n,2)==0 x=[-x(n/2:-1:1),x];w=[w(n/2:-1:1),w];else x=[-x((n-1)/2:-1:1),0,x];w=[w((n+1)/2:-1:2),w];end;

xx=[];ww=[];
for ii=1:m;% on separe en m morceaux l'intervalle [a  b]
aa=a+(ii-1)*(b-a)/m;bb=a+ii*(b-a)/m;    
xx=[xx,(bb+aa)/2+(bb-aa)/2.*x];ww=[ww,w.*((bb-aa)/2)];
end;

else; % methode des rectangles d'ordre n*m
n=n*m;xx=a+((b-a)/n)*([1:n]-.5);ww=(b-a)/n.*ones(1,n);   
end;
end;  % trapezes ou (gauss ou rectangles)

else; % discontinuitees en xd
if b<a;
if nargin>5;[xx,ww]=retgauss(b,a,n,-m*bornes,xd,d);else;[xx,ww]=retgauss(b,a,n,-m*bornes,xd);end;
xx=fliplr(xx);ww=-fliplr(ww);return;end;% donc maintenant b>=a ...    
% 'mise en forme' de xd
xd=retelimine(sort(xd));
if nargin>5;% discontinuitees definies modulo d
xd=mod(xd,d);
n1=ceil((a-max(xd))/d);
n2=floor((b-min(xd))/d);
xdd=[];
for ii=n1:n2;xdd=[xdd,xd+ii*d];end;
xd=xdd;    
end;    
xd=xd((xd>a)&(xd<b));  
xd=retelimine([a,sort(xd),b]);
% calcul des xx et ww    
xdd=xd(2:end)-xd(1:end-1);mm=ceil(m*xdd/sum(xdd)); % repartition proportionnelle aux longueurs   
xx=[];ww=[];
for ii=1:length(xdd);
[xg,wg]=retgauss(xd(ii),xd(ii+1),n,mm(ii));
xx=[xx,xg];ww=[ww,wg];
end;
end;

if bornes>0;
if a~=xx(1);xx=[a,xx];ww=[0,ww];end;
if b~=xx(end);xx=[xx,b];ww=[ww,0];end;
end;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù
function zz=legendre1(z,n); 
if n==1;zz=ones(size(z));return;end;
if n==2;zz=z;return;end;
zz1=ones(size(z));zz2=z;
for ii=1:n-2;zz=((2*ii+1)*z.*zz2-ii*zz1)/(ii+1);zz1=zz2;zz2=zz;end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù
function zz=legendre(z,n);z=z(:).';
zz=zeros(max(n),length(z));zz1=ones(1,length(z));zz2=z;
switch class(z);case 'retprecis';zz=retprecis(zz);zz1=retprecis(zz1);end;    
zz(1,:)=ones(1,length(z));zz(2,:)=z;
for ii=1:max(n)-2;zz(ii+2,:)=((2*ii+1)*z.*zz2-ii*zz1)/(ii+1);zz1=zz2;zz2=zz(ii+2,:);end;
zz=zz(n,:);
