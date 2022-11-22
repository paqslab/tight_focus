%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; %close all; clc
%run time sim 10s with L=2^12, and N=2^7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%initial parameters
L=2^12; % main domain size      
N=2^7; % number of points for the aperture radius
z=0; %observation plane (0 at the focal plane)  
m0=1/2; %coordinate shift 
xmax=0.9; %half the size of the cropped output in wavelengths 
%(final size of images is 2xmax)
n1=1.0;  % refractive index before the lens (air)
n2=1.518; % refractive index after the lens (immersion oil) 
f=1800; % effective focal length of the lens in um (10^-6 m)
NA=1.4; %effective numerical aperture  
lambda=1.064;  %wavelength in um 
R=f*NA/n2;  %aperture radius in um;   
dx_inf=R/N;  %spatial resolution at aperture um/px or um/point 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Spatial coordinate system before lens (x_inf, y_inf)
d_inf=linspace(-L/2+m0,L/2-m0,L)*dx_inf; %spatial axis shifted by m0.
[x_inf,y_inf]=meshgrid(d_inf,-d_inf); %mesh, y inverted in Matlab images
[theta_inf,rho_inf] = cart2pol(x_inf,y_inf); %auxiliary polar coordinates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% angular 
k0=2*pi/lambda;
k1=n2*k0;
dk=k0*NA/N;    %dk=(1/2)*(k1/(f))*dx0;
kx=x_inf*k1/f;
ky=y_inf*k1/f;
dxf=2*pi/(L*dk); %equivalent to dxf=N*lambda/(L*NA); 
%conversion factor at output plane in um/px
kM1=k1*ones(L,L);       %%%array with magnitude of k1
kz1=sqrt(kM1.^2-kx.^2-ky.^2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmax2=round(xmax*lambda/dxf); %xmax in px at output plane

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%correction for m0 shift
x=linspace(-L/2,L/2,L); %%%auxiliary coordinate for correction
y=x; [X,Y]=meshgrid(x,y);
PhaseShift=m0*2*pi*X/L+m0*2*pi*Y/L; %correction phase
%Center of Fourier transform displaced (vert and hor) by m0.
%shift of \delta kx=\delta ky= m0*(2*pi/L). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%input field   initial field transverse
w0=f*NA/n2;  % The beam waist is equal to the aperture radius. 
E_incx=exp(-(x_inf.^2+y_inf.^2)/(w0^2));  %normalized amplitude E_x0   A_0
phi_incx=double(zeros(L,L));  %%%constant
E_incy=E_incx;
phi_incy = phi_incx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%round soft aperture (spatial coordinates)
Theta=0.5*(1+tanh((1.5/1)*(N-sqrt((x_inf/dx_inf).^2+(y_inf/dx_inf).^2)))); 
%angular space
%kmax=R*k1/f;
%Theta=0.5*(1+tanh((1.5/dk)*(kmax-sqrt((kx).^2+(ky).^2))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%polarization state
%cx=ones(L,L); cy=zeros(L,L);  % lin horizontal
cx=zeros(L,L); cy=ones(L,L);  % lin vertical
%cx=1*ones(L,L); cy=1*ones(L,L);   % lin diag 
%cx=ones(L,L); cy=1i*ones(L,L); %%right circ
%cx=1i*ones(L,L); cy=ones(L,L); %%left circ
%cx=kx./sqrt(kx.^2+ky.^2);  cy=ky./sqrt(kx.^2+ky.^2); %%radial  
%cx=-ky./sqrt(kx.^2+ky.^2);  cy= kx./sqrt(kx.^2+ky.^2);%%azimuthal  
%phi_aux=mod(atan2(ky,kx),2*pi); s=8; cx=cos((s/2)*phi_aux); cy=sin((s/2)*phi_aux); %flower



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%E_inf  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input x
%x component
E_infxx=sqrt(n1/n2)*sqrt(kz1./kM1).*Theta.*E_incx.*exp(1i*phi_incx).*...
    (ky.^2+kx.^2.*kz1./kM1)./(kx.^2+ky.^2);
%y component
E_infxy=sqrt(n1/n2)*sqrt(kz1./kM1).*Theta.*E_incx.*exp(1i*phi_incx).*...
    (-kx.*ky+kx.*ky.*kz1./kM1)./(kx.^2+ky.^2);
%z component
E_infxz=sqrt(n1/n2)*sqrt(kz1./kM1).*Theta.*E_incx.*exp(1i*phi_incx).*...
    (-(kx.^2+ky.^2).*kx./kM1)./(kx.^2+ky.^2);

%y input
%x component
E_infyx=sqrt(n1/n2)*sqrt(kz1./kM1).*Theta.*E_incy.*exp(1i*phi_incy).*...
    (-ky.*kx+kx.*ky.*kz1./kM1)./(kx.^2+ky.^2);
%y component
E_infyy=sqrt(n1/n2)*sqrt(kz1./kM1).*Theta.*E_incy.*exp(1i*phi_incy).*...
    (kx.^2+ky.^2.*kz1./kM1)./(kx.^2+ky.^2);
%z component
E_infyz=sqrt(n1/n2)*sqrt(kz1./kM1).*Theta.*E_incy.*exp(1i*phi_incy).*...
    (-(kx.^2+ky.^2).*ky./kM1)./(kx.^2+ky.^2);

%factors to assemble E_inf
CF2=-1i*f*exp(-1i*f*2*pi*n2/lambda)*exp(1i*kz1.*z)./(2*pi*kz1);  

%%E_inf
Fieldx=CF2.*(cx.*E_infxx +cy.*E_infyx); 
Fieldy=CF2.*(cx.*E_infxy +cy.*E_infyy);
Fieldz=CF2.*(cx.*E_infxz +cy.*E_infyz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Integration
Ex0=ifftshift(ifft2(fftshift((Fieldx)))).*exp(1i*PhaseShift);
Ex=Ex0((L/2+1-xmax2):(L/2+1+xmax2),(L/2+1-xmax2):(L/2+1+xmax2)); %cropped 
%x=y=0 at the center

Ey0=ifftshift(ifft2((fftshift(Fieldy)))).*exp(1i*PhaseShift);
Ey=Ey0((L/2+1-xmax2):(L/2+1+xmax2),(L/2+1-xmax2):(L/2+1+xmax2)); %cropped

Ez0=ifftshift(ifft2((fftshift(Fieldz)))).*exp(1i*PhaseShift);
Ez=Ez0((L/2+1-xmax2):(L/2+1+xmax2),(L/2+1-xmax2):(L/2+1+xmax2)); %cropped

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%intensity ratios
c0y=max(max(abs(Ey).^2))/max(max(abs(Ex).^2));
c0z=max(max(abs(Ez).^2))/max(max(abs(Ex).^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plots

a='jet';
b='gray';

figure(1) %amplitude Ex
Exn=abs(Ex)./max(max(abs(Ex)));
imagesc(abs(Exn)); colormap(b); axis image
set(gca,'xtick',[]); set(gca,'ytick',[])

figure(2) %phase Ex
imagesc(angle(Ex)+pi); axis image 
colormap(a); caxis([0 2*pi]); 
set(gca,'xtick',[]); set(gca,'ytick',[])

figure(3) %amplitude Ey
Eyn=abs(Ey)./max(max(abs(Ey)));
imagesc(Eyn); colormap(b); axis image
set(gca,'xtick',[]); set(gca,'ytick',[])

figure(4) %phase Ey
imagesc(angle(Ey)+pi); axis image
colormap(a); caxis([0 2*pi])
set(gca,'xtick',[]); set(gca,'ytick',[])

figure(5) %amplitude Ez
Ezn=abs(Ez)./max(max((abs(Ez))));
imagesc(Ezn); axis image
colormap(b); 
set(gca,'xtick',[]); set(gca,'ytick',[])

figure(6) %phase Ez
imagesc(angle(Ez)+pi); axis image
%imagesc(mod(atan2(imag(Ez),real(Ez)),2*pi)); axis image
colormap(a); caxis([0 2*pi])
set(gca,'xtick',[]); set(gca,'ytick',[])
