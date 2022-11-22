%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; %
%close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%initial parameters
L=2^10; % main domain size      
N=2^6; % number of points for the aperture radius
z0=0; %interface position (focus at z=0)  
z=linspace(-2,2,41); %axial domain
m0=1/2; %coordinate shift 
xmax=4; %half the size of the cropped output in wavelengths 
%(final size of images is 2xmax)
n1=1.518;  %refractive index for incident field (immersion oil)
n2=1.518; % refractive index for reflected field  (immersion oil) 
n3=1; % refractive index for transmitted field (after interface)
e1=n2^2; 
e2=n3^2;
m1=1; %permeability
m2=1; %permeability
f=1800; % effective focal length of the lens in um (10^-6 m)
NA=1.3; %effective numerical aperture  
lambda=1.064;  %wavelength in um 
l=1.064;  %wavelength   um 
R=f*NA/n2;  %aperture radius in um;   
dx_inf=R/N;  %spatial dimension at aperture um/px or um/point 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Spatial coordinate system before lens (x_inf, y_inf)
d_inf=linspace(-L/2+m0,L/2-m0,L)*dx_inf; %spatial axis shifted by m0.
[x_inf,y_inf]=meshgrid(d_inf,-d_inf); %mesh, y inverted in Matlab images
[theta_inf,rho_inf] = cart2pol(x_inf,y_inf); %auxiliary polar coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% angular 
k0=2*pi/lambda;
k1=n2*k0;
k2=n3*k0;
dk=k0*NA/N;    %dk=(1/2)*(k1/(f))*dx0;
kx=x_inf*k1/f;
ky=y_inf*k1/f;
dxf=2*pi/(L*dk); %equivalent to dxf=N*lambda/(L*NA); 
%conversion factor at output plane in um/px
kM1=k1*ones(L,L);       %%%array with magnitude of k1
kz1=sqrt(kM1.^2-kx.^2-ky.^2); 
kM2=k2*ones(L,L);  
kz2=sqrt(kM2.^2-kx.^2-ky.^2); 

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
phi_incx=double(zeros(L,L));  %%%phase phi_x0 fase init;    phi_0
E_incy=E_incx;
phi_incy = phi_incx;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%round soft aperture (spatial coordinates)
Theta=0.5*(1+tanh((1.5/1)*(N-sqrt((x_inf/dx_inf).^2+(y_inf/dx_inf).^2)))); 
%angular space
%kmax=R*k1/f;
%Theta=0.5*(1+tanh((1.5/dk)*(kmax-sqrt((kx).^2+(ky).^2))));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%polarization 
cx=ones(L,L); cy=zeros(L,L);  % lin horizontal
%cx=zeros(L,L); cy=ones(L,L);  % lin vertical
%cx=1*ones(L,L); cy=1*ones(L,L);   % lin diag 
%cx=ones(L,L); cy=1i*ones(L,L); %%right circ
%cx=1i*ones(L,L); cy=ones(L,L); %%left circ
%cx=(kx)./sqrt(kx.^2+ky.^2);  cy=(ky)./sqrt(kx.^2+ky.^2); %%radial  cx=cos(phi) cy=sin(phi)
%cx=-ky./sqrt(kx.^2+ky.^2);  cy= kx./sqrt(kx.^2+ky.^2);%%azimuthal  cx=cos(phi) cy=sin(phi)
%phik=atan2(ky,kx); s12=8; cx=cos((s12/2)*phik); cy=sin((s12/2)*phik);%flower
%phik=atan2(ky,kx); s12=-8; cx=cos((s12/2)*phik); cy=sin((s12/2)*phik);% spider

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Fresnel coeficients
rs=(m2*kz1-m1*kz2)./(m2*kz1+m1*kz2); 
rp=(e2*kz1-e1*kz2)./(e2*kz1+e1*kz2);
ts=(2*m2*kz1)./(m2*kz1+m1*kz2); 
tp=sqrt((m2*e1)/(m1*e2))*(2*e2*kz1)./(e2*kz1+e1*kz2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Factors for E_inf
CF1=-(1./kz1).*sqrt(kz1./kM1).*(1./(kx.^2+ky.^2))*1i*f*exp(-1i*f*2*pi*n2/lambda); %Factors for E_inf_r and f
CF2=-(1./kz2).*sqrt(kz1./kM1).*(1./(kx.^2+ky.^2))*1i*f*exp(-1i*f*2*pi*n2/lambda).*(kz2./kz1); %Factors for E_inf _t

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%Focused propagating before interface
%for E_inf
% x component
MTx1=(ky.^2 +kx.^2.*kz1./kM1).*CF1; %lens transmision factor, incident x
MTx2=(-ky.*kx +(ky.*kx.*kz1)./kM1).*CF1; %lens transmission factor incident y
%y component
MTy1=(-kx.*ky+(kx.*ky.*kz1)./kM1).*CF1;%lens transmision factor, incident x  
MTy2=(kx.^2 +ky.^2.*kz1./kM1).*CF1;%lens transmision factor, incident y
%z component
MTz1=(-(kx.^2+ky.^2).*kx./kM1).*CF1;%lens transmision factor, incident x
MTz2=(-(kx.^2+ky.^2).*ky./kM1).*CF1; %lens transmision factor, incident y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Reflected E_inf_r
%x component
MTx1r=(rs.*ky.^2 -rp.*kx.^2.*kz1./kM1).*CF1; % incident x
MTx2r=(-rs.*ky.*kx  -rp.*(ky.*kx.*kz1)./kM1).*CF1;% incident y
%y component
MTy1r=(-rs.*kx.*ky-rp.*(kx.*ky.*kz1)./kM1).*CF1;% incident x
MTy2r=(rs.*kx.^2 -rp.*ky.^2.*kz1./kM1).*CF1;% incident y
%z component
MTz1r=(-rp.*(kx.^2+ky.^2).*kx./kM1).*CF1;%incident x
MTz2r=(-rp.*(kx.^2+ky.^2).*ky./kM1).*CF1; % incident y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%Transmited E_inf_t
%x component
MTx1t=(ts.*ky.^2+tp.*kx.^2.*kz2./kM2).*CF2; % x input
MTx2t=(-ts.*ky.*kx+tp.*(ky.*kx.*kz2)./kM2).*CF2;% y input
%y component
MTy1t=(-ts.*kx.*ky+tp.*(kx.*ky.*kz2).*(1./kM2)).*CF2;% x input
MTy2t=(ts.*kx.^2+tp.*ky.^2.*kz2./kM2).*CF2;% y input
%z component
MTz1t=(-tp.*(kx.^2+ky.^2).*kx./kM2).*CF2;%x input
MTz2t=(-tp.*(kx.^2+ky.^2).*ky./kM2).*CF2; % y input

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%start 

for j=1:length(z); 
Hzf=exp(1i*kz1.*z(j)).*Theta;  %propagating
Hzr=exp(-1i*kz1.*z(j)+1i*2*kz1*z0).*Theta; %reflected
Hzt=exp(1i*kz2.*z(j)+1i*(kz1-kz2)*z0).*Theta; %transmitted



if z(j)<=z0;  % propagating and reflected up to the boundary at z<z0
%%%%%%%%%% focal
% x component E_inf
Fieldxxf=(cx.*MTx1.*E_incx.*exp(1i*phi_incx)+cy.*MTx2.*E_incy.*exp(1i*phi_incy)).*Hzf;
%integration x component
Ex_inff=ifftshift(ifft2(fftshift((Fieldxxf)))).*exp(1i*PhaseShift);
Exf=Ex_inff((L/2+1-xmax2):(L/2+1+xmax2),(L/2+1-xmax2):(L/2+1+xmax2)); %cropped

%y component E_inf
Fieldxyf=(cx.*MTy1.*E_incx.*exp(1i*phi_incx)+cy.*MTy2.*E_incy.*exp(1i*phi_incy)).*Hzf;
%integration y component
Ey_inff=ifftshift(ifft2((fftshift(Fieldxyf)))).*exp(1i*PhaseShift);
Eyf=Ey_inff((L/2+1-xmax2):(L/2+1+xmax2),(L/2+1-xmax2):(L/2+1+xmax2));

%z component E_inf
Fieldxzf=(cx.*MTz1.*E_incx.*exp(1i*phi_incx)+cy.*MTz2.*E_incy.*exp(1i*phi_incy)).*Hzf;
%integration z component
Ez0f=ifftshift(ifft2((fftshift(Fieldxzf)))).*exp(1i*PhaseShift);
Ezf=Ez0f((L/2+1-xmax2):(L/2+1+xmax2),(L/2+1-xmax2):(L/2+1+xmax2));
%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% reflected
%x component E_inf_r
Fieldxxr=(cx.*MTx1r.*E_incx.*exp(1i*phi_incx)+cy.*MTx2r.*E_incy.*exp(1i*phi_incy)).*Hzr;
%integration
Ex_infr=ifftshift(ifft2(fftshift((Fieldxxr)))).*exp(1i*PhaseShift);
Exr=Ex_infr((L/2+1-xmax2):(L/2+1+xmax2),(L/2+1-xmax2):(L/2+1+xmax2)); 

%y component E_inf_r
Fieldxyr=(cx.*MTy1r.*E_incx.*exp(1i*phi_incx)+cy.*MTy2r.*E_incy.*exp(1i*phi_incy)).*Hzr;
%integration
Ey_infr=ifftshift(ifft2((fftshift(Fieldxyr)))).*exp(1i*PhaseShift);
Eyr=Ey_infr((L/2+1-xmax2):(L/2+1+xmax2),(L/2+1-xmax2):(L/2+1+xmax2));

%z component E_inf_r
Fieldxzr=(cx.*MTz1r.*E_incx.*exp(1i*phi_incx)+cy.*MTz2r.*E_incy.*exp(1i*phi_incy)).*Hzr;
%integration
Ez0r=ifftshift(ifft2((fftshift(Fieldxzr)))).*exp(1i*PhaseShift);
Ezr=Ez0r((L/2+1-xmax2):(L/2+1+xmax2),(L/2+1-xmax2):(L/2+1+xmax2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ex(:,:,j)=Exf+Exr;   % total field at the material n_2 incident+reflected
Ey(:,:,j)=Eyf+Eyr;
Ez(:,:,j)=Ezf+Ezr;

%%%%%%%%%%%
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%transmitted
else
% x component E_inf_t    
Fieldxxt=(cx.*MTx1t.*E_incx.*exp(1i*phi_incx)+cy.*MTx2t.*E_incy.*exp(1i*phi_incy)).*Hzt;
%integration
Ex_inft=ifftshift(ifft2(fftshift((Fieldxxt)))).*exp(1i*PhaseShift);
Ext=Ex_inft((L/2+1-xmax2):(L/2+1+xmax2),(L/2+1-xmax2):(L/2+1+xmax2)); %cropped 

%y component E_inf_t
Fieldxyt=(cx.*MTy1t.*E_incx.*exp(1i*phi_incx)+cy.*MTy2t.*E_incy.*exp(1i*phi_incy)).*Hzt;
%integration
Ey_inft=ifftshift(ifft2((fftshift(Fieldxyt)))).*exp(1i*PhaseShift);
Eyt=Ey_inft((L/2+1-xmax2):(L/2+1+xmax2),(L/2+1-xmax2):(L/2+1+xmax2)); %cropped

%z component E_inf_t
Fieldxzt=(cx.*MTz1t.*E_incx.*exp(1i*phi_incx)+cy.*MTz2t.*E_incy.*exp(1i*phi_incy)).*Hzt;
%integration
Ez0t=ifftshift(ifft2((fftshift(Fieldxzt)))).*exp(1i*PhaseShift); 
Ezt=Ez0t((L/2+1-xmax2):(L/2+1+xmax2),(L/2+1-xmax2):(L/2+1+xmax2)); %cropped

Ex(:,:,j)=Ext;  %field beyond the boundary, transmitted
Ey(:,:,j)=Eyt;
Ez(:,:,j)=Ezt;

end

end




%%%plot yz cross section through the center at xmax2+1
figure(1)
EEx=permute(Ex(xmax2+1,:,:),[3 2 1]);
EEy=permute(Ey(xmax2+1,:,:),[3 2 1]);
EEz=permute(Ez(xmax2+1,:,:),[3 2 1]);
I=abs(EEx).^2+abs(EEy).^2+abs(EEz).^2;
Imax=max(max(I));
In=I/Imax;
imagesc((In).^(1/4)); colormap(jet); colorbar