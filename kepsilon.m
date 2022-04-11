clearvars;
close all;
clc;
n=400; %number of grid points
L=2; %Width of channel
vi(1:n)=1/395;%viscosity
del_y=L/n;%stencile size
%l=L/2;
%for i=1:n
%vit(i)=((l*(0.14-0.08*(1-((i*del_y)./l).^2)-0.06*(1-((i*del_y)./l).^4))).^2)*(i*del_y);%Mixing length model
vit(1:n)=0.025;
%end
ep(1:n)=90; %Guess value of ep(dissipation rate)
for i=1:n
kn(i)=((vit(i)*ep(i))./0.09).^0.5;% Guess value of TKE
end

nm=100;% dummy Value
N=1;
u=1;
%while (nm>10.^-4)

for i=1:4000
    i
%Equation1 Velocity
up=u;
u = equation1(n,L,vit,vi);
for i=1:n
     if i==1
        dudy(i)=((2*u(i))./(del_y));
        elseif i==n
        dudy(i)=(2*u(i))./(del_y);
        else
        dudy(i)=(u(i+1)-u(i-1))./(2*del_y);
        end
end   

%Equation2 Turbulent Kienetic Energy
sig_k=1;
kg=kn;
kn=equation2(n,L,u,vit,ep,vi,sig_k,kg);

if(isnan(kn))
    break;
end

%Equation3 Dissipation
sig_ep=1.3;
c1=1.44;
c2=1.92;
epg=ep;
ep= equation3(n,L,u,vit,epg,vi,sig_ep,c1,c2,kn);


rf=0.98;% Relaxation Factor
%ep=((1-rf).*ep)+rf.*epg; 
%kn=((1-rf).*kn)+rf.*kg;
%u=((1-rf)*u)+rf*up;
if(N>100)
nm=norm(u-up); % To calculate Norm
vitn=vit;
vit=0.09*(kn.^2)./ep; %Turbulent viscocity
%vit = ((1-rf)*vitn)+rf*vit;
end
N=N+1; %Iteration Count
end
for i=1:n
      if i==1
          dkn(i)=(kn(i+2)-2*kn(i+1)+kn(i))./(del_y.^2);
      elseif i==n
          dkn(i)=(kn(i)-2*kn(i-1)+kn(i-2))./(del_y.^2); 
      else
          dkn(i)=(kn(i+1)-2*kn(i)+kn(i-1))./(del_y.^2); 
      end
end
for i=1:n
vis_diff(i) = vi(i)*dkn(i)./395;
turb_diff(i) = vit(i)*dkn(i)./395;
end       
%plotting
% Read DNS data [half-channel is given (till centerline)]
load y_dns.dat
load u_dns.dat
load u2_dns.dat
load v2_dns.dat
load w2_dns.dat
load uv_dns.dat
load dns_data.dat
for i=1:1:n
    epsilon_dns = 395.*dns_data(:,2);
    production_dns = dns_data(:,3);
    press_diff_dns = dns_data(:,4);
    turb_diff_dns = dns_data(:,5);
    visc_diff_dns = dns_data(:,6);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nu=1/395;
ustar=1;
rho = 1;
kappa=0.41;
% k-epsilon model constants
c_mu=0.09;
c1_eps=1.44;
c2_eps=1.92;
sigma_k=1;
sigma_eps=1.3;

% Residual error limit
residue_limit = 10^(-4);
k_dns=0.5*(u2_dns+v2_dns+w2_dns);
eps_dns=dns_data(:,2)*ustar^4/nu; % eps is normalized by ustar^4/nu

x=(0.5*del_y):del_y:L/2;

figure(1)
plot(u_dns,y_dns,'bo');
xlabel('U'); ylabel('y/h'); title('U-velocity');
hold on
plot(u(1:n/2),x,'o--m');
legend('DNS','Model','Location','Best'); legend boxoff;

figure(2)
plot(y_dns,k_dns,'bo');
xlabel('y/h'); ylabel('k'); title('Turbulence kinetic energy');
hold on
plot(x,kn(1:n/2),'o--m');
legend('DNS','Model','Location','Best'); legend boxoff;

figure(3)
plot(y_dns,-uv_dns,'bo')
xlabel('y/h'); ylabel('-<uv>'); title('Turbulence shear stress');
hold on
 for i=2:n-1
 uv(i)=(vit(i)*dudy(i));
 end
 uv(1)=(vit(1)*dudy(1));
 uv(n)=(vit(n)*dudy(n));

plot(x,uv(1:n/2),'o--m');
legend('DNS','Model','Location','Best'); legend boxoff;

figure(4)
plot(y_dns,dns_data(:,3)/nu,'bo');
xlabel('y/h'); ylabel('P_k'); title('Production rate of k');
hold on
for i=2:n-1
pk(i)=(vit(i)*(dudy(i).^2));
end
pk(1)=(vit(1)*(dudy(1).^2));
pk(n)=(vit(n)*(dudy(n).^2));
plot(x,pk(1:n/2),'o--m');
legend('DNS','Model','Location','Best'); legend boxoff;

figure(5)
for i=1:97
    vdns(i)=0.09*(k_dns(i).^2)./eps_dns(i);
   end
plot(y_dns,vdns,'bo');
xlabel('y/h'); ylabel('P_k'); title('Turbulent Viscosity');
hold on
plot(x,vit(1:n/2),'o--m');
legend('DNS','Model','Location','Best'); legend boxoff;

figure(6)
plot(y_dns,eps_dns,'bo')
xlabel('y/h'); ylabel('\epsilon'); title('Dissipation rate of k');
hold on
plot(x,ep(1:n/2),'o--m');
legend('DNS','Model','Location','Best'); legend boxoff;

figure(7);
plot(y_dns,visc_diff_dns,'bo')
xlabel('y/h'); ylabel('Viscous diffusion'); title('Viscous diffusion');
hold on
plot(x,vis_diff(1:n/2),'o--m');
legend('DNS','Model','Location','Best');legend boxoff;

figure(8);
plot(y_dns,turb_diff_dns,'bo')
xlabel('y/h'); ylabel('Turbulent diffusion'); title('Turbulent diffusion of KE');
hold on
plot(x,turb_diff(1:n/2),'o--m');
legend('DNS','Model','Location','Best');legend boxoff;



