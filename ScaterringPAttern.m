
%================ScatteringPattern======================%

%By:Abdiladi Azemi ETH Zurich

%The Matlab Code for Calculatinf the Scattering Pattern E(?)=-?_(n=1)^N I'nH0(g?) e^(/jhz)

%The program essentially calculates In

%=========================================================%


clear; format compact;


thearray=1; %The array 

%Phio and Thetao define the axis of propagation.

THETAO=pi/2.0;
LAMBDA=1.0;
EO=1.0;
R=1.125*LAMBDA;
K=2.0*pi/LAMBDA;
AA=0.05/K;
S=1.0/K;
H=K*cos(THETAO);
G=sqrt(K^2-H^2);

%Define Wire locations

if thearray ==1
    NN=15;
    PHIO=40;
    X=zeros(1,NN);
    Y=linspace(-1,1,NN);
elseif thearry ==2
    NN=30;
    PHIO=0.0;
    PHI2=linspace(-pi/2,pi/2,NN);
    X=R*cos(PHI2);
    Y=R*sin(PHI2);
end
%Calculate RHO

[Mx,My]=meshgrid(X,Y);
RHO = sqrt((Mx-Mx').^2+(My-My').^2)+diag(AA*ones(1,NN));

%Construction Matrix A

A=besselh(0,2,G*RHO);

%Construction Matrix B

ALPHA = X*sin(THETAO)*cos(PHIO)+Y*sin(THETAO)*sin(PHIO);
B = EO*exp( -i*K*ALPHA).';

%Solve for Matrix
A=inv(A);
I=A*B;

%Finally, Calculate the Scaterring Patern
PHI=linspace(0,pi,128);
ALP=cos(PHI.')*X+sin(PHI.')*Y;
E=abs(exp(i*G*ALP)*I);

%Plot RHI

if thearray ==2
    figure(1),plot(PHI*180/pi,E)
    title('Scattering Pattern')
    xlabel('\phi (Degrees)')
    ylabel('E(\phi)')
    grid on
elseif thearray ==1
    
    %ang= PHI*180/pi/90;
    %field=fftshift(E);
    
    figure(2),plot(PHI*180/pi,E)
    title('Figure-A.Azemi')
    xlabel('\phi (Degrees)')
    ylabel('E(\phi)')
    grid on
end
