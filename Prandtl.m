%Analytical determination of aerodynamic coefficients for a tapered wing

%defining wing geometry
b = 0.3;                %wingspan
Cr = 0.04;              %root chord length
Ct = 0.02;              %wing tip length
AR = (2*b)/(Cr + Ct);   %aspect ratio

%discretizing the span to get distinct algebraic equations
N = 30;                                        %number of discrete points
[j,theta] = meshgrid(1:N, linspace(0,pi,N));   %j = 1 to N, theta = location along span
y = -(b/2)*cos(theta);                         %coordinate transform as applied to Prandtl's equation
c = Cr + (2/b)*(Ct - Cr)*abs(y);               %length of chord at each distinct location 

aoa = deg2rad(4);                             %geometric angle of attack in radians
zero_lift_aoa = deg2rad(-6);                   %zero lift angle of attack
AOA = (aoa - zero_lift_aoa)*ones(N,1);         %RHS matrix to store angles of attack

%determining coefficients of algebraic equations
C = (2*b./(pi.*c)).*sin(j.*theta) + j.*sin(j.*theta)./sin(theta);   
for i = 1:N                                     
    C(1,i) = i*i;         %applying L'Hopital's rule at theta = 0 to avoid singularity    
end

%Results
A = C\AOA;              %determining fourier coefficients and solving the set of equations
CL = A(1)*pi*AR;        %calculating lift coefficient 
sigma = j(1,:)*(A.*A);      	
s = sum(sigma);             %calculating summation term
Cdi = pi*AR*s;        %induced drag coefficient 
disp(CL);
disp(Cdi); 
e = (CL*CL)/(pi*AR*Cdi);