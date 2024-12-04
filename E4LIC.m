%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           (p(x)u')' + k^2u = f(x)                                                 %%
%%% u[x=0.2] = exp(-1/25)-sin(1/5)[pu'][x=0.2] = -2*0.2*exp(-0.04)*P2-P1*cos(0.2)    %%
%%% u[x=0.6] = cos(3/5)-exp(-9/25)[pu'][x=0.6] = P3*(-sin(0.6)) + 2*0.6*exp(-0.36)*P2 %%
%%%     p(x) = 1 = P1 [0,0.2] //  p(x) = 2 = P2 ]0.2,0.6] p(x) = 3 = P3 ]0.6,1]       %%
%%%   Exact solution u(x)=sin(x)  [0,0.2] // exp(-x^2)  ]0.2,0.6]  // cos(x) ]0.6,1]  %%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
format short E;
x0 = 0;
xf = 1;
k = 30;
%N = input('The number of sub-intervals N = ');
for p = 1:7
    N(p) = 55*2^(p-1);
    h(p) = (xf-x0)/N(p);
for j = 1:N(p)+1
    x(j) = x0 + (j-1)*h(p);
end
p1 = -1/sqrt(3);
p2 =  1/sqrt(3);
    xi(1) = x0 + 0.5*(1 + p1)*h(p);
    xi(2) = x0 + 0.5*(1 + p2)*h(p);
for j = 2:N(p)
     xi(2*j-1) =   x(j) + 0.5*(1 + p1)*h(p);
     xi(2*j)   =   x(j) + 0.5*(1 + p2)*h(p);
end
a = (1/2)*h(p)*(1+p1);
b = (1/2)*h(p)*(1+p2);
A = zeros(4*N(p)+2);
A(1,1) = 1;
A(1,2) = 1;
A(4*N(p)+2,4*N(p)+1) = 1;
A(4*N(p)+2,4*N(p)+2) = 1;
P1=1;
P2=2;
P3=3;
for j = 1:2*(N(p)/10)
     A(4*j-2,4*j-3) = k;
     A(4*j-2,4*j-2) = k*a;
     A(4*j-2,4*j-1) = 2*P1   +k*a^2;
     A(4*j-2,4*j)   = 6*P1*a +k*a^3;
     A(4*j-1,4*j-3) = k;
     A(4*j-1,4*j-2) = k*b;
     A(4*j-1,4*j-1) = 2*P1   +k*b^2;
     A(4*j-1,4*j)   = 6*P1*b +k*b^3;
     A(4*j,4*j-3)   = -1;
     A(4*j,4*j-2)   = -h(p);
     A(4*j,4*j-1)   = -h(p)^2;
     A(4*j,4*j)     = -h(p)^3;
     A(4*j,4*j+1)   = 1;
     A(4*j+1,4*j-3) = 0;
     A(4*j+1,4*j-2) = -1;
     A(4*j+1,4*j-1) = -2*h(p);
     A(4*j+1,4*j)   = -3*h(p)^2;
     A(4*j+1,4*j+2) = 1;
end
     A(4*j+1,4*j-2) = -P1;
     A(4*j+1,4*j-1) = -2*P1*h(p);
     A(4*j+1,4*j)   = -3*P1*h(p)^2;
     A(4*j+1,4*j+2) = P2;
for j =N(p)/5+1:3*N(p)/5
     A(4*j-2,4*j-3) = k;
     A(4*j-2,4*j-2) = k*a;
     A(4*j-2,4*j-1) = 2*P2   +k*a^2;
     A(4*j-2,4*j)   = 6*P2*a +k*a^3;
     A(4*j-1,4*j-3) = k;
     A(4*j-1,4*j-2) = k*b;
     A(4*j-1,4*j-1) = 2*P2   +k*b^2;
     A(4*j-1,4*j)   = 6*P2*b +k*b^3;
     A(4*j,4*j-3)   = -1;
     A(4*j,4*j-2)   = -h(p);
     A(4*j,4*j-1)   = -h(p)^2;
     A(4*j,4*j)     = -h(p)^3;
     A(4*j,4*j+1)   = 1;
     A(4*j+1,4*j-3) = 0;
     A(4*j+1,4*j-2) = -1;
     A(4*j+1,4*j-1) = -2*h(p);
     A(4*j+1,4*j)   = -3*h(p)^2;
     A(4*j+1,4*j+2) = 1;
end 
     A(4*j+1,4*j-2) = -P2;
     A(4*j+1,4*j-1) = -2*P2*h(p);
     A(4*j+1,4*j)   = -3*P2*h(p)^2;
     A(4*j+1,4*j+2) = P3;
for  j = 3*2*(N(p)/10)+1:N(p)
     A(4*j-2,4*j-3) = k;
     A(4*j-2,4*j-2) = k*a;
     A(4*j-2,4*j-1) = 2*P3   +k*a^2;
     A(4*j-2,4*j)   = 6*P3*a +k*a^3;
     A(4*j-1,4*j-3) = k;
     A(4*j-1,4*j-2) = k*b;
     A(4*j-1,4*j-1) = 2*P3   +k*b^2;
     A(4*j-1,4*j)   = 6*P3*b +k*b^3;
     A(4*j,4*j-3)   = -1;
     A(4*j,4*j-2)   = -h(p);
     A(4*j,4*j-1)   = -h(p)^2;
     A(4*j,4*j)     = -h(p)^3;
     A(4*j,4*j+1)   = 1;
     A(4*j+1,4*j-3) = 0;
     A(4*j+1,4*j-2) = -1;
     A(4*j+1,4*j-1) = -2*h(p);
     A(4*j+1,4*j)   = -3*h(p)^2;
     A(4*j+1,4*j+2) = 1;
 end 
     
b1 = zeros(4*N(p)+2,1);
b1(1) = 1;
b1(4*N(p)+2) = cos(1)-sin(1);
for j = 1:N(p)/5
    b1(4*j-2) = (k-P1)*sin(xi(2*j-1));
    b1(4*j-1) = (k-P1)*sin(xi(2*j));
   
end

b1(4*j)   = exp(-0.04)-sin(0.2);
b1(4*j+1) = -2*0.2*exp(-0.04)*P2-P1*cos(0.2);
for j = N(p)/5+1:3*(N(p)/5)
    b1(4*j-2) = (-2*P2+4*P2*(xi(2*j-1))^2+k)*exp(-(xi(2*j-1))^2);
    b1(4*j-1) = (-2*P2+4*P2*(xi(2*j))^2+k)*exp(-(xi(2*j))^2);
    
end

b1(4*j) = cos(3/5)-exp(-9/25);
b1(4*j+1) = P3*(-sin(0.6)) + 2*0.6*exp(-0.36)*P2;
for j = 3*2*(N(p)/10)+1:N(p)
    b1(4*j-2) = (k-P3)*cos(xi(2*j-1));
    b1(4*j-1) = (k-P3)*cos(xi(2*j));
    
end
s = A\b1;

PK = [h(p)/11 2*h(p)/11 3*h(p)/11 4*h(p)/11 5*h(p)/11 6*h(p)/11 7*h(p)/11 8*h(p)/11 9*h(p)/11 10*h(p)/11];
   
for j = 1:10
    for k = 1:N(p)
        U1(j,k) = s(4*k-3) + s(4*k-2)*PK(j)+ s(4*k-1)*PK(j)^2 + s(4*k)*PK(j)^3;
    end
end
Z1 = U1';
Z2 = Z1(:);

for j = 1:10
    for k = 1:N(p)/5
        exact1(j,k) = sin((x(k)+PK(j)));
    end
end
for j = 1:10
    for k = N(p)/5+1:3*N(p)/5
        exact1(j,k) = exp(-((x(k)+PK(j))^2));
    end
end
for j = 1:10
    for k = 3*N(p)/5+1:N(p)
        exact1(j,k) = cos((x(k)+PK(j)));
    end
end
Z3 = exact1';
Z4 = Z3(:);

error(p) = max(abs(Z2-Z4))
end
%%%%% Order of Convergence %%%%%%%
for j = 1:p-1
order(j) = log(error(j)/error(j+1))/log(h(j)/h(j+1))
end

% figure(1)
% plot(x,Asol)
% xlabel('x')
% ylabel('Y')
% title('Numerically Computed Soln')
%     
% figure(2)
% plot(x,exact)
% xlabel('x')
% ylabel('y')
% title('Exact Soln')
%  
% figure(3)
% plot(x,exact,'o',x,Asol)
