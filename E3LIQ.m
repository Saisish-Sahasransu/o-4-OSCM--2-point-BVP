%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           (p(x)u')' + k^2u = f(x)                           %%
%%%       u[x=0.5] = cos(0.5)-sin(1)                            %%
%%%   [pu'][x=0.5] = -4*sin(0.5)-4*cos(1)                       %%
%%%     p(x) = 2 = P1 [0,0.5] //  p(x) = 4 = P2 ]0.5,1]         %%
%%%   Exact solution u(x)=sin(2*x)  [0,0.5] // cos(x)  ]0.5,1]  %%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
format short E;
x0 = 0;
xf = 1;
kmin = 100;
kmax = 100;
%N = input('The number of sub-intervals N = ');
for p = 1:5
    N(p) = 4*p;
h(p) = (xf-x0)/N(p);
for j = 1:N(p)+1
    x(j) = x0 + (j-1)*h(p);
end
p1 =-sqrt(3/5);
p2 = 0;
p3= sqrt(3/5);
    xi(1) = x0 + 0.5*(1 + p1)*h(p);
    xi(2) = x0 + 0.5*(1 + p2)*h(p);
    xi(3) = x0 + 0.5*(1 + p3)*h(p);
for j = 2:N(p) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     xi(3*j-2) = x(j) + 0.5*(1 + p1)*h(p);
     xi(3*j-1) = x(j) + 0.5*(1 + p2)*h(p);
     xi(3*j)   = x(j) + 0.5*(1 + p3)*h(p);
end
a1=2;
b1=4;
a = (1/2)*h(p)*(1+p1);
b = (1/2)*h(p)*(1+p2);
c = (1/2)*h(p)*(1+p3);
A = zeros(5*N(p)+2);

 A(1,1) = 1;
 A(1,2) = 1;
 A(5*N(p)+2,5*N(p)+2) = 1;
 A(5*N(p)+2,5*N(p)+1) = 1;
 for j = 1:N(p)/2
     A(5*j-3,5*j-4) = kmin;
     A(5*j-3,5*j-3) = a*kmin;
     A(5*j-3,5*j-2) = 2*a1 + kmin*a^2;
     A(5*j-3,5*j-1) = 6*a*a1 + kmin*a^3;
     A(5*j-3,5*j)   = kmin*a^4 + 12*a^2*a1;
     A(5*j-2,5*j-4) = kmin;
     A(5*j-2,5*j-3) = b*kmin;
     A(5*j-2,5*j-2) = 2*a1 + kmin*b^2;
     A(5*j-2,5*j-1) = 6*b*a1 + kmin*b^3;
     A(5*j-2,5*j)   = kmin*b^4 + 12*b^2*a1;
     A(5*j-1,5*j-4) = kmin;
     A(5*j-1,5*j-3) = c*kmin;
     A(5*j-1,5*j-2) = 2*a1 + kmin*c^2;
     A(5*j-1,5*j-1) = 6*c*a1 + kmin*c^3;
     A(5*j-1,5*j)   = kmin*c^4 + 12*c^2*a1;
     A(5*j,5*j-4)   = -1;
     A(5*j,5*j-3)   = -h(p);
     A(5*j,5*j-2)   = -h(p)^2;
     A(5*j,5*j-1)   = -h(p)^3;
     A(5*j,5*j)     = -h(p)^4;
     A(5*j,5*j+1)   = 1;
     A(5*j+1,5*j-4) = 0;
     A(5*j+1,5*j-3) = -1;
     A(5*j+1,5*j-2) = -2*h(p);
     A(5*j+1,5*j-1) = -3*h(p)^2;
     A(5*j+1,5*j)   = -4*h(p)^3;
     A(5*j+1,5*j+2) = 1;
 end
for j = N(p)/2+1:N(p)
     A(5*j-3,5*j-4) = kmax;
     A(5*j-3,5*j-3) = a*kmax;
     A(5*j-3,5*j-2) = 2*b1 + kmax*a^2;
     A(5*j-3,5*j-1) = 6*a*b1 + kmax*a^3;
     A(5*j-3,5*j)   = kmax*a^4 + 12*a^2*b1;
     A(5*j-2,5*j-4) = kmax;
     A(5*j-2,5*j-3) = b*kmax;
     A(5*j-2,5*j-2) = 2*b1 + kmax*b^2;
     A(5*j-2,5*j-1) = 6*b*b1 + kmax*b^3;
     A(5*j-2,5*j)   = kmax*b^4 + 12*b^2*b1;
     A(5*j-1,5*j-4) = kmax;
     A(5*j-1,5*j-3) = c*kmax;
     A(5*j-1,5*j-2) = 2*b1 + kmax*c^2;
     A(5*j-1,5*j-1) = 6*c*b1 + kmax*c^3;
     A(5*j-1,5*j)   = kmax*c^4 + 12*c^2*b1;
     A(5*j,5*j-4)   = -1;
     A(5*j,5*j-3)   = -h(p);
     A(5*j,5*j-2)   = -h(p)^2;
     A(5*j,5*j-1)   = -h(p)^3;
     A(5*j,5*j)     = -h(p)^4;
     A(5*j,5*j+1)   = 1;
     A(5*j+1,5*j-4) = 0;
     A(5*j+1,5*j-3) = -1;
     A(5*j+1,5*j-2) = -2*h(p);
     A(5*j+1,5*j-1) = -3*h(p)^2;
     A(5*j+1,5*j)   = -4*h(p)^3;
     A(5*j+1,5*j+2) = 1;
end
     A(5*N(p)/2+1,5*N(p)/2-4) = 0;
     A(5*N(p)/2+1,5*N(p)/2-3) = -a1;
     A(5*N(p)/2+1,5*N(p)/2-2) = -2*a1*h(p);
     A(5*N(p)/2+1,5*N(p)/2-1) = -3*a1*h(p)^2;
     A(5*N(p)/2+1,5*N(p)/2)   = -4*a1*h(p)^3;
     A(5*N(p)/2+1,5*N(p)/2+2) = b1;
b = zeros(5*N(p)+2,1);
b(1) = 2;
b(5*N(p)+2) = cos(1)-sin(1);
for j = 1:N(p)/2
    b(5*j-3) = (-4*a1+kmin)*sin(2*xi(3*j-2));
    b(5*j-2) = (-4*a1+kmin)*sin(2*xi(3*j-1));
    b(5*j-1) = (-4*a1+kmin)*sin(2*xi(3*j));
end
    b(5*N(p)/2) = cos(0.5)-sin(1);
    b(5*N(p)/2+1) = -4*sin(0.5)-4*cos(1);
for j = N(p)/2+1:N(p)
    b(5*j-3) = (-b1+kmax)*cos(xi(3*j-2));
    b(5*j-2) = (-b1+kmax)*cos(xi(3*j-1));
    b(5*j-1) = (-b1+kmax)*cos(xi(3*j));
end
s= A\b;
P = [h(p)/11 2*h(p)/11 3*h(p)/11 4*h(p)/11 5*h(p)/11 6*h(p)/11 7*h(p)/11 8*h(p)/11 9*h(p)/11 10*h(p)/11];
   
for j = 1:10
    for k = 1:N(p)
        U1(j,k) = s(5*k-4) + s(5*k-3)*P(j)+ s(5*k-2)*P(j)^2 + s(5*k-1)*P(j)^3+s(5*k)*P(j)^4;
    end
end
Z1 = U1';
Z2 = Z1(:);


for j = 1:10
   
    for k = 1:N(p)/2
        exact(j,k) = sin(2*(x(k)+P(j)));
    end
    for k = N(p)/2+1:N(p)
        exact(j,k) = cos((x(k)+P(j)));
    end 
end

Z3 = exact';
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

    