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
k=100;
%N = input('The number of sub-intervals N = ');
for p = 1:5
    N(p) = 5*2^(p-1);
    h(p) = (xf-x0)/N(p);
for j = 1:N(p)+1
    x(j) = x0 + (j-1)*h(p);
end
p1 =-sqrt(3/5);
p2 = 0;
p3 = sqrt(3/5);
    xi(1) = x0 + 0.5*(1 + p1)*h(p);
    xi(2) = x0 + 0.5*(1 + p2)*h(p);
    xi(3) = x0 + 0.5*(1 + p3)*h(p);
for j = 2:N(p) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     xi(3*j-2) = x(j) + 0.5*(1 + p1)*h(p);
     xi(3*j-1) = x(j) + 0.5*(1 + p2)*h(p);
     xi(3*j)   = x(j) + 0.5*(1 + p3)*h(p);
end
a = (1/2)*h(p)*(1+p1);
b = (1/2)*h(p)*(1+p2);
c = (1/2)*h(p)*(1+p3);
A = zeros(5*N(p)+2);
P1=1;
P2=2;
P3=3;
 A(1,1) = 1;
 A(1,2) = 1;
 A(5*N(p)+2,5*N(p)+1) = 1;
 A(5*N(p)+2,5*N(p)+2) = 1;
 for j = 1:N(p)/5
     A(5*j-3,5*j-4) = k;
     A(5*j-3,5*j-3) = a*k;
     A(5*j-3,5*j-2) = 2*P1 + k*a^2;
     A(5*j-3,5*j-1) = 6*a*P1 + k*a^3;
     A(5*j-3,5*j)   = k*a^4 + 12*a^2*P1;
     A(5*j-2,5*j-4) = k;
     A(5*j-2,5*j-3) = b*k;
     A(5*j-2,5*j-2) = 2*P1 + k*b^2;
     A(5*j-2,5*j-1) = 6*b*P1 + k*b^3;
     A(5*j-2,5*j)   = k*b^4 + 12*b^2*P1;
     A(5*j-1,5*j-4) = k;
     A(5*j-1,5*j-3) = c*k;
     A(5*j-1,5*j-2) = 2*P1 + k*c^2;
     A(5*j-1,5*j-1) = 6*c*P1 + k*c^3;
     A(5*j-1,5*j)   = k*c^4 + 12*c^2*P1;
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
 
     A(5*j+1,5*j-3) = -P1;
     A(5*j+1,5*j-2) = -2*P1*h(p);
     A(5*j+1,5*j-1) = -3*P1*h(p)^2;
     A(5*j+1,5*j)   = -4*P1*h(p)^3;
     A(5*j+1,5*j+2) = P2;
for j = N(p)/5+1:3*N(p)/5
     A(5*j-3,5*j-4) = k;
     A(5*j-3,5*j-3) = a*k;
     A(5*j-3,5*j-2) = 2*P2 + k*a^2;
     A(5*j-3,5*j-1) = 6*a*P2 + k*a^3;
     A(5*j-3,5*j)   = k*a^4 + 12*a^2*P2;
     A(5*j-2,5*j-4) = k;
     A(5*j-2,5*j-3) = b*k;
     A(5*j-2,5*j-2) = 2*P2 + k*b^2;
     A(5*j-2,5*j-1) = 6*b*P2 + k*b^3;
     A(5*j-2,5*j)   = k*b^4 + 12*b^2*P2;
     A(5*j-1,5*j-4) = k;
     A(5*j-1,5*j-3) = c*k;
     A(5*j-1,5*j-2) = 2*P2 + k*c^2;
     A(5*j-1,5*j-1) = 6*c*P2 + k*c^3;
     A(5*j-1,5*j)   = k*c^4 + 12*c^2*P2;
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
     A(5*j+1,5*j-3) = -P2;
     A(5*j+1,5*j-2) = -2*P2*h(p);
     A(5*j+1,5*j-1) = -3*P2*h(p)^2;
     A(5*j+1,5*j)   = -4*P2*h(p)^3;
     A(5*j+1,5*j+2) = P3;
     
for j = 3*N(p)/5+1:N(p)
     A(5*j-3,5*j-4) = k;
     A(5*j-3,5*j-3) = a*k;
     A(5*j-3,5*j-2) = 2*P3 + k*a^2;
     A(5*j-3,5*j-1) = 6*a*P3 + k*a^3;
     A(5*j-3,5*j)   = k*a^4 + 12*a^2*P3;
     A(5*j-2,5*j-4) = k;
     A(5*j-2,5*j-3) = b*k;
     A(5*j-2,5*j-2) = 2*P3+ k*b^2;
     A(5*j-2,5*j-1) = 6*b*P3 + k*b^3;
     A(5*j-2,5*j)   = k*b^4 + 12*b^2*P3;
     A(5*j-1,5*j-4) = k;
     A(5*j-1,5*j-3) = c*k;
     A(5*j-1,5*j-2) = 2*P3+ k*c^2;
     A(5*j-1,5*j-1) = 6*c*P3 + k*c^3;
     A(5*j-1,5*j)   = k*c^4 + 12*c^2*P3;
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
b = zeros(5*N(p)+2,1);
b(1) = 1;
b(5*N(p)+2) = cos(1)-sin(1);
for j = 1:N(p)/5
    b(5*j-3) = (k-P1)*sin(xi(3*j-2));
    b(5*j-2) = (k-P1)*sin(xi(3*j-1));
    b(5*j-1) = (k-P1)*sin(xi(3*j));
end
    b(5*j) = exp(-1/25)-sin(1/5);
    b(5*j+1) = -2*0.2*exp(-0.04)*P2-P1*cos(0.2);
for j = (N(p)/5+1):3*N(p)/5
    b(5*j-3) = (-2*P2+4*P2*(xi(3*j-2))^2+k)*exp(-(xi(3*j-2))^2);
    b(5*j-2) = (-2*P2+4*P2*(xi(3*j-1))^2+k)*exp(-(xi(3*j-1))^2);
    b(5*j-1) = (-2*P2+4*P2*(xi(3*j))^2+k)*exp(-(xi(3*j))^2);
end
b(5*j) = cos(3/5)-exp(-9/25);
b(5*j+1) = P3*(-sin(0.6)) + 2*0.6*exp(-0.36)*P2;
for j = 3*2*(N(p)/10)+1:N(p)
    b(5*j-3) = (k-P3)*cos(xi(3*j-2));
    b(5*j-2) = (k-P3)*cos(xi(3*j-1));
    b(5*j-1) = (k-P3)*cos(xi(3*j));
    
end
sol = A\b;

%%%%%%%%%%%%%%%%solution%%%%%%%%%%%%%%%
    for j = 1:N(p)+1
        Asol(j) = sol(5*j-4);
    end
    for j = 1:N(p)/5
        exact(j) = sin(x(j));
    end
    for j = N(p)/5+1:3*N(p)/5
        exact(j) = exp(-(x(j))^2);
    end 
     for j = 3*N(p)/5+1:N(p)+1
        exact(j) = cos(x(j));
    end 
%%%%%%%%%%%%%%%%Derrivative%%%%%%%%%%%%%%%  
    for j = 1:N(p)+1
        Asol1(j) = sol(5*j-3);
    end
    for j = 1:N(p)/5
        exact1(j) = cos(x(j));
    end
    for j = N(p)/5+1:3*N(p)/5
        exact1(j) =-2*x(j)*exp(-(x(j))^2);
    end 
    for j = 3*N(p)/5+1:N(p)+1
        exact1(j) = -sin(x(j));
    end
     
   
     error(p)  = max(abs(exact-Asol))
     error1(p) = max(abs(exact1-Asol1))
end
     
 %%%%% Order of Convergence %%%%%%%
for j = 1:p-1
order(j) = log(error(j)/error(j+1))/log(h(j)/h(j+1))
end
for j = 1:p-1
order1(j) = log(error1(j)/error1(j+1))/log(h(j)/h(j+1))
end

%%% Order of Convergence Graphically %%%%%
figure(1)
plot(log(h),log(error),'*',log(h),log(error));
xlabel('log(h_i)')
ylabel('log(e_i)')
title('Order of Convergence Graphically')

figure(2)
plot(x,Asol)
xlabel('x')
ylabel('Y')
title('Numerically Computed Soln')
    
figure(3)
plot(x,exact)
xlabel('x')
ylabel('y')
title('Exact Soln')
