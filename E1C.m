%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           u'' + k^2u = f(x)                                 %%
%%%   u(-1)=0, u'(1)=-pi                                        %%
%%%   Exact Solution u(x)=sin(pi*x)                             %%
%%%                                                             %%
%%%                                                             %%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
format short E;
x0 = -1;
xf = 1;
kmin = 1;
kmax = 100;
%N = input('The number of sub-intervals N = ');
for p = 1:4
    N(p) = 4*2^p;
h(p) = (xf-x0)/N(p);
for j = 1:N(p)+1
    x(j) = x0 + (j-1)*h(p);
end
p1 = -1/sqrt(3);
p2 = 1/sqrt(3);
    xi(1) = x0 + 0.5*(1 + p1)*h(p);
    xi(2) = x0 + 0.5*(1 + p2)*h(p);
for j = 2:N(p)
     xi(2*j-1) = x(j) + 0.5*(1 + p1)*h(p);
     xi(2*j) = x(j) + 0.5*(1 + p2)*h(p);
end
a = (1/2)*h(p)*(1+p1);
b = (1/2)*h(p)*(1+p2);
A = zeros(4*N(p)+2);

 A(1,1) = 1;
 A(4*N(p)+2,4*N(p)+2) = 1;
 for j = 1:N(p)/2
     A(4*j-2,4*j-3) = kmin;
     A(4*j-2,4*j-2) = a*kmin;
     A(4*j-2,4*j-1) = 2 + kmin*a^2;
     A(4*j-2,4*j) = 6*a + kmin*a^3;
     A(4*j-1,4*j-3) = kmin;
     A(4*j-1,4*j-2) = b*kmin;
     A(4*j-1,4*j-1) = 2 + kmin*b^2;
     A(4*j-1,4*j) = 6*b + kmin*b^3;
     A(4*j,4*j-3) = -1;
     A(4*j,4*j-2) = -h(p);
     A(4*j,4*j-1) = -h(p)^2;
     A(4*j,4*j) = -h(p)^3;
     A(4*j,4*j+1) = 1;
     A(4*j+1,4*j-3) = 0;
     A(4*j+1,4*j-2) = -1;
     A(4*j+1,4*j-1) = -2*h(p);
     A(4*j+1,4*j) = -3*h(p)^2;
     A(4*j+1,4*j+2) = 1;
 end
for j = N(p)/2+1:N(p)
     A(4*j-2,4*j-3) = kmax;
     A(4*j-2,4*j-2) = a*kmax;
     A(4*j-2,4*j-1) = 2 + kmax*a^2;
     A(4*j-2,4*j) = 6*a + kmax*a^3;
     A(4*j-1,4*j-3) = kmax;
     A(4*j-1,4*j-2) = b*kmax;
     A(4*j-1,4*j-1) = 2 + kmax*b^2;
     A(4*j-1,4*j) = 6*b + kmax*b^3;
     A(4*j,4*j-3) = -1;
     A(4*j,4*j-2) = -h(p);
     A(4*j,4*j-1) = -h(p)^2;
     A(4*j,4*j) = -h(p)^3;
     A(4*j,4*j+1) = 1;
     A(4*j+1,4*j-3) = 0;
     A(4*j+1,4*j-2) = -1;
     A(4*j+1,4*j-1) = -2*h(p);
     A(4*j+1,4*j) = -3*h(p)^2;
     A(4*j+1,4*j+2) = 1;
end
b = zeros(4*N(p)+2,1);
b(1) = 0;
b(4*N(p)+2) = -pi;
for j = 1:N(p)/2
    b(4*j-2) = (-pi^2+kmin)*sin(pi*xi(2*j-1));
    b(4*j-1) = (-pi^2+kmin)*sin(pi*xi(2*j));
    b(4*j) = 0;
    b(4*j+1) = 0;
end
for j = N(p)/2+1:N(p)
    b(4*j-2) = (-pi^2+kmax)*sin(pi*xi(2*j-1));
    b(4*j-1) = (-pi^2+kmax)*sin(pi*xi(2*j));
    b(4*j) = 0;
    b(4*j+1) = 0;
end
s = A\b;
   
    for j = 1:N(p)+1
        exact(j) = sin(pi*x(j));
    end
    
    for j = 1:N(p)+1
        Asol(j) = s(4*j-3);
    end
    for j = 1:N(p)+1
        exact1(j) = pi*cos(pi*x(j));
    end
    
    for j = 1:N(p)+1
        Asol1(j) = s(4*j-2);
    end
    error(p) = max(abs(exact-Asol))
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