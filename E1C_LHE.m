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
   P1 = -sqrt(3/5);
   P2 = 0;
   P3 = sqrt(3/5);
   c = [0.5*(1+P1)*h(p) 0.5*(1+P2)*h(p) 0.5*(1+P3)*h(p)];
   
   for k = 1:3
         for j = 1:N(p)
          U1(k,j) = s(4*j-3) + s(4*j-2)*(c(k)) + s(4*j-1)*(c(k))^2 + s(4*j)*(c(k))^3;
          U2(k,j) = s(4*j-2) + 2*s(4*j-1)*c(k) + 3*s(4*j)*(c(k))^2;
         end
   end
    Z = U1(:);  
    Z1 = U2(:);
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
    for k = 1:3
          for j = 1:N(p)/2
              exact2(k,j) = sin(pi*(x(j)+c(k)));
              exact3(k,j) = pi*cos(pi*(x(j)+c(k)));
          end
      end
      for k = 1:3
          for j = N(p)/2+1:N(p)
              exact2(k,j) = sin(pi*(x(j)+c(k)));
              exact3(k,j) = pi*cos(pi*(x(j)+c(k)));
          end
      end
Z2 = exact2(:);
Z3 = exact3(:);
sum1 = 0;
sum2 = 0;

for j = 1:3*N(p)
   sum1 = sum1  + (Z(j)-Z2(j))^2 + (Z1(j)-Z3(j))^2;
   sum2 = sum2  + (Z(j)-Z2(j))^2 ;
end

    
     error(p) = max(abs(exact-Asol))
     error1(p) = max(abs(exact1-Asol1))
     error2(p) = sqrt(h(p)*sum1)
     error3(p) = sqrt(h(p)*sum2)
end
     
 %%%%% Order of Convergence %%%%%%%
for j = 1:p-1
order(j) = log(error(j)/error(j+1))/log(h(j)/h(j+1))
end
for j = 1:p-1
order1(j) = log(error1(j)/error1(j+1))/log(h(j)/h(j+1))
end
for j = 1:p-1
order2(j) = log(error2(j)/error2(j+1))/log(h(j)/h(j+1))
end
for j = 1:p-1
order3(j) = log(error3(j)/error3(j+1))/log(h(j)/h(j+1))
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