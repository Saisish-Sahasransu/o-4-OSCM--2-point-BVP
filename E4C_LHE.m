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
for p = 1:6
    N(p) = 10*2^(p-1);
    h(p) = (xf-x0)/N(p);
for j = 1:N(p)+1
    x(j) = x0 + (j-1)*h(p);
end
p1 = -1/sqrt(3);
p2 = 1/sqrt(3);
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

b1(4*(2*(N(p)/10)))   = exp(-1/25)-sin(1/5);
b1(4*(2*(N(p)/10))+1) = -2*0.2*exp(-0.04)*P2-P1*cos(0.2);
for j = N(p)/5+1:3*(N(p)/5)
    b1(4*j-2) = (-2*P2+4*P2*(xi(2*j-1))^2+k)*exp(-(xi(2*j-1))^2);
    b1(4*j-1) = (-2*P2+4*P2*(xi(2*j))^2+k)*exp(-(xi(2*j))^2);
    
end

b1(4*(3*2*(N(p)/10))) = cos(3/5)-exp(-9/25);
b1(4*(3*2*(N(p)/10))+1) = P3*(-sin(0.6)) + 2*0.6*exp(-0.36)*P2;
for j = 3*2*(N(p)/10)+1:N(p)
    b1(4*j-2) = (k-P3)*cos(xi(2*j-1));
    b1(4*j-1) = (k-P3)*cos(xi(2*j));
    
end
s = A\b1;
   PP1 = -sqrt(3/5);
   PP2 = 0;
   PP3 = sqrt(3/5);
   c = [0.5*(1+PP1)*h(p) 0.5*(1+PP2)*h(p) 0.5*(1+PP3)*h(p)];
   
   for k = 1:3
         for j = 1:N(p)
          U1(k,j) = s(4*j-3) + s(4*j-2)*(c(k)) + s(4*j-1)*(c(k))^2 + s(4*j)*(c(k))^3;
          U2(k,j) = s(4*j-2) + 2*s(4*j-1)*c(k) + 3*s(4*j)*(c(k))^2;
         end
   end
    Z = U1(:);  
    Z1 = U2(:);
    for j = 1:N(p)/5
        exact(j) = sin(x(j));
    end
    for j = N(p)/5+1:3*N(p)/5
        exact(j) = exp(-(x(j))^2);
    end
    for j = 3*N(p)/5+1:N(p)+1
        exact(j) = cos(x(j));
    end
    
    for j = 1:N(p)+1
        Asol(j) = s(4*j-3);
    end
     for j = 1:N(p)/5
        exact1(j) = cos(x(j));
    end
    for j = N(p)/5+1:3*N(p)/5
        exact1(j) = -2*x(j)*exp(-(x(j))^2);
    end
    for j = 3*N(p)/5+1:N(p)+1
        exact1(j) = -sin(x(j));
    end
    
    
    for j = 1:N(p)+1
        Asol1(j) = s(4*j-2);
    end
    for k = 1:3
          for j = 1:N(p)/5
              exact2(k,j) = sin((x(j)+c(k)));
              exact3(k,j) = cos((x(j)+c(k)));
          end
      end
     for k = 1:3
          for j = N(p)/5+1:3*N(p)/5
              exact2(k,j) = exp(-(x(j)+c(k))^2);
              exact3(k,j) = -2*(x(j)+c(k))*exp(-(x(j)+c(k))^2);
          end
     end
      for k = 1:3
          for j = 3*N(p)/5+1:N(p)
              exact2(k,j) = cos((x(j)+c(k)));
              exact3(k,j) = -sin((x(j)+c(k)));
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

    
     error(p)  = max(abs(exact-Asol))
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