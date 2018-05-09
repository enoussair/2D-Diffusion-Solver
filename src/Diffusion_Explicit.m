clc
clear all
%Splits up the domain of interest. 2 pi is the length
%and width of the rectangle
N = 1000;
d_x = 2*pi/(N-1);
d_y = 2*pi/(N-1);

a_x = -pi; %left edge
b_x = pi; %right edge
a_y = -pi; %bottom edge
b_y = pi; %top edge

x = a_x:d_x:b_x; 
y = a_y:d_y:b_y;

%Time domain
t = 4; %seconds
d_t = d_x^2/(8*N-1);

%The dirchelet boundary conditions
f_a = (x-a_x).^2.*cos(pi*x/(a_x)); %top boundary 
g_a = x.*(x-a_x).^2; %bottom edge

%The initial solution at t = 0
U_init = zeros(N,N);
U_init(N,:) = f_a;
U_init(1,:) = g_a; 
U_init(:,N) = g_a(N) + (y-a_y)/(b_y - a_y)*(f_a(N) - g_a(N));

%Explicit Method
U_solut = U_init;
lam = d_t/d_x^2;

for i = 1:d_t:4
    U_init = U_solut;
for i = 2:N-1
    for j = 1:N-1
        if i == 2 && j > 1
             U_solut(i,j) = lam*(U_init(i,j-1) + U_init(i, j+1) + U_init(i+1,j) + U_init(i-1,j)) + (1-4*lam)*(U_init(i,j)); 
        elseif j == 1 && i > 2
             U_solut(i,j) = lam*(2*U_init(i, j+1) + U_init(i+1,j) + U_init(i-1,j)) + (1-4*lam)*(U_init(i,j)); 
        elseif i == N-1 && j > 1
             U_solut(i,j) = lam*(U_init(i,j-1) + U_init(i, j+1) + U_init(i+1,j) + U_init(i-1,j)) + (1-4*lam)*(U_init(i,j));  
        elseif i > 2 && j == N-1
             U_solut(i,j) = lam*(U_init(i,j-1) + U_init(i, j+1) + U_init(i+1,j) + U_init(i-1,j)) + (1-4*lam)*(U_init(i,j));       
        end
    end
end
end
mesh(x,y,U_solut)