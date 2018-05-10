clc
clear all
%Splits up the domain of interest. 2 pi is the length
%and width of the rectangle
N = 40;
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
d_t = d_x^2/4;

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
lam = .5*d_t/d_x^2;

%initial error value, set high to make sure condition is not met automatically
error = 10;

for t = 0:d_t:10
    U_init = U_solut;
for i = 2:N-1
    for j = 1:N-1
        %Neumannn conditions with ghost node
        if i >= 2 && j == 1 && i <= N-1
             U_solut(i,j) = lam*(2* U_init(i, j+1) + U_init(i+1,j) + U_init(i-1,j)) + (1-4*lam)*(U_init(i,j)); 
        %interior points
        else
             U_solut(i,j) = lam*(U_init(i,j-1)+U_init(i, j+1) + U_init(i+1,j) + U_init(i-1,j)) + (1-4*lam)*(U_init(i,j)); 
        end
    end

end
        %This finds the relative error
        error = abs(mean(mean(U_solut)) - mean(mean(U_init)))/abs(mean(mean(U_solut)));
       %This shows the time interval it took for the solution to change very little for each time step
       if error < .00001
            disp(t)
       break
       end
    surf(x,y,U_solut)
    drawnow
end
