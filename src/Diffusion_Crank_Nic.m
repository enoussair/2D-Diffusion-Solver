%Splits up the domain of interest. 2 pi is the length
%and width of the rectangle
N = 60;
d_x = 2*pi/(N-1);
d_y = 2*pi/(N-1);

a_x = -pi; %left edge
b_x = pi; %right edge
a_y = -pi; %bottom edge
b_y = pi; %top edge

x = a_x:d_x:b_x; 
y = a_y:d_y:b_y;

%Time domain

d_t = d_x^2/4;

%The dirchelet boundary conditions
f_a = (x-a_x).^2.*cos(pi*x/(a_x)); %top boundary 
g_a = x.*(x-a_x).^2; %bottom edge

%The initial solution at t = 0
U_init = zeros(N,N);
U_init(N,:) = f_a;
U_init(1,:) = g_a; 
U_init(:,N) = g_a(N) + (y-a_y)/(b_y - a_y)*(f_a(N) - g_a(N));

%Setting up the diag matrix
lam = .5*d_t/(d_x)^2;
LHS = zeros((N-1)*(N-2),(N-1)*(N-2));

r = 1 + 4 * lam;

%Determines where Neumann conditions are present, which happen on left
%boundary condition at every N-1 node
check = 0;

%This value will contain the norm of the vector at specific Node value
grid = 0;

%This sets up with 5 Band System
for j = 1:(N-2)*(N-1) %rows
    for k = 1:(N-1)*(N-2) %columns
        if (j == k) %begins creating the diagnols
            LHS(j,k) = r; %Creates the main diagnol
            if (j+N-1) <= (N-1)*(N-2)%Checks to see if it can add to diagnol N-1 ahead of main diagnol
                LHS(j,k+N-1) = -lam;
            end
            if (k ~=1) && k-(N-1) >0  %checks to see if it can add to diagnol N-1 behing main diagnol
                LHS(j,k-(N-1)) = -lam;
            end 
            if (k ~=1)
                LHS(j,j-1) = -lam;
            end
            if k < (N-1)*(N-2)
                LHS(j,j+1) = - lam;
            end
        end
    end
        %This part takes into account Neumann conditions on left edge by
        %checking if the code is at the left edge, which would modify the
        %coefficients around the main diagnol
            if j == 1+check*(N-1)
            LHS(j,j+1) = -2*lam;
            if j > 1
                LHS(j,j-1) = 0;
                LHS(j-1,j) = 0;
            end
            check = check + 1; %increments the check value once a Neumann row has been found
        end
end


%Creates the RHS and implements the solution by solving LHS * U_unknowns =
%RHS
RHS = zeros((N-1)*(N-2),1);

%Initial initializtion of U_solut. U_solut will be compared to U_init,
%which will play the role of the previous solution in the upcoming for loop
%for the sake of comparing the new solutions to the previous ones to
%determine convegence.
U_solut = U_init; 

%The error variable will be used to determine the final error
error = 100; %percent (relative error, assumed to be 100 at first)

%This keeps track of how many iterations into iter
iter = 1;

for t = 0:d_t:4
n = 1;
for i = 2:(N-1)
    for j = 1:(N-1)
        %This if statement finds the solutions within the interior points
        if j > 1 && i > 2 && j < (N-1) && i < (N-1)
            RHS(n) = lam*(U_solut(i,j-1) + U_solut(i, j+1) + U_solut(i+1,j) + U_solut(i-1,j)) + (1-4*lam)*(U_solut(i,j)); 
        end 
        
        %This code checks to see whether the code has reached one of the
        %nodes near the Dirchelet boundary conditions
        
        %This makes sure that the column is higher than one, meaning that
        %the Neumann condition is not at play
        if j > 1
            %If the code is currently looking at the second row, it
            %means the bottom boundary condition must be added to the RHS
            %twice since we know it's a constant
            if (i == 2) && j <N-1
                RHS(n) = lam*(U_solut(i,j-1) + U_solut(i, j+1) + U_solut(i+1,j) + U_solut(i-1,j)) + (1-4*lam)*(U_solut(i,j)) + lam * U_solut(i-1,j);
            
            %If the code is at the second to last row but before the last column, it means the top boundary must be added    
            elseif i == N-1 && j < N-1
                RHS(n) = lam*(U_solut(i,j-1) + U_solut(i, j+1) + U_solut(i+1,j) + U_solut(i-1,j)) + (1-4*lam)*(U_solut(i,j)) + lam * U_solut(i+1,j);
            
            %If the code is at the bottom right corner edge, then two boundary conditions will be added: The right
            %edge boundary and the bottom edge boundary
            elseif j == N-1 && i == 2
            	RHS(n) = lam*(U_solut(i,j-1) + U_solut(i, j+1) + U_solut(i+1,j) + 2*U_solut(i-1,j)) + (1-4*lam)*(U_solut(i,j)) + lam * U_solut(i,j+1);
            
            %If the code is at the second to last column, but not at a corner,
            %then the right boundary must be added
            elseif j == N-1 && i > 2 && i < N-1
            	RHS(n) = lam*(U_solut(i,j-1) + 2*U_solut(i, j+1) + U_solut(i+1,j) + U_solut(i-1,j)) + (1-4*lam)*(U_solut(i,j));
            %If the code is at the top right corner
            %then the right boundary must be added along with the top boundary
            elseif j == N-1 && i == N-2
            	RHS(n) = lam*(U_solut(i,j-1) + 2*U_solut(i, j+1) + 2*U_solut(i+1,j) + U_solut(i-1,j)) + (1-4*lam)*(U_solut(i,j));
            end
        end
        
        %This portion uses the Neumann condition, dU/dx = 0 |x=-pi
        
        %This checks if the code is at the bottom left corner
        if j == 1 && i == 2
            RHS(n) = lam*(U_solut(i+1,j) + 2*U_solut(i-1,j)) + (1-4*lam)*(U_solut(i,j)) +2* lam * U_solut(i,j+1);
        
        %This is for when the code is within the left boundary condition
        elseif j == 1 && i > 2 && i < N-1
            RHS(n) = lam*(U_solut(i+1,j) + U_solut(i-1,j)) + (1-4*lam)*(U_solut(i,j)) +2* lam * U_solut(i,j+1);
            
        %This checks to see if it is at the top left corner    
        elseif j == 1 && i == N-1
            RHS(n) = lam*(2*U_solut(i+1,j) + U_solut(i-1,j)) + (1-4*lam)*(U_solut(i,j)) +2* lam * U_solut(i,j+1);
        end 
      n = n + 1;
    end
end



    %Solution for each time step

    U_solved = linsolve(LHS,RHS);
    %This breaks the solution, which is initially a vector, into a matrix
    %for easy meshing
    U_new = vec2mat(U_solved,N-1);
   
    %This injects the solution into the U_solut matrix, which already
    %contains the set boundary conditions.
    U_solut(2:N-1,1:N-1) = U_new;
    
    %This finds the relative error
    error = abs(mean(mean(U_solut)) - mean(mean(U_init)))/abs(mean(mean(U_solut)));
    
    %Once the error between the current and previous solution has been
    %found, the previous solution, U_init, becomes the new solution and the new solution,U_solut, goes
    %on to to change again.
    U_init = U_solut;
    iter = iter+1
    surf(x,y,U_solut)
    drawnow
    
    %If the average error is low, it means the equation is converging to a value and the for loop stops and the current time is shown

end
 grid = norm(U_solut,2)
