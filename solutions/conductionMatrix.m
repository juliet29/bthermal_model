function A = conductionMatrix(lambda,dt,dx,Nx)

A = zeros(Nx,Nx);

% First and last row 
A(1,1)     = 1-lambda;
A(1,2)     = lambda;
A(Nx,Nx-1) = lambda;
A(Nx,Nx)   = 1-lambda;

% All other rows 
for i = 2:Nx-1
    A(i,i)   = 1-2*lambda;
    A(i,i-1) = lambda;
    A(i,i+1) = lambda; 
end