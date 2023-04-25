function [J] = JacobianFDSecondOrder(f,x,solverOptions,dx)
    arguments
        f
        x
        solverOptions                       =[];
        dx                                  =1e-4;
    end
%JACOBIANFDSECONDORDER Finds Jacobian via finite differencing
%   Uses a fixed displacement dx in each of the DOFs with the central
%   finite difference to get a second-order accurate estimate of the
%   Jacobian, which is the gradient of a column vector of functions.
%
%   INPUTS
%   f - m*1 function vector that takes n dimensional input
%   x - n*1 initial guess of free variables
%   solverOptions - struct for optional passed variables to f
%   [dx] - Optional perturbation value for the finite differencing
%
%   OUTPUTS
%   J - The m*n Jacobian of f

  % Size of DOF
  n = length(x);
  % Size of f
  m = size(f(x,solverOptions),1);
  % initialize shifted x value which we will use in the finite difference
  delx = zeros(n,1);
  % Size of mapped range
  J = zeros(m,n);
  
  for i = 1:n
    % Difference vector dx * e_i
    delx(i) = dx;

    % Populate J
    J(:,i) = (f(x+delx,solverOptions)-f(x-delx,solverOptions))/2/dx;
  
    % reset shifted x value
    delx(i) = 0;
  end
end