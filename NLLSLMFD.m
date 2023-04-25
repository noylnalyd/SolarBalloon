function [x, res, conv, iter, xs, ress, convs] = NLLSLMFD(f, guess, tol, max_iter, solverOptions, omega, gamma, g, rrefCheck)
    arguments
        f
        guess
        tol                                 =1e-3;
        max_iter(1,1) {mustBePositive}      =30;
        solverOptions (1,1) struct          =[];
        omega (1,1) double                  =1;
        gamma (1,1) double                  =.5;
        g (1,1) double                      =.85;
        rrefCheck (1,1)                     =false;
    end
%NLLSLMFD Nonlinear Least-Squares Levenburg-Marquardt Finite Difference
%MIMO Solver
%   Method to solve a function vector by linearizing at each iteration
%   about a guess input. The algorithm used is Levenburg-Marquardt which
%   is gradient descent (GD) until the constant c is insignificant
%   and then switches to newton-raphson (NR). The only difference between
%   the two methods is the step size: GD uses an exponentially decreasing
%   value gamma*g^k, NR uses the full magnitude of omega*J^{-1}f
%
%   In the case where there are more/less functions than DOFs (i.e. m!=n),
%   this method uses a least-squares approach based on the Gauss-Newton
%   method. Here:
%       <latex>
%       \begin{align*}
%       f(x) &= 0\\
%       J(x_{k}) (x_{k+1}-x_{k}) &= - c \cdot  f(x_{k})\\
%       J^{T}(x_{n})J(x_{n})(x_{n+1}-x_{n}) &= - c \cdot J^{T}(x_{k}) f(x_{k})\\
%       x_{k+1} &= x_{k} + \left(J^{T}(x_{k})J(x_{k}) \right)^{-1} \left(- c \cdot J^{T}(x_{n}) f(x_{k})\right)\\
%       \end{align*}
%       </latex>
%
%   INPUTS
%   f - m*1 function vector that takes n dimensional input
%   guess - n*1 initial guess
%   [tol] - residual and convergence tolerance
%   [max_iter] - max number of iterations permitted
%   [solverOptions] - Optional struct object that contains passed values to f
%   [omega] - SOR filter / relaxation parameter
%   [gamma] - initial radius of step
%   [g] - g exponentiation parameter of step size
%   [rrefCheck] - bool to check matrices for full rankness (Set false for
%           large problems!)
%
%   OUTPUTS
%   x - n*1 solution
%   res - 1*1 final residual
%   conv - 1*1 final iterative convergence
%   [iters] - number of iterations til convergence / divergence
%   [xs] - iter*1 cell array of n*1 solutions
%   [ress] - iter*1 residuals
%   [convs] - iter*1 convergences

    n = size(guess,1);
    m = size(f(guess,solverOptions),1);

    isAllIters =  nargout > 4;

    isLS = m~=n;

    x = guess; % n*1, upheld solution
    iter = 1;  % initialize iteration counter
    res = 0;
    conv = NaN;
    ress(1) = res;
    convs(1) = conv;
    if isAllIters
        xs = cell(1); % kth solutions starting after initial guess
    else
        xs = NaN;
    end
    
    while iter < max_iter
        % save solution from previous iteration
        xold = x; % n * 1

        % evaluate the Jacobian
        J = JacobianFDSecondOrder(f,x,solverOptions); % m * n
        
        % check Jacobian validity
        if sum(sum(isnan(J(:,:))))>0 || norm(J,1)<1e-5
            warning("Final Jacobian has NANs or indeterminate.")
            return;
        end
        % deeper check
        if rrefCheck && rank(J,1e-7) < m
            warning("Final Jacobian not full rank.")
            return;
        end
        
        if isLS
            LHS = transpose(J)*J;
            RHS = transpose(J)*f(x,solverOptions);
        else
            LHS = J;
            RHS = f(x,solverOptions);
        end
        
        % solve for update vector
        % If not well-conditioned, use this as final solution
        if cond(LHS) < 1e10
            delx = -LHS \ RHS; % n
        else
            warning('Poorly conditioned last iteration in NLLSLMFD')
            return;
        end
        % Norm of delx
        ndelx = norm(delx,2);
        % if too large, switch to GD. Else, already NR
        if gamma < ndelx
            % scale to gamma length
            delx = delx * gamma / ndelx;
        else
            % scale by NR relaxation
            delx = delx * omega;
        end

        x = xold+delx;

        % calculate residual and iterative convergence
        res = norm(f(x,solverOptions),2);
        conv = ndelx;
        if isAllIters
            ress(iter) = res;
            convs(iter) = conv;
            xs{iter} = x;
        end

        % check for convergence
        if (conv < tol && iter>1) && ( isLS || res<tol )
          return;
        end
        
        % check for divergence
        if res > 1e10
          return;
        end

        gamma = gamma * g; % scale down gamma step
        iter=iter+1;  % advance iteration counter         
    end
end