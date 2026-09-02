% settings is a structure that allows easily switching between cases tested
% in the simulations
%
% Author: Antoine Falisse
% Date: 9/9/2019
% settings(1): derivative supplier identifier
%   1: Algorithmic differentiation / Recorder  
%   2: Algorithmic differentiation / ADOL-C  
%   3: Finite differences
% settings(2): Hessian identifier
%   1: Approximated Hessian
%   2: Exact Hessian
% settings(3):  linear solver identifier
%   1: mumps
%   2: ma27
%   3: ma57
%   4: ma77
%   5: ma86
%   6: ma97
% settings(4): initial guess identifier
%   1: quasi-random initial guess  
%   2: data-informed initial guess (constant values for muscle variables)
%   3: data-informed initial guess (experimental values for muscle variables)
settings = [
    % A. Impact of the derivative supplier
    % Recorder - Approximated Hessian - mumps
    1,1,1,1;     %5    1
    1,1,1,2;     %6    2
    1,1,1,3;     %23   3
    % ADOL-C - Approximated Hessian - mumps
    2,1,1,1;     %21   4
    2,1,1,2;     %20   5
    2,1,1,3;     %22   6
    % FD - Approximated Hessian - mumps
    3,1,1,1;     %43   7
    3,1,1,2;     %44   8
    3,1,1,3;     %45   9
    % B. Impact of the linear solver
    % Recorder - Approximated Hessian - ma27 
    1,1,2,1;     %67   10
    1,1,2,2;     %70   11
    1,1,2,3;     %73   12
    % Recorder - Approximated Hessian - ma77 
    1,1,4,1;     %68   13
    1,1,4,2;     %71   14
    1,1,4,3;     %74   15
    % Recorder - Approximated Hessian - ma86 
    1,1,5,1;     %69   16
    1,1,5,2;     %72   17
    1,1,5,3;     %75   18    
    ];   
