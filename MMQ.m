function [X,NPA] = MMQ(x,dT)
%-------------------------------------------------------------------------
% Least Squares Algorithm
%   (Luckett, 1975; Sachdev, 1976)

% a) Description of input variables:
%   x : Time samples matrix (samples x number of channels)
%   dT : intervalo de amostragem (s)
%
% b) Description of output variables:
%   X : Matrix with the calculated phasors (samples x number of channels)
%   NPA : Number of points in the sampling window

%--------------------------------------------------------------------------
F0 = 60;        % Fundamental frequency
JAN = 1;        % Window size (in cycles)
NREGDC = 2;     % Number of DC regressors [ 1  t  t^2  t^3  ... ]
NREGAC = 1;     % Number of AC regressors [ sin(wot)  cos(wot)   sin(3wot)  cos(3wot) ... ]
W0 = 2*pi*F0;
NPA = JAN/(F0*dT);
THETAC = W0*dT;         % Characteristic angle
[NTOT,NC] = size(x);
X = zeros(NTOT,NC);

% Calculation of the PSEUDO-INVERSE
R = [];         % Regressor matrix
A = [];         % Pseudo-inverse
AC = [];
DC = [];

for i = 1:NPA
    
    vAC = [];
    vDC = [];
    
    %Determine AC regressors
    for j = 1:NREGAC
        vAC = [  vAC, sin((2*j-1)*W0*i*dT), cos((2*j-1)*W0*i*dT) ] ;
    end
    %Determine DC regressors
    for j = 1:NREGDC
        vDC = [  vDC, (dT)^(j-1) ] ;
    end
    
    AC = [AC; vAC];
    DC = [DC; vDC];
    
end

R = [AC DC];

A = pinv(R);

% Calculate the phasors
for i = NPA:NTOT
    
    THETA = A*x((i-NPA+1):i,:);
    
    X(i,:) = THETA(1,:) + sqrt(-1)*THETA(2,:);
    
    % Perform angular correction
    X(i,:) = X(i,:)*exp(-sqrt(-1)*(i-1)*THETAC);
    
end
