function [ E1, E2, E0 ] = calc_ce_mse(p, Dd, Dp, H, A, numPaths )
%
% Analytical MSE computation
% Dd, Dp are data and pilot matrices respectively
% A: transmitter matrix
% 
% Author: Shahab Ehsanfar


%% Initialization
L=numPaths;
Dd(:,1) = zeros(p.K,1);

F = get_DFT_matrix(p);

N = p.M*p.K;
dd = Dd(:);
dp = Dp(:);

% Create the Fourier matrix for truncation
Ft = F;
for i = L+1:N
   Ft(:,i) = zeros(N,1);
end

Xd = (sqrt(N)*F*A*dd);
Xp = (sqrt(N)*F*A*dp);

U = Ft*F';

Z = U'*U;

Xp1 = 1./Xp;

%% Arrangement 1 : E[ H^H  Xd^H  D  Xd  H]
D = diag(Xp1)'*(Z)*diag(Xp1);

D = roundn(D,-12); % Truncate to 12 digits to make real part symmetric

psi = diag(Xd)*H;

E1 = psi'*D*psi/N   % Works
% proof expressions
% ...
% ...
E11 = trace(D*psi*psi')/N  % Works
E22 = trace(D* trace(psi*psi')/N )/(N*sqrt(N)*2)  % Does not work
% ...
% ...
% end of proof
E2 = trace(D* (psi'*psi)/N )/N % Does not work


%% Arrangement 2 : E[ Xd^H  H^H  Xp^H  U^H  U  Xp  H  Xd]
phi = U*diag(Xp1)*diag(H);
SigmaPhi = phi'*phi;

E3 = Xd'*SigmaPhi*Xd/N  % Works
% proof expressions
% ...
% ...
E33 = trace(SigmaPhi* Xd*Xd' )/N  % Works
E44 = trace(SigmaPhi* trace(Xd*Xd')/N )/(N*sqrt(N))  % Does not work
% ...
% ...
% end of proof
E4 = trace(SigmaPhi* (Xd'*Xd)/N )/N  % Does not work



%% Measure the error from simulation
xp0 =  do_modulate(p, Dp);
xd0 = do_modulate(p, Dd);
InterF = H.*(fft(xd0)./fft(xp0));
interf = ifft(InterF);
interf(numPaths+1:end)=0;
InterFh = fft(interf);

E0 = sum(abs(InterFh).^2)./N % Measured error

E2 = E22;

end

