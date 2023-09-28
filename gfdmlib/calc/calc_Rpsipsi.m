function [R_PsiPsi] = calc_Rpsipsi(method, nT, nR, conf, Dd, Dd1_indices)
%
% method:
%       1: Estimate the channel using only one pilot subsymbol
%       2: Estimate the channel through Harmonic mean of the pilot subsymbols
%       
%
%
% X:    FFT of transmit signal
% Y:    FFT of received signal
% Dp:   K*M Matrix of Pilot signals with zero on data positions
% numPaths: Number of paths in FSC situation
%
% Author: Shahab Ehsanfar


J = conf.J;
P = conf.P;
numPaths = conf.numPaths;
noblk = 2; %conf.noblk;
N = conf.N;
Ncp = conf.Ncp;
F_L = conf.F_L;
A = conf.A;
F_cpblk= conf.F_cpblk;
FcpIcpA= conf.FcpIcpA;

switch method
    case 'R_PsiPsi_Z'
        % initializations
        
        FA = conf.FA;
        F = conf.F;
        M = conf.M;
        K = conf.K;
        F_L = conf.F_L_N;
        
        F_M = dftmtx(M)/sqrt(M);
        %m = 1;
        %Z_m = kron(F_M(m,:),eye(K));
        Z = kron(F_M,eye(K));
        
        % Initialize matrices
        sigma2d = false(N,nT);
        R_hh_i = zeros(numPaths,numPaths,nT,nR);
        Upsilon = zeros(N,N,nT,nR);              
        R_XdXd = zeros(N,N,nT);
        innerR_PsiPsi_i = zeros(N,N,nT,nR);   
        R_PsiPsi_i = zeros(N,N,nT,nR);       
        R_PsiPsi_iR = zeros(N,N,nR);
                       
        % interference calculation
        for iR = 1:nR
            for iT = 1:nT
                R_hh_i(:,:,iT,iR) = diag(P(:,iT,iR));
                Upsilon(:,:,iT,iR) = N*F_L*R_hh_i(:,:,iT,iR)*F_L';
                
                Dd_i = Dd(:,:,iT);
                
                sigma2d(:,iT) = [(abs( Dd_i(:) ) > 0)]; %(abs( repmat(Dd_i(:),noblk,1) ) > 0) + 0;
             %   sigma2d_blk(:,iT) = (abs( Dd_blk0(:) ) > 0) + 0;
                
             %   IcpA = blkdiag(A(end-Ncp+1:end,:), I_cp*A); %kron(eye(noblk),I_cp*A);
             
%                 sigma2d(:,iT) = (abs( Dd_i(:) ) > 0) + 0;                
%                 R_XdXd(:,:,iT) = F_cpblk*IcpA*diag([ones(N,1) ;sigma2d(:,iT)])*IcpA'*F_cpblk';   
                
%                 FcpIcpAsigmad = FcpIcpA .* repmat(sqrt(sigma2d(:,iT)+0).',[NcpNcpN 1]);             
%                 test2 = FcpIcpAsigmad*FcpIcpAsigmad';
               
                R_XdXd(:,:,iT) = FA(:,sigma2d(:,iT))*FA(:,sigma2d(:,iT))';
            
                
                innerR_PsiPsi_i(:,:,iT,iR) = Upsilon(:,:,iT,iR) .* R_XdXd(:,:,iT);
                
                ZF = Z*F';
                R_PsiPsi_i(:,:,iT,iR) = ZF*innerR_PsiPsi_i(:,:,iT,iR)*ZF';
                
                R_PsiPsi_iR(:,:,iR) = R_PsiPsi_iR(:,:,iR) + R_PsiPsi_i(:,:,iT,iR);
            end
        end
        
        % Generate the large block diagonal R_PsiPsi from the individual
        % blocks
        R_PsiPsi_str = char(zeros(nR,17 + numel(num2str(nR)) ));
        for iR = 1:nR
            R_PsiPsi_str(iR,1:17+numel(num2str(iR))) = sprintf('R_PsiPsi_iR(:,:,%d)',iR);
            
            if iR == 1
                R_PsiPsi_blks = R_PsiPsi_str(iR,1:17+numel(num2str(iR)));
            else
                R_PsiPsi_blks = [R_PsiPsi_blks ',' R_PsiPsi_str(iR,1:17+numel(num2str(iR)))];
            end
        end
        R_PsiPsi = eval(['blkdiag(' R_PsiPsi_blks ');']);

%%
    case 'R_PsiPsi_CPZm'
        % initializations
        
        NNcp = N+Ncp;
        NcpNcpN = Ncp+NNcp;
        M = conf.M;
        K = conf.K;
        F_M = dftmtx(M)/sqrt(M);
        
        R_PsiPsi = zeros(K*nR,K*nR,M);
        for m = 1:M
            Z_m = kron(F_M(m,:),eye(K));%;repmat(diag(F_M(m,:)),[1 48]);%
            
            sigma2d = false(noblk*N,nT);
            %sigma2d = zeros(N,nT);
            %         sigma2d_blk = zeros(1*N,nT);
            R_hh_i = zeros(numPaths,numPaths,nT,nR);
            Upsilon = zeros(NcpNcpN,NcpNcpN,nT,nR);
            
            R_XdXd = zeros(NcpNcpN,NcpNcpN,nT);
            innerR_PsiPsi_i = zeros(NcpNcpN,NcpNcpN,nT,nR);
            
            
            R_PsiPsi_i = zeros(K,K,nT,nR);
            
            
            R_PsiPsi_iR = zeros(K,K,nR);
            
            
            % interference calculation
            for iR = 1:nR
                for iT = 1:nT
                    R_hh_i(:,:,iT,iR) = diag(P(:,iT,iR));
                    Upsilon(:,:,iT,iR) = NcpNcpN*F_L*R_hh_i(:,:,iT,iR)*F_L';
                    
                    Dd_i = Dd(:,:,iT);
                    
                    sigma2d(:,iT) = [Dd1_indices; (abs( Dd_i(:) ) > 0)]; %(abs( repmat(Dd_i(:),noblk,1) ) > 0) + 0;
                    %   sigma2d_blk(:,iT) = (abs( Dd_blk0(:) ) > 0) + 0;
                    
                    %   IcpA = blkdiag(A(end-Ncp+1:end,:), I_cp*A); %kron(eye(noblk),I_cp*A);
                    
                    %                 sigma2d(:,iT) = (abs( Dd_i(:) ) > 0) + 0;
                    %                 R_XdXd(:,:,iT) = F_cpblk*IcpA*diag([ones(N,1) ;sigma2d(:,iT)])*IcpA'*F_cpblk';
                    
                    %                 FcpIcpAsigmad = FcpIcpA .* repmat(sqrt(sigma2d(:,iT)+0).',[NcpNcpN 1]);
                    %                 test2 = FcpIcpAsigmad*FcpIcpAsigmad';
                    
                    R_XdXd(:,:,iT) = FcpIcpA(:,sigma2d(:,iT))*FcpIcpA(:,sigma2d(:,iT))';
                    
                    
                    innerR_PsiPsi_i(:,:,iT,iR) = Upsilon(:,:,iT,iR) .* R_XdXd(:,:,iT);
                    
                    ZJF1 = Z_m*conf.JF2;
                    ZJF2 = ZJF1; %kron(F_M(1,:),eye(K))*conf.JF2;
                    R_PsiPsi_i(:,:,iT,iR) = ZJF1*innerR_PsiPsi_i(:,:,iT,iR)*ZJF2';
                    
                    R_PsiPsi_iR(:,:,iR) = R_PsiPsi_iR(:,:,iR) + R_PsiPsi_i(:,:,iT,iR);
                end
            end
            
            % Generate the large block diagonal R_PsiPsi from the individual
            % blocks
            R_PsiPsi_str = char(zeros(nR,17 + numel(num2str(nR)) ));
            for iR = 1:nR
                R_PsiPsi_str(iR,1:17+numel(num2str(iR))) = sprintf('R_PsiPsi_iR(:,:,%d)',iR);
                
                if iR == 1
                    R_PsiPsi_blks = R_PsiPsi_str(iR,1:17+numel(num2str(iR)));
                else
                    R_PsiPsi_blks = [R_PsiPsi_blks ',' R_PsiPsi_str(iR,1:17+numel(num2str(iR)))];
                end
            end
            R_PsiPsi(:,:,m) = eval(['blkdiag(' R_PsiPsi_blks ');']);
        end
    
    
%%    
        case 'R_PsiPsi_CP'
        % initializations
        
        NNcp = N+Ncp;
        NcpNcpN = Ncp+NNcp;
        
        sigma2d = false(noblk*N,nT);
        %sigma2d = zeros(N,nT);
%         sigma2d_blk = zeros(1*N,nT);
        R_hh_i = zeros(numPaths,numPaths,nT,nR);
        Upsilon = zeros(NcpNcpN,NcpNcpN,nT,nR);        
        
        I_cp = eye(N);        
        I_cp = [I_cp(end-Ncp+1:end,:); I_cp];
        
        R_XdXd = zeros(NcpNcpN,NcpNcpN,nT);
        innerR_PsiPsi_i = zeros(NcpNcpN,NcpNcpN,nT,nR); 
        
        
        R_PsiPsi_i = zeros(NNcp,NNcp,nT,nR);
        
        
        R_PsiPsi_iR = zeros(NNcp,NNcp,nR);
                       
        
        % interference calculation
        for iR = 1:nR
            for iT = 1:nT
                R_hh_i(:,:,iT,iR) = diag(P(:,iT,iR));
                Upsilon(:,:,iT,iR) = NcpNcpN*F_L*R_hh_i(:,:,iT,iR)*F_L';
                
                Dd_i = Dd(:,:,iT);
                
                sigma2d(:,iT) = [Dd1_indices; (abs( Dd_i(:) ) > 0)]; %(abs( repmat(Dd_i(:),noblk,1) ) > 0) + 0;
             %   sigma2d_blk(:,iT) = (abs( Dd_blk0(:) ) > 0) + 0;
                
             %   IcpA = blkdiag(A(end-Ncp+1:end,:), I_cp*A); %kron(eye(noblk),I_cp*A);
             
%                 sigma2d(:,iT) = (abs( Dd_i(:) ) > 0) + 0;                
%                 R_XdXd(:,:,iT) = F_cpblk*IcpA*diag([ones(N,1) ;sigma2d(:,iT)])*IcpA'*F_cpblk';   
                
%                 FcpIcpAsigmad = FcpIcpA .* repmat(sqrt(sigma2d(:,iT)+0).',[NcpNcpN 1]);             
%                 test2 = FcpIcpAsigmad*FcpIcpAsigmad';
               
                R_XdXd(:,:,iT) = FcpIcpA(:,sigma2d(:,iT))*FcpIcpA(:,sigma2d(:,iT))';
            
                
                innerR_PsiPsi_i(:,:,iT,iR) = Upsilon(:,:,iT,iR) .* R_XdXd(:,:,iT);
                
                JF = J*F_cpblk';
                R_PsiPsi_i(:,:,iT,iR) = JF*innerR_PsiPsi_i(:,:,iT,iR)*JF';
                
                R_PsiPsi_iR(:,:,iR) = R_PsiPsi_iR(:,:,iR) + R_PsiPsi_i(:,:,iT,iR);
            end
        end
        
        % Generate the large block diagonal R_PsiPsi from the individual
        % blocks
        R_PsiPsi_str = char(zeros(nR,17 + numel(num2str(nR)) ));
        for iR = 1:nR
            R_PsiPsi_str(iR,1:17+numel(num2str(iR))) = sprintf('R_PsiPsi_iR(:,:,%d)',iR);
            
            if iR == 1
                R_PsiPsi_blks = R_PsiPsi_str(iR,1:17+numel(num2str(iR)));
            else
                R_PsiPsi_blks = [R_PsiPsi_blks ',' R_PsiPsi_str(iR,1:17+numel(num2str(iR)))];
            end
        end
        R_PsiPsi = eval(['blkdiag(' R_PsiPsi_blks ');']); 
 %%
        case 'R_PsiPsi_CPZ'
        % initializations
        Z = conf.Z;
        
        NNcp = N+Ncp;
        NcpNcpN = Ncp+NNcp;
        
        sigma2d = false(noblk*N,nT);
        %sigma2d = zeros(N,nT);
%         sigma2d_blk = zeros(1*N,nT);
        R_hh_i = zeros(numPaths,numPaths,nT,nR);
        Upsilon = zeros(NcpNcpN,NcpNcpN,nT,nR);        
        
%         I_cp = eye(N);        
%         I_cp = [I_cp(end-Ncp+1:end,:); I_cp];
        
        R_XdXd = zeros(NcpNcpN,NcpNcpN,nT);
        innerR_PsiPsi_i = zeros(NcpNcpN,NcpNcpN,nT,nR); 
        
        
        R_PsiPsi_i = zeros(N,N,nT,nR);
        
        
        R_PsiPsi_iR = zeros(N,N,nR);
                       
        
        % interference calculation
        for iR = 1:nR
            for iT = 1:nT
                R_hh_i(:,:,iT,iR) = diag(P(:,iT,iR));
                Upsilon(:,:,iT,iR) = NcpNcpN*F_L*R_hh_i(:,:,iT,iR)*F_L';
                
                Dd_i = Dd(:,:,iT);
                
                sigma2d(:,iT) = [Dd1_indices; (abs( Dd_i(:) ) > 0)]; %(abs( repmat(Dd_i(:),noblk,1) ) > 0) + 0;
             %   sigma2d_blk(:,iT) = (abs( Dd_blk0(:) ) > 0) + 0;
                
             %   IcpA = blkdiag(A(end-Ncp+1:end,:), I_cp*A); %kron(eye(noblk),I_cp*A);
             
%                 sigma2d(:,iT) = (abs( Dd_i(:) ) > 0) + 0;                
%                 R_XdXd(:,:,iT) = F_cpblk*IcpA*diag([ones(N,1) ;sigma2d(:,iT)])*IcpA'*F_cpblk';   
                
%                 FcpIcpAsigmad = FcpIcpA .* repmat(sqrt(sigma2d(:,iT)+0).',[NcpNcpN 1]);             
%                 test2 = FcpIcpAsigmad*FcpIcpAsigmad';
               
                R_XdXd(:,:,iT) = FcpIcpA(:,sigma2d(:,iT))*FcpIcpA(:,sigma2d(:,iT))';
            
                
                innerR_PsiPsi_i(:,:,iT,iR) = Upsilon(:,:,iT,iR) .* R_XdXd(:,:,iT);
                
                ZJF = Z*conf.JF2;
                R_PsiPsi_i(:,:,iT,iR) = ZJF*innerR_PsiPsi_i(:,:,iT,iR)*ZJF';
                
                R_PsiPsi_iR(:,:,iR) = R_PsiPsi_iR(:,:,iR) + R_PsiPsi_i(:,:,iT,iR);
            end
        end
        
        % Generate the large block diagonal R_PsiPsi from the individual
        % blocks
        R_PsiPsi_str = char(zeros(nR,17 + numel(num2str(nR)) ));
        for iR = 1:nR
            R_PsiPsi_str(iR,1:17+numel(num2str(iR))) = sprintf('R_PsiPsi_iR(:,:,%d)',iR);
            
            if iR == 1
                R_PsiPsi_blks = R_PsiPsi_str(iR,1:17+numel(num2str(iR)));
            else
                R_PsiPsi_blks = [R_PsiPsi_blks ',' R_PsiPsi_str(iR,1:17+numel(num2str(iR)))];
            end
        end
        R_PsiPsi = eval(['blkdiag(' R_PsiPsi_blks ');']);
        
 %%       
        case 'R_PsiPsi_CP_Dhat'
        % initializations
        
        NNcp = N+Ncp;
        NcpNcpN = Ncp+NNcp;
        
        sigma2d = zeros(noblk*N,nT);
%         sigma2d_blk = zeros(1*N,nT);
        R_hh_i = zeros(numPaths,numPaths,nT,nR);
        Upsilon = zeros(NcpNcpN,NcpNcpN,nT,nR);        
        
%         I_cp = eye(N);        
%         I_cp = [I_cp(end-Ncp+1:end,:); I_cp];
        
        R_XdXd = zeros(NcpNcpN,NcpNcpN,nT);
        innerR_PsiPsi_i = zeros(NcpNcpN,NcpNcpN,nT,nR); 
        
        
        R_PsiPsi_i = zeros(NNcp,NNcp,nT,nR);
        
        
        R_PsiPsi_iR = zeros(NNcp,NNcp,nR);
        
        
        
        
        % interference calculation
        for iR = 1:nR
            if iR == 1
                for iT = 1:nT
                    R_hh_i(:,:,iT,iR) = diag(P(:,iT,iR));
                    Upsilon(:,:,iT,iR) = NcpNcpN*F_L*R_hh_i(:,:,iT,iR)*F_L';
                    
                    Dd_i = Dd(:,:,iT);
                    
                    %sigma2d(:,iT) = [Dd1_indices; (Dd_i(:))]; %(abs( repmat(Dd_i(:),noblk,1) ) > 0) + 0;
                    sigma2d(:,iT) = [Dd1_indices; Dd_i(:)]; %(abs( repmat(Dd_i(:),noblk,1) ) > 0) + 0;
                    
                    %   sigma2d_blk(:,iT) = (abs( Dd_blk0(:) ) > 0) + 0;
                    
                    %   IcpA = blkdiag(A(end-Ncp+1:end,:), I_cp*A); %kron(eye(noblk),I_cp*A);
                    
                    %    R_XdXd(:,:,iT) = FcpIcpA*diag([ones(N,1) ;sigma2d(:,iT)])*FcpIcpA';
                    %    R_XdXd(:,:,iT) = F_cpblk*IcpA*diag([sigma2d_blk(:,iT) ;sigma2d(:,iT)])*IcpA'*F_cpblk';
                    
                    FcpIcpAsigmad = FcpIcpA .* repmat(sqrt(sigma2d(:,iT)).',[NcpNcpN 1]);
                    
                    R_XdXd(:,:,iT) = FcpIcpAsigmad*FcpIcpAsigmad';
                    
                    innerR_PsiPsi_i(:,:,iT,iR) = Upsilon(:,:,iT,iR) .* R_XdXd(:,:,iT);
                    
                    JF = J*F_cpblk';
                    R_PsiPsi_i(:,:,iT,iR) = JF*innerR_PsiPsi_i(:,:,iT,iR)*JF';
                    
                    R_PsiPsi_iR(:,:,iR) = R_PsiPsi_iR(:,:,iR) + R_PsiPsi_i(:,:,iT,iR);
                end
            else
                R_PsiPsi_iR(:,:,iR) = R_PsiPsi_iR(:,:,iR-1);
            end
        end
        
        % Generate the large block diagonal R_PsiPsi from the individual
        % blocks
        R_PsiPsi_str = char(zeros(nR,17 + numel(num2str(nR)) ));
        for iR = 1:nR
            R_PsiPsi_str(iR,1:17+numel(num2str(iR))) = sprintf('R_PsiPsi_iR(:,:,%d)',iR);
            
            if iR == 1
                R_PsiPsi_blks = R_PsiPsi_str(iR,1:17+numel(num2str(iR)));
            else
                R_PsiPsi_blks = [R_PsiPsi_blks ',' R_PsiPsi_str(iR,1:17+numel(num2str(iR)))];
            end
        end
        R_PsiPsi = eval(['blkdiag(' R_PsiPsi_blks ');']); 
 %%                       
end

