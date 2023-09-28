function [Dout, Dp] = do_pilot_symbols_time(p,Din,Delta_k,seq)
% Inserts pilot symbols into the first subsymbol of all subcarriers
% The pilot sequence generation is based on QPSK modulation accourding to
% 3GPP TS 36.211 Physical layer procedures, Release 8.9
%
% Din: D matrix after S/P
%
% Dout: Data with Pilot subsymbols
% Dp: Pilot subsymbols with zero on data positions
%
% Author: Shahab Ehsanfar

if nargin < 4
   seq = 'ZCH';
   if nargin < 3
       Delta_k = 1;
   end
end

if isfield(p,'mode') 
    if strcmp(p.mode,'precode')
        msetp = 1;
    elseif p.Mp == 2 
        msetp = [1 p.M];
    end
elseif p.Mp == 2 
    msetp = [1 p.M];
else
    msetp = 1;
end

% % Construct a PN object
% h = commsrc.pn('Shift', 0);
% 
% % Set the output for p.M PN bits
% set(h, 'NumBitsOut', p.K);

if strcmp(seq,'Rand')
    ZC1 = exp(2i*pi*randn(1,floor(p.Kon/Delta_k)+1)); %lteZadoffChuSeq(1, p.K+1);      %ones(1,p.K);%
elseif strcmp(seq,'Orthog')
    Seq = lteZadoffChuSeq(1, floor(p.K/Delta_k)+1); %exp(2i*pi*randn(1,floor(p.K/Delta_k)));
    Seq = Seq(1:floor(p.K/Delta_k)).';
    [U,~,~] = svd(Seq'*Seq);   
elseif strcmp(seq,'Ones')
    ZC1 = ones(1,floor(p.Kon/Delta_k)+1);
else
    if mod(floor(p.Kon/Delta_k),2) == 0
        ZC1 = lteZadoffChuSeq(1, floor(p.Kon/Delta_k)+1);
        %ZC1 = circshift(lteZadoffChuSeq(1, floor(p.K/Delta_k)+1),floor(64*rand(1)));
    else
        ZC1 = lteZadoffChuSeq(1, floor(p.Kon/Delta_k)+2);
        %ZC1 = circshift(lteZadoffChuSeq(1, floor(p.K/Delta_k)+2),floor(64*rand(1)));
    end
end

kset = get_kset(p)+1;
Dd = Din; 
for p_index = msetp        
%     ZCsize = p.K+1; ZCseq = exp(-1j*pi*13*(0:ZCsize-1).*((0:ZCsize-1)+1)/ZCsize).';
%     if p_index == msetp(1)
%         Din(1:Delta_k:end,p_index) = real(ZCseq(1:p.K)); 
%     else
%         Din(1:Delta_k:end,p_index) = imag(ZCseq(1:p.K)); %real(circshift(ZCseq(1:p.K),1,1)); 
%     end
    if strcmp(seq,'Orthog')
        Din(1:Delta_k:end,p_index) = U(1:floor(p.K/Delta_k),p_index+4).*sqrt(p.K);       
    else
       % Din(kset(1:Delta_k:end),p_index) = ones(size(ZC1(1:floor(p.Kon/Delta_k))));
        Din(kset(1:Delta_k:end),p_index) = ZC1(1:ceil(p.Kon/Delta_k));
        %Din(kset(1:Delta_k:end),p_index) = exp(2i*pi*randn(1,floor(p.Kon/Delta_k))); 
    end
%    Din(1:Delta_k:end,p_index) = exp(2i*pi*randn(1,floor(p.K/Delta_k))); %
    Dd(1:Delta_k:end,p_index) = 0;
end

% ZC1 = lteZadoffChuSeq(1, p.K+1);     
% ZC2 = circshift(lteZadoffChuSeq(1, p.K+1),p.K/2,1);
% Din(:,1) = ZC1(1:p.K);    
% Din(:,p.M) = ZC2(1:p.K);    





% if nargin == 4
%     Dp = zeros(p.K,p.M);
%     Dp(:,msetp) = Din(:,msetp);
%     
%     Dd = Din;
%     Dd(:,msetp) = 0;
%     Dout = Dd;
%     
%     %for k = 2:2:p.K
%         Dout(2:2:p.K,:) = circshift(Dout(2:2:p.K,:),-1,2);
%         Dp(2:2:p.K,:) = circshift(Dp(2:2:p.K,:),-1,2);
%     %end   
% else
    Dout = Dd;
    Dp = zeros(p.K,p.M);
    Dp(kset(1:Delta_k:end),msetp) = Din(kset(1:Delta_k:end),msetp);
% end


