close all;
addpath('../')
outdir = 'MyOut';
if ~exist(outdir,'dir')
    mkdir(outdir)
end

Ns=41;  % Number of source locations
Nr=11;  % Number of receivers
Not=6;

Xgrid=linspace(-5,5,Ns);%source location
Xr=linspace(-5,5,Nr);%receiver location
% OT=linspace(0,3,Not);%origin time
L=50;

Xrx = linspace(0,5,Nr);
XrL = L - 0.4*Xrx;

%%
t=-1:0.01:3;
nt=length(t);

freq=5;%Hz
w=2*pi*freq;
c=3;%km/s
k=w/c;

%compute arrival times from epicenter as Ta
%Ta=sqrt( Xr.^2 + L^2 )/c;
Ta = sqrt( Xrx.^2 + XrL.^2)/c;

Data=zeros(Nr,nt);
st=0.03;%time scale of gaussian in second


%P=[1 1 1];
%Xs=[0 1 3];
Xs = linspace(0,3,30);
P = ones(size(Xs));
Vr = 1/0.3*ones(size(Xs)).*(Xs<1)+2/0.3*ones(size(Xs)).*(Xs>=1);
OT = cumsum((Xs(2)-Xs(1))./Vr);
%OT=[0 0.3  0.6];
Nsub=length(P);
%phi=[0 0 0];
myphi = load('Phirand');%pi*rand(size(Xs));
phi = myphi.phi;
%phi = 2*pi*rand(size(Xs))-pi;
smoothphi = smooth(phi);
%%
%%
SNR = 20;
%Data=forward(Nr,nt,Nsub, P, Xs, OT, Xr, L, c, Ta, phi, st,w,t);
Data=forward(Nr,nt,Nsub, P, Xs, OT, Xrx, XrL, c, Ta, phi, st,w,t);
DataSmooth=forward(Nr,nt,Nsub, P, Xs, OT, Xrx, XrL, c, Ta, zeros(size(phi)), st,w,t);
Data2=forward(Nr,nt,Nsub, P, Xs, OT, Xrx, XrL, c, Ta, smoothphi, st,w,t);
%Data = DataSmooth;
for i = 1:Nr
    Data(i,:) = awgn(Data(i,:),SNR);
end

%%
% figure(2);clf;
% 
% plot(XrL,Xrx, 'o'); hold on;
% plot(zeros(size(Xs)),Xs,'rp')
% xlim([0 max(XrL)+2])
% ylim([min(Xrx)-2 max(Xrx)+2])
% xlabel('Transverse Distance')
% 
% figure(1);clf;
% set(gcf,'Position',[360 1 754 697])
% subplot(4,1,1)
% plot(Xs,phi/pi); hold on;
% plot(Xs,smoothphi/pi);
% xlabel('Source location')
% title('Random source phase')
% subplot(4,1,2)
% title('Data with Zero Phase')
% hold on;
% for kk=1:Nr
%     plot(t,DataSmooth(kk,:)+2*kk,'k');
% end
% box on;
% ylim([-1 2*Nr+5])
% xlabel('Time (s)')
% subplot(4,1,3)
% title('Data with Random Phase')
% hold on;
% for kk=1:Nr
%     plot(t,Data(kk,:)+2*kk,'k');
% end
% box on;
% ylim([-1 2*Nr+5])
% xlabel('Time (s)')
% 
% subplot(4,1,4)
% title('Data with Smoothed Random Phase')
% hold on;
% for kk=1:Nr
%     plot(t,Data2(kk,:)+2*kk,'k');
% end
% box on;
% ylim([-1 2*Nr+5])
% xlabel('Time (s)')
% return


%% inverse for phi of 2nd and 3rd subevent
nDiv = 1;
dtij = zeros(Ns,Nr);
for si = 1:Ns
      x_xi = Xgrid(si);
      tij = sqrt( ( x_xi-Xrx ).^2 + ( XrL ).^2 )/c;
      dtij(si,:) = Ta-tij;
end

lowF = 1.0;
highF = 7.0;

DataFilt = Data;
dt = t(2)-t(1);
fnyq = 1/dt/2;
[B,A] = butter(4,[lowF highF]./fnyq);
for st=1:Nr
    DataFilt(st,:) = filter(B,A,Data(st,:));
end

nfft = 2^nextpow2(nt);
fspace0 = 1/dt * (0:(nfft/2))/nfft;
nf = length(fspace0);

% Bin the frequencies
df = fspace0(2)-fspace0(1);
fL = 0.0;
fH = 7.0;
ffilt = find(fspace0 >= fL & fspace0 <=fH);
dffilt = 1;
ffilt = ffilt(1:dffilt:end);
df = dffilt*df;

binpop = ceil(0.1/df);
overflow = binpop - mod(length(ffilt),binpop);
if overflow ~= 0
   ffilt = ffilt(1):(ffilt(end)+overflow); 
end
fspace = fspace0(ffilt);
nf = length(fspace); % number of frequencies
nfbin = nf/binpop;

DataSpec = zeros(Nr,nf);
DataSpecraw = zeros(Nr,nfft/2+1);
for i = 1:Nr
    spec = fft(DataFilt(i,:),nfft);
    spec = spec(1:nfft/2+1);
    DataSpecraw(i,:) = spec;
    DataSpec(i,:) = spec(ffilt);
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%             Prepare and Perform Inversion              %
%              for each Discrete Frequency               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% Redefine values based on population of subarrays 
% (currently same as total stations)
np = Nr;%sum(DivPop);         % number of stations in inversion population
ncomb = Ns*nDiv;          % total number of model parameters
pl = sqrt(nDiv)*ones(1,Ns);  

% Sparsity parameter
Orders = [-2;-1;0;1;2];
factors = [1;5];
Lambdas = zeros(length(Orders)*length(factors),1);
for i1 = 1:length(Orders)
    for i2 = 1:length(factors)
        Lambdas((i1-1)*length(factors)+i2) = factors(i2)*10^(Orders(i1));
    end
end
nLam = length(Lambdas);

cvx_solver_settings('cvx_slvitr',2);
%cvx_solver_settings -clear
tic
for fbin = 1:nfbin
    % Output models
    moutTemp = zeros(nLam,ncomb*binpop);
%    mout = zeros(binpop,ncomb);
    mm   = zeros(nLam,binpop,nDiv,Ns);

    % Displacement vector for entire population
    uom = zeros(np*binpop,1);

    % Synthetic spectra 
    syntmp = zeros(np*binpop,nLam);

    % Spectral Power for each source
    specPowerF = zeros(Ns,nDiv,nLam);
    
    findices = ((fbin-1)*binpop+1):(fbin*binpop);
    f0s = fspace(findices); % frequency

    % Fill data vectors for frequency
    u = reshape(DataSpec(:,findices),np*binpop,1);

    % Create kernels for each source location and station
    K1 = zeros(binpop*np,binpop*Ns);
    Kf = zeros(binpop*np,binpop*ncomb);

    for d = 1:nDiv
        % find the station indices within the subarray
        popu = 1:length(Xs);%((sum(DivPop(1:d))+1):(sum(DivPop(1:d+1))));
        for i = 1:ns 
            for fi = 1:binpop
                % Fill kernels
                K1((fi-1)*np+popu,(fi-1)*Ns+i) = (exp(2i*pi*f0s(fi)*(dtij(i,popu)')));
                Kf((fi-1)*np+popu,(fi-1)*nDiv*Ns+(d-1)*Ns+i) = (exp(2i*pi*f0s(fi)*(dtij(i,popu)')));
            end
        end
    end
    
    parfor f = 1:nLam       % parallelized over frequency
        % Perform the Inversion
        lambda = Lambdas(f);                         % Sparsity prior weight
        m = GroupLassoBin(u,Kf,pl,lambda,Ns,ncomb,binpop);

        syntmp(:,f) = Kf*m;
        tmpspecPowerF = zeros(Ns,nDiv);
        moutTemp(f,:) = m
        for fi = 1:binpop
            fsource = ((fi-1)*ncomb+1:fi*ncomb);
            mtmp = m(fsource);

            % Calculate power and synthetics at each frequency from the subevents
            mmtmp = zeros(Ns,nDiv);
            tmpspecPower = zeros(Ns,nDiv);
            for d = 1:nDiv
                popu = ((sum(DivPop(1:d))+1):(sum(DivPop(1:d+1))));
                Ktemp = K1((fi-1)*np+popu,((fi-1)*Ns+1):fi*Ns);
                for s = 1:Ns
                    mmtmp(s,d) = mtmp((d-1)*Ns+s);
                    tmp = Ktemp(:,s)*mmtmp(s,d);
                    tmpspecPower(s,d) =  sum(real(tmp).*real(tmp));
                end
            end
            tmpspecPowerF = tmpspecPowerF + tmpspecPower;
        end
        specPowerF(:,:,f) = tmpspecPowerF;

    %%
    end
toc
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%               Reorganize model array                   %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
    for la = 1:nLam
        for fi = 1:binpop
            for d = 1:nDiv
                rang = ((fi-1)*ncomb+(d-1)*Ns+1: (fi-1)*ncomb+d*Ns);
                mm(la,fi,d,:) = moutTemp(la,rang);
            end
        end
    end

    %% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %                       Save Info                        %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
    info.f0s = f0s;
    info.x_ev = Xs;
    info.t_ev = OT;
    info.m_ev = P;
    info.x = Xgrid;
    info.nx = Ns;
    info.ns = length(Xs);
    info.lowF = lowF;
    info.highF = highF;
    info.binpop = binpop;
    info.fspace = fspace;
    info.t = t;
    save([outdir,sprintf('InversionOutput_%d.mat',fbin)],'spec','DataSpec','syntmp','specPowerF','mm','Lambdas','fspace','info','-v7.3');
end
poolobj = gcp('nocreate');
delete(poolobj);    

%%


function Pred=forward(Nr,nt,Nsub, P, Xs, OT, Xrx, XrL, c, Ta, phi, st,w,t)
Pred=zeros(Nr,nt);
for ii=1:Nsub
    A=P(ii);
    xs=Xs(ii);
    ot=OT(ii);
    T = ot + sqrt( (Xrx-xs).^2 + XrL.^2 )/c - Ta;
    for kk=1:Nr
        Pred(kk,:)=Pred(kk,:)+A*cos(w*(t-T(kk))-phi(ii)).*exp(-(t-T(kk)).^2/2/st^2);
    end
end
end