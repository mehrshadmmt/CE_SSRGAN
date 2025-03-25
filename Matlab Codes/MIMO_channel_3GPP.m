%%  3GPP TR 38.901 release 17
% generate random pathdelays, pathgains, AoA, AoD, ZoA, ZoD
clc
clear
close all

Nc = 1;

Nt = 16;
M = 144;
IRS = 144;

samplingrate=1e8;
num_fre = 1;
num_sta = 50;
num_ffading = 20000;

load ch_sta_mtx2 % channel statistics for training channel data
tic
for l=1:1
hh=zeros(Nt,M,num_sta*num_ffading);
for m=1:num_sta
fprintf('m=%d\n',m);
sta_mtx=channel_statistic(:,:,m);
LOSangle=60;
philosAoA=LOSangle;
philosAoD=LOSangle;
thetalosZoA=LOSangle;
thetalosZoD=LOSangle;
cdl = nrCDLChannel;
cdl.CarrierFrequency=28e9;
cdl.TransmitAntennaArray.Size = [Nt 1 1 1 1];
cdl.TransmitAntennaArray.ElementSpacing = [0.5 0.5 1.0 1.0];
cdl.ReceiveAntennaArray.Size = [M 1 1 1 1];
cdl.ReceiveAntennaArray.ElementSpacing = [0.5 0.5 1.0 1.0];
fc=cdl.CarrierFrequency/1e9;
cdl.MaximumDopplerShift = 0;
cdl.ChannelFiltering=false;    
cdl.DelayProfile='Custom'; 
% cdl.DelayProfile='CDL-D'; 
cdl.PathDelays=sta_mtx(1,:);
cdl.AveragePathGains=sta_mtx(2,:);
cdl.AnglesAoA=sta_mtx(3,:);
cdl.AnglesAoD=sta_mtx(4,:);
cdl.AnglesZoA=sta_mtx(5,:);
cdl.AnglesZoD=sta_mtx(6,:);
for n=1:num_ffading
    % if mod( n , 1000 ) == 0
    %     fprintf('n=%d\n',n);
    % end
    cdl.Seed = (l-1)*num_ffading+n;
    [pathgains,sampletimes]=step(cdl);
    pathgains=sqrt(M)*pathgains;
    for nt=1:Nt
        for nr=1:M
            pathpower=pathgains(:,:,nt,nr);
            h=zeros(1,1024);
            I=floor(cdl.PathDelays*samplingrate)+1;
            I_uniq=unique(I);
            Power_sum=[];
            for i=1:length(I_uniq)
                Power_sum(i)=sum(pathpower(I==I_uniq(i)));
            end
            h(I_uniq)=Power_sum;
            h_ntnr=h(1:Nc);
            h_fre=fft(h_ntnr);            
            fre=32;
            pathgains_fre(:,nt,nr)=h_fre; %frequency domain channel matrix
        end
    end
    for j=1:num_fre
        pathgains_fre2(:,:,j)=pathgains_fre(j,:,:);
    end
   hh(:,:,(m-1)*num_ffading+n)=pathgains_fre2;
    release(cdl);
end
end
save(['hh',num2str(IRS)],'hh','-v7.3')

end
toc
