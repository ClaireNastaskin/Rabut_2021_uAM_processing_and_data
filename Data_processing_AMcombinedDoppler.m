% Processing script for uAM_Doppler data acquired with
% RunSetUp_uAM_Doppler_acquire.m acquisition script on a Verasonics scanner
% and beamformed with RunSetUp_uAM_Doppler_reconAll.m script
%
% This script loads the beamformed data and will process them to give two
% set of data: the AM signal (nonlinear signal only from the medium) and
% the Power Doppler signal. The Power Doppler signal is obtained applying a
% clutter filtering (SVD filtering: Demene et al, 2015). The cutoff value
% of this filter might need to be adjusted depending of the set of data.
%
% Rabut et al., Appl. Phys. Lett. 118, 244102 (2021); https://doi.org/10.1063/5.0050807

clear, clc, format compact, close all

%% Data to process
Path = 'D:\Your_path\YourFileName';
cd(Path)
load([Path '\UF.mat'])

%% Input parameters 

n_eig = floor(UF.numFrames/5) ; % SVD cutoff value: might need to vary depending of the ratio blood/tissue motion

acqLength = UF.NbOfBlocs;
lat_resol = 2*UF.imaging_aperture;
axial_resol = 2*ceil(ceil(UF.Depth(2)/UF.Lambda)-floor(UF.Depth(1)/UF.Lambda));
ensemble_length = UF.numFrames;

%% AM calculation and Doppler filtering


for ii = 1:acqLength   
    
    fid = fopen([sprintf('IQ_AMneg_%.3d.bin', ii)], 'r');
    IQ_neg = fread(fid, 'double');
    fclose(fid);
    IQ_temp_neg = reshape(IQ_neg, [],lat_resol*2, ensemble_length);
    IQ_complex_neg = IQ_temp_neg(:,1:lat_resol,:)+1i*IQ_temp_neg(:,(lat_resol)+1:lat_resol*2,:);
     
    fid = fopen([sprintf('IQ_AMpos_%.3d.bin', ii)], 'r');
    IQ_pos = fread(fid, 'double');
    fclose(fid);
    IQ_temp_pos = reshape(IQ_pos, [],lat_resol*2, ensemble_length);
    IQ_complex_pos = IQ_temp_pos(:,1:lat_resol,:)+1i*IQ_temp_pos(:,(lat_resol)+1:lat_resol*2,:);
    
    clear IQBmode
    IQBmode = squeeze(IQ_complex_pos-IQ_complex_neg); 
    
    % AM_signal
    clear IQ_AM
    IQ_AM= IQ_complex_neg + IQ_complex_pos;
    AM_signal(:,:,ii) = squeeze(mean(abs(IQ_AM),3)) ;
    
    % Doppler signal
    
    IQ_signal = IQBmode ; 
    [nz, nx, nt] = size(IQ_signal);    
    IQ_signal = reshape(IQ_signal, [nz*nx, nt]);
    cov_matrix = IQ_signal'*IQ_signal;
    [Eig_vect, Eig_val]= eig(cov_matrix);
    Eig_vect=fliplr(Eig_vect);
    Eig_val=rot90(Eig_val,2);
    M_A = IQ_signal*Eig_vect;
    skipped_eig_val = 1:n_eig; 
    IQF_tissu = M_A(:,skipped_eig_val)*Eig_vect(:,skipped_eig_val)';
    IQF_tissu = reshape(IQF_tissu, [nz, nx, nt]);
    IQ_signal = reshape(IQ_signal, [nz, nx, nt]);
    IQF_corrected = IQ_signal-IQF_tissu;
    Dop(:,:,ii) = mean(abs(IQF_corrected(:,:,:)).^2,3);
    
    ii
end
return
