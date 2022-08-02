function [mean_Kappa_error,std_Kappa_error,PP_FAs,IT_FAs,Kappa_error_results] = GenDataAndAnalyse(Scheme,ratios,T1s,T2,TRs,Noise_pcs,repeats,Reordering,Phase_Resolution, Partial_Fourier, Matrix_Size,Segment_Factor)
% Function handles iterating through scheme(s) to simulate and hands over to analysing
% data in Analysis_function_core.
%
%       INPUT:
%           Scheme - Name of scheme(s) to simulate
%           ratios - Equal to 0 to produce a 'parameter' map of PP_FA vs IT_FA OR Array of fixed
%                    ratios of PP_FA to IT_FA e.g. [5 10 20 30]
%           T1s - Array of T1 values to simulate (s)
%           T2 - T2 value to simulate (s)
%           TRs - Array of TR values to simulate (s)
%           Noise_pcs - Array of fixed noise levels to simulate
%           repeats - Number of Monte Carlo simulations to repeat and
%           average over
%           Partial_Fourier - fraction of phase encodes acquired
%           Matrix_Size - Size of the scan matrix to be zero-filled up to
%           UseSyntheticData - Are we simulating with syntheticdata? (1)
%           yes or (0) no.
%
%       OUTPUT:
%           mean_Kappa_error - Mean of repeats for the error on the reconstructed
%           pre-pulse flip angle (B1 map)
%           std_Kappa_error - standard deviation of the error on the reconstructed
%           pre-pulse flip angle (B1 map)
%           PP_FAs - Array of pre-pulse flip angles used to simulate
%           IT_FAs - Array of image train flip angles
%           Kappa_error_results - cell containg results from ALL simulated
%           schemes

% James Kent. 2021. Using B. Hargreaves EPG Functions.
Velocity = 0; % Velocity of coherent flow (m/s)
Angle = 0; % Angle of coherent flow (rad)
Diff_co = 0; % Diffusion co-efficient (m^2/s) (isotropic)
Kappa_error_results = cell(2,1);
Train_Size = round(Phase_Resolution*Partial_Fourier*Matrix_Size);

if strcmp(Scheme,'satTFL')
    seq_n = 1;
    [ITs,PP_FAs,IT_FAs,Segment_Sizes] = epg_satTFL(TRs,T1s,T2,ratios,Velocity,Angle,Diff_co,Train_Size,Segment_Factor);
    disp('Simulation Successful')
    [mean_Kappa_error,std_Kappa_error] = Analysis_function_core(seq_n, ITs, Noise_pcs, repeats,Reordering,Phase_Resolution,Partial_Fourier,Matrix_Size,Segment_Sizes);
    
elseif strcmp(Scheme,'ShortTR')
    seq_n = 2;
    [ITs,PP_FAs,IT_FAs,Segment_Sizes] = epg_ShortTR(TRs,T1s,T2,ratios,Velocity,Angle,Diff_co,Train_Size,Segment_Factor);
    disp('Simulation Successful')
    [mean_Kappa_error,std_Kappa_error] = Analysis_function_core(seq_n, ITs, Noise_pcs, repeats,Reordering,Phase_Resolution,Partial_Fourier,Matrix_Size,Segment_Sizes);
elseif strcmp(Scheme,'Sandwich')
    seq_n = 3;

        [ITs,PP_FAs,IT_FAs,Segment_Sizes] = epg_Sandwich(TRs,T1s,T2,ratios,Velocity,Angle,Diff_co,Train_Size,Segment_Factor);
        disp('Simulation Successful')
        [mean_Kappa_error,std_Kappa_error] = Analysis_function_core(seq_n, ITs, Noise_pcs, repeats,Reordering,Phase_Resolution,Partial_Fourier,Matrix_Size,Segment_Sizes);
        
    
else
    error('ABORTED: Scheme not recognised, please input either ''satTFL'', ''ShortTR'', ''Sandwich'' OR ''ALL''.')
end
end

function [mean_Kappa_error,std_Kappa_error] = Analysis_function_core(seq_n, ITs, Noise_pcs, repeats,Reordering,Phase_Resolution, Partial_Fourier, Matrix_Size,Segment_Sizes)
% Function handles the actually analysis i.e. Reordering, zero-filling,
% generating synthetic noise, Fourier transform, calculating alpha and bias
% correction (if switched on) and returns alpha.
%
% Only processes one set of image trains at a time (E.g. IT1 and IT2/IT3
% for a specific scheme and for multiTx it is ran for however modes there are)

% Extract Image Trains
IT1 = ITs{1};
IT2 = ITs{2};

% Re-order Image Trains (zero centred)
% Collected as 1,2,3,...,30,31,32 in the image train.

    Reordered_IT1 = zeros(size(IT1));
    Reordered_IT2 = zeros(size(IT1));
    
    % Determine which line in the readout is the centre (i.e. not the centre of Matrix_Size for partial Fourier)
    Centre_Line = ceil(Phase_Resolution*(Matrix_Size/2 - (Matrix_Size*(1-Partial_Fourier)))); % (zero-indexed)
    if strcmp(Reordering,'CentricOut')
        % We need to determine the ordering of the data, taking into account
        % any partial Fourier or reduced phase resolution.
        A = Centre_Line+1:1:(round(Phase_Resolution*Partial_Fourier*Matrix_Size)-1);
        B = (Centre_Line)-1:-1:0;
        N = min(numel(A),numel(B));
        Reorder = [Centre_Line,reshape([B(1:N);A(1:N)],1,[]),A(N+1:end),B(N+1:end)] +1; % Define Re-ordering (+1 accounts for non-zero indexing of Matlab). Since it is centric we have to interleave the reordering
    elseif strcmp(Reordering,'CentricIn')
        
        %         A = (round(Phase_Resolution*Partial_Fourier*Matrix_Size)-1):-1:Centre_Line+1;
        %         B = 0:+1:(Centre_Line)-1;
        %         N = min(numel(A),numel(B));
        %         Reorder = [reshape([B(1:N);A(1:N)],1,[]),A(N+1:end),B(N+1:end),Centre_Line] +1; % Define Re-ordering (+1 accounts for non-zero indexing of Matlab). Since it is centric we have to interleave the reordering
        %
    elseif strcmp(Reordering,'Linear')
        Reorder = (0:1:round(Phase_Resolution*Partial_Fourier*Matrix_Size)-1) +1; % Define Re-ordering (+1 accounts for non-zero indexing of Matlab)
    else
        disp('Reordering scheme not recognised. Reordering not performed.');
    end
    
    for Seg_n = 1:size(Segment_Sizes,2)
        % Reorder taking into account any segmentation that has occured due
        % to pulse sequence design
        Reordered_Seg(1,sum(Segment_Sizes(1:Seg_n - 1))+1 : sum(Segment_Sizes(1:Seg_n))) = Reorder(Seg_n:size(Segment_Sizes,2):end);
    end
    if size(Reorder,2) ~= size(Reordered_Seg,2)
        error('ERROR: Segment reordering failed.')
    else
        Reorder = Reordered_Seg;
    end
    clearvars Reordered_Seg A B N
    
    for IT_n = 1:Phase_Resolution*Partial_Fourier*Matrix_Size % (Train size)
        
        % Perform the reordering of the data
        Reordered_IT1(:,:,:,:,Reorder(IT_n),:) = IT1(:,:,:,:,IT_n,:);
        Reordered_IT2(:,:,:,:,Reorder(IT_n),:) = IT2(:,:,:,:,IT_n,:);
        
    end
clearvars IT1 IT2 IT3



% Add zero-mean complex Gaussian noise independently to each reordered IT
IT1_Noisy = Reordered_IT1 + bsxfun(@times,randn([size(Reordered_IT1,1:6),size(Noise_pcs,2),repeats]),reshape(Noise_pcs,[1 1 1 1 1 1 size(Noise_pcs,2) 1])) + 1i.*bsxfun(@times,randn([size(Reordered_IT1,1:6),size(Noise_pcs,2),repeats]),reshape(Noise_pcs,[1 1 1 1 1 1 size(Noise_pcs,2) 1]));
IT2_Noisy = Reordered_IT2 + bsxfun(@times,randn([size(Reordered_IT2,1:6),size(Noise_pcs,2),repeats]),reshape(Noise_pcs,[1 1 1 1 1 1 size(Noise_pcs,2) 1])) + 1i.*bsxfun(@times,randn([size(Reordered_IT2,1:6),size(Noise_pcs,2),repeats]),reshape(Noise_pcs,[1 1 1 1 1 1 size(Noise_pcs,2) 1]));
clearvars Reordered_IT1 Reordered_IT2 Reordered_IT3


% Zero-fill Phase-Encode, accounts for partial Fourier by unevenly filling
if Partial_Fourier ~= 1 || Phase_Resolution ~= 1
    Pad_Size = Matrix_Size-size(IT1_Noisy,5);
    PadTop_Size = (Matrix_Size/2)-Centre_Line;
    PadBot_Size = Pad_Size-PadTop_Size;
    PadTop = zeros([size(IT1_Noisy,1:4),PadTop_Size,size(IT1_Noisy,6:8)]);
    PadBot = zeros([size(IT1_Noisy,1:4),PadBot_Size,size(IT1_Noisy,6:8)]);
    IT1_Noisy = cat(5,cat(5,PadTop, IT1_Noisy),PadBot);
    IT2_Noisy = cat(5,cat(5,PadTop, IT2_Noisy),PadBot);
end
clearvars PadTop PadBot

% Fourier transform image train, fft in col and lin dimension:
for f = [5 6]
    IT1_Noisy =  ifftshift(ifft(fftshift(IT1_Noisy,f),[],f),f);
    IT2_Noisy =  ifftshift(ifft(fftshift(IT2_Noisy,f),[],f),f);
end
FT_IT1 = IT1_Noisy;
FT_IT2 = IT2_Noisy;
clearvars IT1_Noisy IT2_Noisy

FT_IT1 = permute(FT_IT1,[1 2 3 4 7 8 5 6]); % Remove dimension of length 1
FT_IT2 = permute(FT_IT2,[1 2 3 4 7 8 5 6]);
% Take centre of PSF
    Max_Val_IT1 = FT_IT1(:,:,:,:,:,:,round(Matrix_Size/2 +1),:); % Store maximum value of IT1
    Max_Val_IT2 = FT_IT2(:,:,:,:,:,:,round(Matrix_Size/2 +1),:); % Store maximum value of IT1
clearvars  FT_IT3 FT_IT1 FT_IT2

Kappa = acos((Max_Val_IT2./Max_Val_IT1)); % Eq. 2 from Chung, S. et al. MRM. 2010. 64(2):439-446.

Kappa = (180/pi)*real(Kappa);

mean_Kappa_error = mean(Kappa,6);
std_Kappa_error = std(Kappa,[],6);

end



