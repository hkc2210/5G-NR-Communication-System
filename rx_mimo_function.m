function [BER_QAM, BER_LDPC, SNR1_Frame, SNR2_Frame, Success_Flag, Time_Frame, ...
          Throughput_QAM, Throughput_LDPC, Capacity_QAM, Capacity_LDPC, ...
          Img_QAM, Img_LDPC, Time_QAM] = ...
    rx_mimo_function(RxStream1, RxStream2, W1, W2, Ref_PSS_Time, ...
                     data_length, img_size, Mat_Interp_Time, maxNumIter, D, CodeRate)

BER_QAM = 0; BER_LDPC = 0; SNR_Frame = 0; Success_Flag = false;
Time_Frame = 0; Throughput_QAM = 0; Throughput_LDPC = 0;
Capacity_QAM = 0; Capacity_LDPC = 0; Img_QAM = []; Img_LDPC = [];
Time_QAM = 0; Time_LDPC_only = 0;

%% Parameters
SyncAdvance = 15;
M = 16;
Nfft = 2048;
Total_Samples = 1228800;
val_DMRS = single(1);

idx_PSS_SC  = [961:1024, 1026:1088];
idx_PBCH_SC = [905:1024, 1026:1145];
idx_SSS_SC  = [961:1024, 1026:1088];

idx_PSS_Sym  = 5;
idx_PBCH_Sym = [6, 8];
idx_SSS_Sym  = 7;
idx_SSB_Sym  = 5:8;

idx_DMRS_SC  = [204:2:1024, 1027:2:1847];
idx_DMRS_Sym = 3:14:560;
idx_Data_SC  = [203:1024, 1026:1847];

SCIdx_DMRS_All = 2:2:1644;
SymIdx_DMRS    = idx_DMRS_Sym;
alt822 = ones(822, 1, 'single');
alt822(2:2:end) = single(-1);

RxStream1 = single(RxStream1(1:1228800));
RxStream2 = single(RxStream2(1:1228800));
Ref_PSS_Time = single(Ref_PSS_Time(:));

persistent Seq_DMRS invPerm NUM_message NUM_VN NUM_CN h_LDPC_rows h_LDPC_cols
persistent Idx_Mat_p SamplesBeforePSS_p idx_lin_p
persistent CodeRate_p

if isempty(Seq_DMRS) || isempty(invPerm)
    Sscr = load('Scrambled.mat');
    Phase_DMRS = reshape(Sscr.Phase_DMRS, length(idx_DMRS_SC), length(idx_DMRS_Sym));
    Seq_DMRS = single(val_DMRS .* exp(1j * Phase_DMRS));

    Sper = load('Permutation.mat');
    Permutation = Sper.Permutation(:);
    invPerm = zeros(length(Permutation), 1);
    invPerm(Permutation) = (1:length(Permutation)).';
end

if isempty(CodeRate_p) || (CodeRate_p ~= CodeRate) % 假如 code rate 是空的或更改 code rate
    h_LDPC_rows = [];
    h_LDPC_cols = [];
    NUM_CN  = []; % the number of rows of H from LDPC
    NUM_VN  = []; % the number of columns of H from LDPC
    NUM_message  = []; % the number of message bits per block
    CodeRate_p = CodeRate;
end

if isempty(h_LDPC_rows)
    switch CodeRate
        case 0.50
            Sld = load('LDPC_11nD2_1296b_R12.mat');
        case 0.67
            Sld = load('LDPC_11nD2_1296b_R23.mat');
        case 0.75
            Sld = load('LDPC_11nD2_1296b_R34.mat');
        case 0.83
            Sld = load('LDPC_11nD2_1296b_R56.mat');
    end

    LDPC = Sld.LDPC;
    Hs = sparse(LDPC.H.x ~= 0);
    [rr, cc] = find(Hs);
    ord = sortrows([rr, cc], [1 2]);
    rr = ord(:,1); cc = ord(:,2);

    h_LDPC_rows = int32(rr - 1);
    h_LDPC_cols = int32(cc - 1);

    NUM_CN = size(Hs, 1);
    NUM_VN = size(Hs, 2);
    NUM_message = NUM_VN - NUM_CN;
end

% ===== (1) Mask + idx_lin =====
if isempty(idx_lin_p)
    Mask_PSS  = false(Nfft, 560);
    Mask_DMRS = false(Nfft, 560);
    Mask_Act  = false(Nfft, 560);

    Mask_PSS(idx_PSS_SC, idx_PSS_Sym) = true;
    Mask_DMRS(idx_DMRS_SC, idx_DMRS_Sym) = true;
    Mask_Act(idx_Data_SC, :) = true;

    Mask_Data = Mask_Act & ~Mask_DMRS & ~Mask_PSS;
    Mask_Data(:, idx_SSB_Sym) = false;
    Mask_Data(idx_PBCH_SC, idx_PBCH_Sym) = true;
    Mask_Data(idx_SSS_SC, idx_SSS_Sym) = true;

    Mask_Data_Act = Mask_Data(idx_Data_SC, :);  % 1644 x 560
    idx_lin_p = find(Mask_Data_Act);            % column-major linear index
end

% ===== (2) OFDM 切窗索引 Idx_Mat + SamplesBeforePSS =====
if isempty(Idx_Mat_p) || isempty(SamplesBeforePSS_p)
    Slot_Pattern = [208; repmat(144, 27, 1)];
    CP_Len_All = repmat(Slot_Pattern, 20, 1);                 % 560x1
    Sym_Start_Pos = cumsum(CP_Len_All + Nfft) - Nfft + 1;     % 560x1

    Idx_Base = (0:Nfft-1).';
    Idx_Mat_p = Idx_Base + Sym_Start_Pos.';                   % 2048 x 560

    symLen = (Nfft + CP_Len_All.');                           % 1 x 560
    SamplesBeforePSS_p = sum(symLen(1:idx_PSS_Sym-1)) + CP_Len_All(idx_PSS_Sym);
end

idx_lin = idx_lin_p;
Idx_Mat = Idx_Mat_p;
SamplesBeforePSS = SamplesBeforePSS_p;

t0 = tic;

%% Sync (PSS) - two independent peaks + same-frame reconcile + simple align

RxLen = min(length(RxStream1), length(RxStream2));
Rx1   = single(RxStream1(1:RxLen));
Rx2   = single(RxStream2(1:RxLen));
Ref   = single(Ref_PSS_Time(:));

% --- 1) independent peak search ---
[Corr1, Lags] = xcorr(Rx1, Ref);
[Corr2, ~]    = xcorr(Rx2, Ref);

[~, i1] = max(abs(Corr1).^2);
[~, i2] = max(abs(Corr2).^2);

PeakLag1 = 1 + Lags(i1) - 1;   % "Ref start index" in Rx1 (1-based)
PeakLag2 = 1 + Lags(i2) - 1;   % "Ref start index" in Rx2

FrameStart1 = PeakLag1 - SamplesBeforePSS + 1 - SyncAdvance;
FrameStart2 = PeakLag2 - SamplesBeforePSS + 1 - SyncAdvance;

% --- 2) force FrameStart2 to be the SAME frame index as FrameStart1 ---
FrameStart2 = FrameStart2 + round(double(FrameStart1 - FrameStart2) / double(Total_Samples)) * Total_Samples;

% --- 3) clamp to valid range (so we can always extract a full frame) ---
FrameStart1 = max(1, min(FrameStart1, RxLen - Total_Samples + 1));
FrameStart2 = max(1, min(FrameStart2, RxLen - Total_Samples + 1));

EndIdx1 = FrameStart1 + Total_Samples - 1;
EndIdx2 = FrameStart2 + Total_Samples - 1;

% --- 4) extract independently ---
Rx1_Synced = Rx1(FrameStart1:EndIdx1);
Rx2_Synced = Rx2(FrameStart2:EndIdx2);

% --- 5) align both frames to Rx1 time origin using integer shift ---
d0 = int32(FrameStart2 - FrameStart1);  % small, usually within +/-10 after reconcile

if d0 > 0
    Rx2_Synced = [Rx2_Synced(1+d0:end); zeros(d0,1,'single')];
elseif d0 < 0
    dd = -d0;
    Rx1_Synced = [Rx1_Synced(1+dd:end); zeros(dd,1,'single')];
end

fprintf('[SYNC] PeakLag1=%d, PeakLag2=%d, FrameStart1=%d, FrameStart2=%d, d0=%d\n', ...
        PeakLag1, PeakLag2, FrameStart1, FrameStart2, d0);



%% OFDM demod (Idx_Mat 已經 persistent)
Rx1_Time = Rx1_Synced(Idx_Mat);
Rx2_Time = Rx2_Synced(Idx_Mat);

Rx1_Freq = fftshift(fft(Rx1_Time) / sqrt(Nfft), 1);
Rx2_Freq = fftshift(fft(Rx2_Time) / sqrt(Nfft), 1);

%% 定義 Active subcarriers (1644x560)
Rx1_Pure = [Rx1_Freq(203:1024, :); Rx1_Freq(1026:1847, :)];
Rx2_Pure = [Rx2_Freq(203:1024, :); Rx2_Freq(1026:1847, :)];

%% 取出 DMRS
Y1_DMRS = Rx1_Pure(SCIdx_DMRS_All, SymIdx_DMRS);
Y2_DMRS = Rx2_Pure(SCIdx_DMRS_All, SymIdx_DMRS);

Y1_tilde = Y1_DMRS ./ Seq_DMRS;
Y2_tilde = Y2_DMRS ./ Seq_DMRS;

%% Channel estimation
H11_dmrs = W1 * Y1_tilde; % LMMSE 頻域內插
H12_dmrs = -(W2 * Y1_tilde);
H21_dmrs = W1 * Y2_tilde;
H22_dmrs = -(W2 * Y2_tilde);

H11 = H11_dmrs * Mat_Interp_Time.'; % 時域線性內插
H12 = H12_dmrs * Mat_Interp_Time.';
H21 = H21_dmrs * Mat_Interp_Time.';
H22 = H22_dmrs * Mat_Interp_Time.';

Y1 = Rx1_Pure;
Y2 = Rx2_Pure;

%% SNR estimation
H11p = H11_dmrs(SCIdx_DMRS_All, :);
H12p = H12_dmrs(SCIdx_DMRS_All, :);
H21p = H21_dmrs(SCIdx_DMRS_All, :);
H22p = H22_dmrs(SCIdx_DMRS_All, :);

Y1_hat_tilde = H11p + alt822 .* H12p;
Y2_hat_tilde = H21p + alt822 .* H22p;

N1 = Y1_tilde - Y1_hat_tilde;
N2 = Y2_tilde - Y2_hat_tilde;

SignalPower_11 = mean(abs(H11p).^2, 'all');   % Tx1 -> Rx1 path power
NoisePower1  = mean(abs(N1).^2,  'all');    % Rx1 noise power
SNR1_Frame = 10*log10( SignalPower_11 / NoisePower1 );

SignalPower_22 = mean(abs(H22p).^2, 'all');   % Tx2 -> Rx2 path power
NoisePower2  = mean(abs(N2).^2,  'all');    % Rx2 noise power
SNR2_Frame = 10*log10( SignalPower_22 / NoisePower2 );

s1 = single(1/sqrt(NoisePower1));
s2 = single(1/sqrt(NoisePower2));


%% LMMSE detector (lmmse_prac，應該可以改寫成更好的形式，分成四個變數讀取會比一個快)
H_mex = [complex(single(H11(idx_lin))).'; ...
         complex(single(H21(idx_lin))).'; ...
         complex(single(H12(idx_lin))).'; ...
         complex(single(H22(idx_lin))).'];

H_mex(1,:) = H_mex(1,:) * s1;  % H11 -> Rx1
H_mex(3,:) = H_mex(3,:) * s1;  % H12 -> Rx1
H_mex(2,:) = H_mex(2,:) * s2;  % H21 -> Rx2
H_mex(4,:) = H_mex(4,:) * s2;  % H22 -> Rx2

Y_mex = [complex(single(Y1(idx_lin))).'; ...
         complex(single(Y2(idx_lin))).'];

Y_mex(1,:) = Y_mex(1,:) * s1;
Y_mex(2,:) = Y_mex(2,:) * s2;

X_mex = lmmse_single(H_mex, Y_mex, single(1.0));

Sym1 = X_mex(1,:).';
Sym2 = X_mex(2,:).';

%% LLR Start
LLR1 = qamdemod(Sym1, M, 'UnitAveragePower', true, 'OutputType', 'approxllr', 'NoiseVariance', NoisePower1);
LLR2 = qamdemod(Sym2, M, 'UnitAveragePower', true, 'OutputType', 'approxllr', 'NoiseVariance', NoisePower2);

LLR1 = LLR1(invPerm);
LLR2 = LLR2(invPerm);

num_blocks = ceil(double(data_length) / double(NUM_message));
coded_len  = num_blocks * NUM_VN;
bits_per_layer = coded_len / 2;

LLR1v = LLR1(1:bits_per_layer);
LLR2v = LLR2(1:bits_per_layer);

LLR_coded = zeros(coded_len, 1);
LLR_coded(1:2:end) = LLR1v;
LLR_coded(2:2:end) = LLR2v;

LLR_Mat = reshape(LLR_coded, NUM_VN, num_blocks);

Bits_QAM = (LLR_Mat(1:NUM_message, :) < 0);
Bits_QAM = double(Bits_QAM(:));
Bits_QAM = Bits_QAM(1:data_length);

Time_QAM = toc(t0);

%% LLR -> Soft decision (using "ldpc_decoder_prac.c")
Bits_LDPC_all = zeros(NUM_message * num_blocks, 1);
for b = 1:num_blocks
    llr_in = double(LLR_Mat(:, b));
    llr_out = ldpc_decoder_prac(llr_in, h_LDPC_rows, h_LDPC_cols, double(NUM_CN), double(NUM_VN), double(maxNumIter)); %% 6 inputs
    Bits_LDPC_all((b-1)*NUM_message + (1:NUM_message)) = double(llr_out(1:NUM_message) < 0);
end
Bits_LDPC = Bits_LDPC_all(1:data_length);

Time_Frame = toc(t0);

%% BER
Img_Ref = imresize(imread("view2.jpg"), img_size(1:2));
Img_Col_u8 = uint8(Img_Ref(:));
TxBitsMat = zeros(numel(Img_Col_u8), 8);

%% Image Reconstruction
for bb = 1:8
    TxBitsMat(:, bb) = bitget(Img_Col_u8, 9-bb);
end
TxBits = double(TxBitsMat(:));
TxBits = TxBits(1:data_length);

BER_QAM  = sum(Bits_QAM  ~= TxBits) / double(data_length);
BER_LDPC = sum(Bits_LDPC ~= TxBits) / double(data_length);
Success_Flag =  BER_LDPC < 0.05;

%% Throughput Calculation
Throughput_QAM  = (double(data_length) * (1 - BER_QAM)  / max(Time_QAM, 1e-12)) / 1e6;
Throughput_LDPC = (double(data_length) * (1 - BER_LDPC) / max(Time_Frame, 1e-12)) / 1e6;

SNR_avg = (SNR1_Frame+SNR2_Frame)/2;
Capacity_Val = log2(1 + 10^(SNR_avg/10));
Capacity_QAM  = Capacity_Val;
Capacity_LDPC = Capacity_Val;

w = 2.^(7:-1:0).';

numBytes = prod(img_size);
needBits = numBytes * 8;

Bq = reshape(Bits_QAM(1:needBits),  numBytes, 8);
Bl = reshape(Bits_LDPC(1:needBits), numBytes, 8);

u8_qam  = uint8(double(Bq) * w);
u8_ldpc = uint8(double(Bl) * w);

Img_QAM  = reshape(u8_qam,  img_size);
Img_LDPC = reshape(u8_ldpc, img_size);

end
