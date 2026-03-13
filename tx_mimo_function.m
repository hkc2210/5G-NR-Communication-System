function [TxStream1, TxStream2, data_length, Img_Resized, Ref_PSS_Time] = tx_mimo_function(CodeRate)
%% 流程
% 定義 mask
% 決定圖片尺寸
% 輸入圖片轉成 8-bit 然後拉成一長條

%% 參數規格
% 16-QAM
% Nfft = 2048
% Total_Samples = 1228800
%
% PSS amp = 1

%% 參數與 mask
M = 16;
NID_cell = 86;
NID2 = mod(NID_cell, 3);
val_DMRS = single(1);
amp_PSS = 1;

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

Mask_PSS  = false(2048, 560);
Mask_DMRS = false(2048, 560);
Mask_Act  = false(2048, 560);

Mask_PSS(idx_PSS_SC, idx_PSS_Sym) = true;
Mask_DMRS(idx_DMRS_SC, idx_DMRS_Sym) = true;
Mask_Act(idx_Data_SC, :) = true;

Mask_Data = Mask_Act & ~Mask_DMRS & ~Mask_PSS;
Mask_Data(:, idx_SSB_Sym) = false;
Mask_Data(idx_PBCH_SC, idx_PBCH_Sym) = true;
Mask_Data(idx_SSS_SC, idx_SSS_Sym) = true;

TxGrid1 = zeros(2048, 560);
TxGrid2 = zeros(2048, 560);

%% Image Data
Img_Original = imread("view2.jpg");
switch CodeRate
    case 0.50
        Struct_LDPC = load('LDPC_11nD2_1296b_R12.mat');
        ImgSize = 383;
    case 0.67
        Struct_LDPC = load('LDPC_11nD2_1296b_R23.mat');
        ImgSize = 442;
    case 0.75
        Struct_LDPC = load('LDPC_11nD2_1296b_R34.mat');
        ImgSize = 469;
    case 0.83
        Struct_LDPC = load('LDPC_11nD2_1296b_R56.mat');
        ImgSize = 494;
end

Img_Resized  = imresize(Img_Original, [ImgSize, ImgSize]); % 494, 469, 442, 383

Img_Col_u8 = uint8(Img_Resized(:));
Img_Bit_Matrix = zeros(numel(Img_Col_u8), 8);
for b = 1:8
    Img_Bit_Matrix(:, b) = bitget(Img_Col_u8, 9-b);
end
raw_bits = Img_Bit_Matrix(:);
data_length = length(raw_bits);

%% LDPC (11n D2 1296, R=1/2, 2/3, 3/4, 5/6)
LDPC = Struct_LDPC.LDPC;

H = sparse(LDPC.H.x ~= 0); % 矩陣本來就是稀疏的，這行是讓資料型態也能稀疏儲存
encCfg = ldpcEncoderConfig(H);
K = encCfg.NumInformationBits; % 區塊 codelength
num_blocks  = ceil(double(data_length) / double(K)); % 區塊數
pad_len     = num_blocks * K - double(data_length); % 要補上的 0 數

bits_padded = [raw_bits; zeros(pad_len, 1)];
bits_matrix = reshape(bits_padded, K, num_blocks);

coded_matrix = ldpcEncode(logical(bits_matrix), encCfg); % 使用內建 function 完成 encode
coded_bits   = coded_matrix(:);

%% Interleaver
bits_layer1 = coded_bits(1:2:end);
bits_layer2 = coded_bits(2:2:end);

%% 依 Mask_Data 決定當前 frame 需要的 bits 並做 16QAM
Ndata = nnz(Mask_Data);
k = log2(M);
Nb_layer = Ndata * k;

b1 = zeros(Nb_layer, 1);
b2 = zeros(Nb_layer, 1);

L1 = min(length(bits_layer1), Nb_layer);
L2 = min(length(bits_layer2), Nb_layer);

b1(1:L1) = bits_layer1(1:L1);
b2(1:L2) = bits_layer2(1:L2);

%% Data permutation
load('Permutation.mat');
b1 = b1(Permutation);
b2 = b2(Permutation);

sym1 = qammod(b1, M, 'InputType', 'bit', 'UnitAveragePower', true);
sym2 = qammod(b2, M, 'InputType', 'bit', 'UnitAveragePower', true);

TxGrid1(Mask_Data) = sym1;
TxGrid2(Mask_Data) = sym2;

%% DMRS (CDM + Scramble)
load('Scrambled.mat');
Phase_DMRS = reshape(Phase_DMRS, length(idx_DMRS_SC), length(idx_DMRS_Sym));
Seq_DMRSScrambled = val_DMRS .* exp(1j * Phase_DMRS);

alt = ones(length(idx_DMRS_SC), 1);
alt(2:2:end) = -1;

TxGrid1(idx_DMRS_SC, idx_DMRS_Sym) = Seq_DMRSScrambled;
TxGrid2(idx_DMRS_SC, idx_DMRS_Sym) = repmat(alt, 1, length(idx_DMRS_Sym)) .* Seq_DMRSScrambled;

%% PSS
Seq_PSS = PSS_Function(NID2);
TxGrid1(idx_PSS_SC, idx_PSS_Sym) = amp_PSS * Seq_PSS;
TxGrid2(idx_PSS_SC, idx_PSS_Sym) = amp_PSS * Seq_PSS;

tx_grid = cat(3, TxGrid1, TxGrid2);

%% IFFT + CP
x1 = ifft(ifftshift(TxGrid1, 1)) * sqrt(2048);
x2 = ifft(ifftshift(TxGrid2, 1)) * sqrt(2048);

CP = 144 * ones(1, 560);
CP(1:28:end) = 208;

symLen = 2048 + CP;

TxStream1 = zeros(1, 1228800);
TxStream2 = zeros(1, 1228800);

start = [0, cumsum(symLen(1:end-1))];

for i = 1:560
    cp = CP(i);
    L  = symLen(i);
    s  = start(i);

    idx_out = (s+1):(s+L);

    TxStream1(idx_out) = [x1(2048-cp+1:2048, i); x1(:, i)].';
    TxStream2(idx_out) = [x2(2048-cp+1:2048, i); x2(:, i)].';
end

scale = 1 / sqrt(2);
TxStream1 = TxStream1 * scale;
TxStream2 = TxStream2 * scale;

%% Ref PSS time
Grid_Ref = zeros(2048, 1);
Grid_Ref(idx_PSS_SC) = amp_PSS * Seq_PSS;
Ref_PSS_Time = ifft(ifftshift(Grid_Ref)) * sqrt(2048);

end