%
%  Comparing Msk between 48layer model and 24layer model and make the new msk for each basin.
%  (Because horizontal msk is little bit different between 24layer model and 48layer model).

clc; close all; clear all;
addpath('../../DATA/BGC_2023Nature/');
addpath('../../DATA/BGC_48layer/');

load OCIM2_CTL_He.mat
msk_24layer = output.M3d;
msk_section_24layer = output.MSKS;
ATL_24layer_sfc = msk_section_24layer.ATL(:,:,1);
PAC_24layer_sfc = msk_section_24layer.PAC(:,:,1);
IND_24layer_sfc = msk_section_24layer.IND(:,:,1);
ARC_24layer_sfc = msk_section_24layer.ARC(:,:,1);
MED_24layer_sfc = msk_section_24layer.MED(:,:,1);
clear output sol;


load OCIM2_CTL_He_48layer.mat
clear TR
msk_48layer = output.M3d;

msk_24layer_sfc = msk_24layer(:,:,1);
msk_48layer_sfc = msk_48layer(:,:,1);


%48layer msk에서 한번 육지이면 계속 육지인지? dry-wet-dry box일수는 없는지?
dmsk_48 = nan(91,180,47);
for i = 1:size(msk_48layer, 3)-1;
    dmsk_48(:,:,i) = msk_48layer(:,:,i+1) - msk_48layer(:,:,i);
end

% 깊어질때 육지였다가(0) wet box(1)인 지점 여부 확인
[row, col, page] = ind2sub(size(dmsk_48), find(dmsk_48 == 1));
dmskpositive= [row, col, page];
% >>>> 없음. 즉, 한번 육지이면 계속 육지. 24layer surface 사용 가능.


% 24 layer msk와 표층에서 어디가 다른지 확인
chmsk = msk_24layer_sfc - msk_48layer_sfc;
ichmsk = find(chmsk (:));
[lat, lon] = ind2sub(size(chmsk), ichmsk);


for i = 1:length(ichmsk);
    valuechmsk(i) = chmsk(ichmsk(i));
end
valuechmsk = valuechmsk'

msk_difference = [lat, lon, valuechmsk];
%msk_difference: 1열: 위도, 2열: 경도, 3열: 24layer msk - 48layer msk
%  1: 24layer에서는 ocean, 48layer에서는 land
% -1: 24layer에서는  land, 48layer에서는 ocean
%  1인 경우 >> 24 layer의 지역별 msk에서 해당 위경도의 1 값을 0으로 바꿈.
% -1인 경우 >> 24 layer의 지역별 msk에서 해당 위경도의 0 값을 1으로 바꿈.

%iATL = find(ATL_24layer_sfc);
%iPAC = find(PAC_24layer_sfc);
%iIND = find(IND_24layer_sfc);
%iARC = find(ARC_24layer_sfc);
%iMED = find(MED_24layer_sfc);

%[ATLlat, ATLlon] = ind2sub(size(ATL_24layer_sfc), iATL);
%[PAClat, PAClon] = ind2sub(size(PAC_24layer_sfc), iPAC);
%[INDlat, INDlon] = ind2sub(size(IND_24layer_sfc), iIND);
%[ARClat, ARClon] = ind2sub(size(ARC_24layer_sfc), iARC);
%[MEDlat, MEDlon] = ind2sub(size(MED_24layer_sfc), iMED);

%msk_ATL = [ATLlat, ATLlon];
%msk_PAC = [PAClat, PAClon];
%msk_IND = [INDlat, INDlon];
%msk_ARC = [ARClat, ARClon];
%msk_MED = [MEDlat, MEDlon];


%for i = 1:length(ichmsk);
   % if      any( ichmsk(i) == iATL(:) );
                  %region_msk {i,1} = 'ATL';
   % elseif any( ichmsk(i) == iPAC(:) );
                  %region_msk {i,1} = 'PAC';
   % elseif any( ichmsk(i) == iIND(:) );
                  %region_msk {i,1} = 'IND';
   % elseif any( ichmsk(i) == iARC(:) );
                  %region_msk {i,1} = 'ARC';
   % elseif any( ichmsk(i) == iMED(:) );
                  %region_msk {i,1} = 'MED';
   % else          region_msk {i,1} = 'Ground';
   % end 
%end

   
%위의 Ground는 pcolor통해서 일일히 확인. 그 결과, 21개중 4개 Arctic, 1개ATL 나머지 MED.
%24 layer의 각 basin별 surface msk를 바꾼 후, 48 layer msk에서 표층 위경도에 해당하면 그 지역.

%ATL layer 바꾸기.
ATL_48layer_sfc = ATL_24layer_sfc;
ATL_48layer_sfc(69, 150) = 1; 

%ARC layer 바꾸기.
ARC_48layer_sfc = ARC_24layer_sfc;
ARC_48layer_sfc(83, 124) = 1;
ARC_48layer_sfc(83, 125) = 1;
ARC_48layer_sfc(84, 124) = 1;
ARC_48layer_sfc(84, 125) = 1;

%MED layer 바꾸기.
MED_48layer_sfc = MED_24layer_sfc;
MED_48layer_sfc(65, 1)   = 1;
MED_48layer_sfc(65, 2)   = 1;
MED_48layer_sfc(65, 177) = 1;
MED_48layer_sfc(65, 178) = 1;
MED_48layer_sfc(65, 179) = 1;
MED_48layer_sfc(65, 180) = 1;

MED_48layer_sfc(63, 1)   = 0;
MED_48layer_sfc(63, 2)   = 0;
MED_48layer_sfc(63, 3)   = 0;
MED_48layer_sfc(63, 4)   = 0;
MED_48layer_sfc(63, 5)   = 0;
MED_48layer_sfc(66, 5)   = 0;
MED_48layer_sfc(63, 177)   = 0;
MED_48layer_sfc(63, 178)   = 0;
MED_48layer_sfc(63, 179)   = 0;
MED_48layer_sfc(63, 180)   = 0;

%PAC layer 설정.
PAC_48layer_sfc = PAC_24layer_sfc;

%IND layer 설정.
IND_48layer_sfc = IND_24layer_sfc;

%24layer로부터 바꿔서 새로 얻은 msk 합하여 기존 msk와 비교
New_msk_48layer_sfc = nan(91,180);
New_msk_48layer_sfc = PAC_48layer_sfc + ATL_48layer_sfc + IND_48layer_sfc + ARC_48layer_sfc + MED_48layer_sfc;

%기존 msk와 surface 자료 비교 >> 동일. 따라서, 바뀐 msk의 1에 해당하는 지역들을 그 지역 msk로 설정.
dnewold = New_msk_48layer_sfc - msk_48layer_sfc;
[drow, dcol] = ind2sub(size(dnewold), find(dnewold ~= 0));



% ATL 지역 설정
ATL_48layer = zeros(size(msk_48layer)); %91x180x48
msk_ATL_sfc = zeros(size(msk_48layer_sfc)); %91x180
for i = 1:size(msk_48layer, 3);
    msk_ATL_sfc = (msk_48layer_sfc == ATL_48layer_sfc);
    ATL_48layer(:,:,i) = msk_48layer(:,:,i) .* msk_ATL_sfc;
end

% PAC 지역 설정
PAC_48layer = zeros(size(msk_48layer)); %91x180x48
msk_PAC_sfc = zeros(size(msk_48layer_sfc)); %91x180
for i = 1:size(msk_48layer, 3);
    msk_PAC_sfc = (msk_48layer_sfc == PAC_48layer_sfc);
    PAC_48layer(:,:,i) = msk_48layer(:,:,i) .* msk_PAC_sfc;
end

% IND 지역 설정
IND_48layer = zeros(size(msk_48layer)); %91x180x48
msk_IND_sfc = zeros(size(msk_48layer_sfc)); %91x180
for i = 1:size(msk_48layer, 3);
    msk_IND_sfc = (msk_48layer_sfc == IND_48layer_sfc);
    IND_48layer(:,:,i) = msk_48layer(:,:,i) .* msk_IND_sfc;
end

% ARC 지역 설정
ARC_48layer = zeros(size(msk_48layer)); %91x180x48
msk_ARC_sfc = zeros(size(msk_48layer_sfc)); %91x180
for i = 1:size(msk_48layer, 3);
    msk_ARC_sfc = (msk_48layer_sfc == ARC_48layer_sfc);
    ARC_48layer(:,:,i) = msk_48layer(:,:,i) .* msk_ARC_sfc;
end

% MED 지역 설정
MED_48layer = zeros(size(msk_48layer)); %91x180x48
msk_MED_sfc = zeros(size(msk_48layer_sfc)); %91x180
for i = 1:size(msk_48layer, 3);
    msk_MED_sfc = (msk_48layer_sfc == MED_48layer_sfc);
    MED_48layer(:,:,i) = msk_48layer(:,:,i) .* msk_MED_sfc;
end


%지역별로 얻은 48layer msk와 기존 msk와 비교
New_msk_48layer = nan(91,180,48);
New_msk_48layer = PAC_48layer + ATL_48layer + IND_48layer + ARC_48layer + MED_48layer;

%기존 msk와 surface 자료 비교
dnewold_full = New_msk_48layer - msk_48layer;
[drowfull, dcolfull] = ind2sub(size(dnewold_full), find(dnewold_full ~= 0));

%msk자료 저장
MSKS.ATL_48layer = ATL_48layer;
MSKS.PAC_48layer = PAC_48layer;
MSKS.IND_48layer = IND_48layer;
MSKS.ARC_48layer = ARC_48layer;
MSKS.MED_48layer = MED_48layer;

%DATA file deirectory에 저장하기
fileName = 'MSKS_48layer.mat'
directory = '../../DATA/BGC_48layer'
filePath = fullfile(directory, fileName);
save(filePath, 'MSKS');






