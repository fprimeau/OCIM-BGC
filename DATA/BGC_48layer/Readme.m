% 코드 및 Rawdata는 Weilei Nature 2023과 동일하게. 그러나 자료 update해서.
% DIP, DIC, ALK, O2 were downloaded from the GLODAPv2.2023으로. 이전에는
% 2021이었던듯.

% DOC는 따로 reference 있음.

%Climatological ocean temperature and silicate are from World Ocean Atlas
%2018.

%MSKS_48layer: 24layer와 조금 다름(21개). 새로 만든 MSKS.

%pme: 새로 만듬, 91x180 size. weilei code에서는 24layer의 iwet x 1 (91x180x24 기준 iwet)로 되어 있음.

%dic_initial: Weilei의 C14방법을 이용한 (DIC - DICant) field, 즉,
%preindustrial_DIC를 interpolation하여 얻은 48layer 값. 즉, 모든 wet grid box에 값이
%있다. 추후 이를 이용하여 model를 1번 run하여 parameter를 얻고, 이후 transient run...? 여기부터 잘
%이해 안감.