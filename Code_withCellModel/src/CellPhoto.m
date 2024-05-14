function [F1opt,fphoto_star] = CellPhoto(Irr)
%CellPhoto computes the optimal (maximum) value of f_photo (alpha_photo in
%CellCNP) by balancing the trade off between photosystems I & II.
%   NOTE: because this is a function of only light, the output field will be
%   constant over all iterations of the parameter optimization.
%   So only need to calculate this once before running the CellCNP model.
%
%   INPUTS:
%       Irr is light level in microeinsteins per m^2 per second

% This function calculates the carbon gathering efficiency of 1 unit of
% the photosynthetic apparatus. It is based on a paper by Geider and
% Talmy. This model involves optimizing photosynthetic allocations
% between carbon fixation and photo capture/electron transport.

% The model assumes the photosynthetic apparatus is divided into two
% compartments: F1 stands for the carbon fixing compartment and
% F2 stands for the light harvesting compartment. Associated with each of
% these compartments are rate constants k1 and k2
% F1 and F2 have as units grams of carbon.

% The maximum photosynthesis rate is: Pmax = min([k1*F1,k2*F2]).

% The actual photosynthesis rate depends on the irradiance and the
% investments F1 and F2 according to f_Photo = Pmax*(1-exp(-alpha_ph *
% PhiM * F2 * Irr / Pmax))

F1opt = NaN(length(Irr),1);
neg_fphoto_star = NaN(length(Irr),1);

for ii= 1:length(Irr)
    I = Irr(ii);
    options=optimset('Display','none','TolX',1e-6);

    % finding the maximum fphoto = solving for min of negative fphoto
    [F1opt(ii),neg_fphoto_star(ii)] = fminbnd(@(F1) -1*fphoto(F1,I),0,1,options);
end

fphoto_star = -1*neg_fphoto_star;

% fphoto is the gross photosynthesis rate [gC/gN/s] in Talmy 2013
% converted to [gC/gC/hr] here
% F1 and F2 are the fraction of the light harvesting apparatus (L)
% allocated to carbon fixation (related to Rubisco) (F1) and electron transport (related more to chlorophyll) (F2)
    function fphoto = fphoto(F1,I) %fphoto(F1,Irr)
        % -------- define constants ---------------
        % k1 is the carbon fixation rate constant and k2 is the maximum rate of
        % electron transfer. Values from Talmy et al 2013 : k1 = 6.1e-4 [gC/gN/sec];
        % k2 = 1.4e-3 [gC/gN/sec]

        % convert k1 and k2 to units of [gC/gC/hour] using a constant N:C ratio
        QF=0.17 ;             % Nitrogen/Carbon ratio of functional components (Geider and La Roche 2002)
                              % note: redfield N:C = 16:106 = 0.151
        k1 = 6.1e-4*QF*3600 ; % carbon fixation rate constant [gC/gC/hour]
        k2 = 1.4e-3*QF*3600 ; % maximum rate of electron transfer [gC/gC/hour]

        alpha_ph = 11.6*QF ;  % Light absorbtion [m^2 /gC] (tuned parameter value from Talmy et al 2013)
        phiM = 1.0e-6 ;  %  maximum quantum efficiency [gC/(mu mol photons)] (Falkowski and Raven 1997)

        % ------ Compute photosynthesis rate -------
        Pmax = min([k1*F1,k2*(1-F1)]);
        fphoto   = Pmax*(1-exp(-alpha_ph*phiM*(1-F1)*(I*3600)/Pmax));
		% parametrization of photosynthesis similar Bouman 2018 - chlorophyll normalized photo-physiological parameters
    end

end %end function CellPhoto(Irr)
