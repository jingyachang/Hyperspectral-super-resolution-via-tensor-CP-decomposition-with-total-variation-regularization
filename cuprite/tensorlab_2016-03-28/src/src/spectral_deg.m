function Pm = spectral_deg(DATA,sensor)

% SPECTRAL_DEG makes spectral degradation matrix from SRI following LANDSAT's spectral degradation
% Pm = SPECTRAL_DEG(DATA) returns Pm from SRI DATA
% 
% INPUT ARGUMENTS:
%     DATA: SRI of size ImxJmxKh
%     sensor: LANDSAT or Quickbird
% OUTPUT ARGUMENTS:
%     Pm: spectral degradation matrix of size KmxKh
% 
% SEE ALSO: SPATIAL_DEG
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr


if sensor=="LANDSAT"
    Km = 6;
    spectral_distrib = linspace(400,2500,size(DATA,3)); %Hyp: uniform distribution

    bands = [450 520; 520 600; 630 690; 760 900; 1550 1770; 2080 2350];
elseif sensor=="Quickbird"
    Km = 4;
    spectral_distrib = linspace(430,860,size(DATA,3)); %Hyp: uniform distribution
    bands = [430 545; 466 620; 590 710; 715 918];
end

Pm = zeros(Km,size(DATA,3));

for k=1:Km
    ind = find(spectral_distrib >= bands(k,1) & spectral_distrib <= bands(k,2));
    Pm(k,ind) = 1/length(ind);
end

end

