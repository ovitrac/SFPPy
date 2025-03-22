% This script shows the prediction of difusivities of toluene using hole Free-Volume theory
% at different temperatures in different polymers (LDPE, PMMA, PS, PVAc, PET).
% gPET corresponds to dry PET (unplasticized) and wPET coresponds to plasticized PET


% Data have been extracted from the script figure8_revSep2019_2025.m used to generate
% Figure 9 for the reference mentioned below.

% REFERENCE
% Zhu Y., Welle, F. and Vitrac O. A blob model to parameterize polymer hole free volumes and solute diffusion",
% *Soft Matter* **2019**, 15(42), 8912-8932. DOI: https://doi.org/10.1039/C9SM01556F
% 
% ABSTRACT
% Solute diffusion in solid polymers has tremendous applications in packaging,
% reservoir, and biomedical technologies but remains poorly understood. Diffusion
% of non-entangled linear solutes with chemically identical patterns (blobs) deviates
% dramatically in polymers in the solid-state (αlin > 1, Macromolecules 2013, 46, 874)
% from their behaviors in the molten state (αlin = 1, Macromolecules, 2007, 40, 3970).
% This work uses the scale invariance of the diffusivities, D, of linear probes
% D(N·Mblob + Manchor,T,Tg) = N−αlin(T,Tg)D(Mblob + Manchor,T,Tg) comprising
% N identical blobs of mass Mblob and possibly one different terminal pattern (anchor of
% mass Manchor) to evaluate the amounts of hole-free volume in seven polymers (aliphatic,
% semi-aromatic and aromatic) over a broad range of temperatures (−70 K ≤ T − Tg ≤ 160 K).
% The new parameterization of the concept of hole-free volumes opens the application of
% the free-volume theory (FVT) developed by Vrentas and Duda to practically any polymer,
% regardless of the availability of free-volume parameters. The quality of the estimations
% was tested with various probes including n-alkanes, 1-alcohols, n-alkyl acetates, and
% n-alkylbenzene. The effects of enthalpic and entropic effects of the blobs and the anchor
% were analyzed and quantified. Blind validation of the reformulated FVT was tested 
% successfully by predicting from first principles the diffusivities of water and toluene
% in amorphous polyethylene terephthalate from 4 °C to 180 °C and in various other polymers.
% The new blob model would open the rational design of additives with controlled diffusivities
% in thermoplastics.


% Constants
R = 8.31;
T0K = 273.15; % K

data = cell2table(...
   { % Tg are in K
    'LDPE'     148.15    1.87e-08     0.615      3     144    40    0    0.5
    'PMMA'     381.15    1.87e-08      0.56      2     252    65    0    0.5 % original better
%    'PMMA'    381.15     4.8e-08      0.56      2     144    40    0    0.5
%    'PS'       373.15    1.87e-08     0.584        1     112    90    0    0.5
    'PS'       373.15     4.8e-08    0.584      2     144    40    0    0.5  % variant better
    'PVAc'     305.15    1.87e-08      0.86      4     142    40    0    0.5
    'gPET'      349.15    1.0205e-08    0.6761    5     252    65   0    0.6153 % dry
    'wPET'      316.15    1.02046e-08   0.6761    5     252    65   0    0.277734 % wet
    }, ...
    'VariableNames',{'polymer','Tg','D0','xi','ref','Ka','Kb','E','r'})

% references
references = {
'Vrentas and Vrentas, 1994'
'Zielinski and Duda, 1992'
'Lutzow \itet al.\rm, 1999'
'Hong, 1995'
% toluene in PET
'Welle,2008'
'Pennarun \itet al.\rm, 2004'
'Welle,2013'
'our work (permeation)'
'our work (sorption)'
};


lookup = @(P,p) data.(p)(find(ismember(data.polymer,P),1,'first')); 
alpha = @(T,Tg,Ka,Kb) 1+ Ka/(T-Tg+Kb);% for T>=Tg 
alphag = @(T,Tg,Ka,Kb,r) 1+Ka/(r*(T-Tg)+Kb); % for T<Tg
deltaT = 2; % (K) sharpness of the transition at Tg
H = @(T,Tg) 1/2*(1+tanh(4/deltaT*(T-Tg))); % Heaviside function
alphaT = @(T,Tg,Ka,Kb,r) H(Tg,T)*alphag(T,Tg,Ka,Kb,r) + H(T,Tg)*alpha(T,Tg,Ka,Kb); % note H(Tg,T) = - H(T,Tg)
betalin = 1;
Plike = @(T,P) (alphaT(T,lookup(P,'Tg'),lookup(P,'Ka'),lookup(P,'Kb'),lookup(P,'r'))+betalin)/0.24;
DscalingPlike = @(T,P) lookup(P,'D0') * ...
    exp(-lookup(P,'E')/(R*T)) * exp(-lookup(P,'xi')*Plike(T,P));

%% example
test = cell2table( {
  %  polymer       T          D        Dpred
  % _________    ______    ________    ____________
     'PS'        451.15       4e-11     3.1418e-11 
     'PS'        433.15     1.5e-11     2.5938e-11 
     'PS'        413.15       2e-12     1.8558e-11 
     'PS'        383.15     1.5e-13     9.9371e-12 
     'PMMA'      433.15       9e-13     1.1549e-12 
     'PMMA'      413.15       2e-13     4.0974e-13 
     'PMMA'      403.15       8e-14     2.0413e-13 
     'LDPE'      343.15       7e-11     2.3129e-11 
     'PVAc'      383.15       2e-12     1.9349e-13 
     'PVAc'      353.15       1e-13     4.4486e-14 
     'PVAc'      313.15       8e-16      3.594e-16 
     'PS'        388.15     7.2e-13     7.7047e-12 
     'PS'        451.15       4e-11     8.3385e-12 
     'PS'        433.15     1.5e-11     4.9264e-12 
     'PS'        413.15       2e-12     1.9212e-12 
     'PS'        383.15     1.5e-13     1.4127e-13 
     'PS'        388.15     7.2e-13     1.4933e-13 
     'PMMA'      433.15       9e-13     4.5605e-12 
     'PMMA'      413.15       2e-13     1.6536e-12 
     'PMMA'      403.15       8e-14       7.79e-13 
     'wPET'      313.15     4.1e-16     2.3712e-16 
     'gPET'      313.15       2e-18     2.3215e-18 
     'gPET'      313.15     3.8e-18     2.3215e-18 
     'gPET'      333.15     2.6e-17     9.3638e-17 
     'wPET'      313.15     6.3e-17     2.3712e-16 
     'wPET'      313.15     2.1e-16     2.3712e-16 
     'wPET'      313.15     4.8e-16     2.3712e-16  % 27
     'wPET'      323.15     4.1e-16     4.0409e-16 
     'wPET'      323.15     6.1e-16     4.0409e-16 
     'gPET'      373.15    3.26e-15     1.2513e-14 
     'gPET'      453.15    2.78e-12     5.4614e-13 
     'gPET'      434.15    1.17e-12     3.2078e-13 
     'gPET'      414.15     1.9e-13     1.5487e-13 
     'gPET'      394.15    2.47e-14     5.7378e-14 
     },'VariableNames',{'polymer','T','Dexp','Dpred'});

ntests = height(test);
Dtest = zeros(ntests,1);
for i=1:ntests
    T = test.T(i);    
    P = test.polymer{i};
    Dtest(i) = DscalingPlike(T,P);
end

figure, plot(test.Dexp,Dtest,'o','MarkerFaceColor','r')
set(gca,'xscale','log','yscale','log')
hr=refline(1,0); set(hr,'color','k')
xlabel('experimental D values (m2/s)')
ylabel('predicted D values (m2/s)')
title('DFV vs. experimental D values')