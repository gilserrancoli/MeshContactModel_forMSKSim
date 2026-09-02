% This function returns the specific tensions in the muscles
% The data come from Uchida et al. (2016).
%
% Author: Antoine Falisse
% Date: 12/19/2018
%         
function sigma = getSpecificTensions(muscleNames)

    sigma_data.glmed1 = 0.74455;
    sigma_data.glmed2 = 0.75395;
    sigma_data.glmed3 = 0.75057;
    sigma_data.glmin1 = 0.75;
    sigma_data.glmin2 = 0.75;
    sigma_data.glmin3 = 0.75116;
    sigma_data.semimem = 0.62524;
    sigma_data.semiten = 0.62121;    
    sigma_data.bflh = 0.62222;
    sigma_data.bfsh = 1.00500;
    sigma_data.sart = 0.74286;
    sigma_data.addmagProx = 0.55217;
    sigma_data.addmagMid = 0.55323;
    sigma_data.addmagDist = 0.54831;
    sigma_data.addmagIsch = 0.55323; % we just took the maximum of the 
        % other adductor magnus bundles
    sigma_data.tfl = 0.75161;
    sigma_data.pect = 0.76000;
    sigma_data.grac = 0.73636;
    sigma_data.glmax1 = 0.75395;
    sigma_data.glmax2 = 0.74455;
    sigma_data.glmax3 = 0.74595;
    sigma_data.iliacus = 1.2477;
    sigma_data.psoas = 1.5041;
    sigma_data.quadfem = 0.74706;
    sigma_data.gem = 0.74545;
    sigma_data.piri = 0.75254;
    sigma_data.recfem = 0.74936;
    sigma_data.vasmed = 0.49961;
    sigma_data.vasint = 0.55263;
    sigma_data.vaslat = 0.50027;
    sigma_data.gasmed = 0.69865;
    sigma_data.gaslat = 0.69694;
    sigma_data.soleus = 0.62703;
    sigma_data.tibpost = 0.62520;
    sigma_data.fdl = 0.5;
    sigma_data.fhl = 0.50313;
    sigma_data.tibant = 0.75417;    
    sigma_data.perbrev = 0.62143;
    sigma_data.perlong = 0.62450;
    sigma_data.pertert = 1;
    sigma_data.edl = 0.75294;
    sigma_data.ehl = 0.73636;    
    sigma_data.ercspn_r = 0.25;
    sigma_data.ercspn_l = 0.25;
    sigma_data.intobl_r = 0.25;
    sigma_data.intobl_l = 0.25;
    sigma_data.extobl_r = 0.25;  
    sigma_data.extobl_l = 0.25;  
    sigma_data.addlong = 0.74643;
    sigma_data.addbrev = 0.75263;
    
    sigma = zeros(length(muscleNames),1);
    for i = 1:length(muscleNames)
        sigma(i,1) = sigma_data.(muscleNames{i});
    end
    
end    

