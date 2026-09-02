% This function returns the percentage of slow twitch fibers in the muscles
% The data come from Uchida et al. (2016).
% We used 0.5 when no data were available.
%
% Author: Antoine Falisse
% Date: 12/19/2018
%     
function pctst = getSlowTwitchRatios(muscleNames)

    pctst_data.glmed1 = 0.55;
    pctst_data.glmed2 = 0.55;
    pctst_data.glmed3 = 0.55;
    pctst_data.glmin1 = 0.55;
    pctst_data.glmin2 = 0.55;
    pctst_data.glmin3 = 0.55;
    pctst_data.semimem = 0.4925;
    pctst_data.semiten = 0.425;    
    pctst_data.bflh = 0.5425;
    pctst_data.bfsh = 0.529;
    pctst_data.sart = 0.50;
    pctst_data.addmagProx = 0.552;
    pctst_data.addmagMid = 0.552;
    pctst_data.addmagDist = 0.552;
    pctst_data.addmagIsch = 0.552;% we just took the same of the 
        % other adductor magnus bundles
    pctst_data.tfl = 0.50;
    pctst_data.pect = 0.50;
    pctst_data.grac = 0.50;
    pctst_data.glmax1 = 0.55;
    pctst_data.glmax2 = 0.55;
    pctst_data.glmax3 = 0.55;
    pctst_data.iliacus = 0.50;
    pctst_data.psoas = 0.50;
    pctst_data.quadfem = 0.50;
    pctst_data.gem = 0.50;
    pctst_data.piri = 0.50;
    pctst_data.recfem = 0.3865;
    pctst_data.vasmed = 0.503;
    pctst_data.vasint = 0.543;
    pctst_data.vaslat = 0.455;
    pctst_data.gasmed = 0.566;
    pctst_data.gaslat = 0.507;
    pctst_data.soleus = 0.803;
    pctst_data.tibpost = 0.60;
    pctst_data.fdl = 0.60;
    pctst_data.fhl = 0.60;
    pctst_data.tibant = 0.70;    
    pctst_data.perbrev = 0.60;
    pctst_data.perlong = 0.60;
    pctst_data.pertert = 0.75;
    pctst_data.edl = 0.75;
    pctst_data.ehl = 0.75;    
    pctst_data.ercspn_r = 0.60;
    pctst_data.ercspn_l = 0.60;
    pctst_data.intobl_r = 0.56;
    pctst_data.intobl_l = 0.56;
    pctst_data.extobl_r = 0.58; 
    pctst_data.extobl_l = 0.58; 
    pctst_data.addlong = 0.50;
    pctst_data.addbrev = 0.50;
    
    pctst = zeros(length(muscleNames),1);
    for i = 1:length(muscleNames)
        pctst(i,1) = pctst_data.(muscleNames{i});
    end
    
end   
