function Institude = cal_dSdt(Institude)
    % dS based on Recovered anomalie
    dS = diff(Institude.Recovered(:,2));
    dS_Uncertainty = zeros(length(dS),1);
    for i=1:length(dS)
        dS_Uncertainty(i) = sqrt(Institude.Recovered(i,2)^2 + Institude.Recovered(i+1,2)^2);
    end
    dS = [Institude.Recovered(2:end,1),dS];
    dS_Uncertainty = [Institude.Recovered(2:end,1),dS_Uncertainty];
    Institude.dS = dS;
    Institude.dS_unc = dS_Uncertainty;
end