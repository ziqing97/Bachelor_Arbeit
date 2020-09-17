function Institude = cal_dSdt(Institude)
    % dS based on Recovered anomalie
    S = Institude.Recovered(:,2);
    dS = zeros(length(S)-2,1);
    for i = 1 : length(S)-2
        dS(i) = (S(i+2) - S(i)) / 2;
    end
    dS_Uncertainty = zeros(length(dS),1);
    for i=1:length(dS)-2
        dS_Uncertainty(i) = sqrt((Institude.Uncertainty(i,2)/2)^2 + (Institude.Recovered(i+2,2)/2)^2);
    end
    dS = [Institude.Recovered(2:end-1,1),dS];
    dS_Uncertainty = [Institude.Recovered(2:end-1,1),dS_Uncertainty];
    Institude.dS_t = Institude.Recovered(2:end-1,1);
    Institude.dS = dS;
    Institude.dS_unc = dS_Uncertainty;
end