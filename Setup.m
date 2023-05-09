Ti=0;Tf=1800;%in seconds. 1800s=30min.
Ntme_discr=100;Ntme=(Tf*Ntme_discr)+1;%is 0.01s i.e. 10ms
t_data=linspace(Ti,Tf,Ntme);

%% Spine neck morphological plasticity
if NeckMorphPlast==0
   r_neck_time_array=0.0564*ones(1,length(t_data));
   l_neck_time_array=0.3*ones(1,length(t_data));
elseif NeckMorphPlast==1
   %r increases but l decreases with time constant tauNM 
   r_neck_time_array=0.0564+((0.1262-0.0564)*(1-exp(-1*t_data./tauNM)));
   l_neck_time_array=0.1250+((0.3-0.1250)*exp(-1*t_data./tauNM));
elseif NeckMorphPlast==2
   %r decreases but l increases with time constant tauNM
   r_neck_time_array=0.03+((0.0564-0.03)*exp(-1*t_data./tauNM));
   l_neck_time_array=0.3+((1-0.3)*(1-exp(-1*t_data./tauNM)));
end
%% Spine neck septin plasticity
if NeckSepPlast==0
    D_neck_time_array=6.7e-3*ones(1,length(t_data));
elseif NeckSepPlast==1
    %D only increases to D_esm=0.1um2.s-1, with time constant tauNS
    D_neck_time_array=(6.7e-3)+((0.1-(6.7e-3))*(1-exp(-1*t_data./tauNS)));
end

%% Phos rates in stimulated spine head and dendritic shaft
N_spine=100;N_dend=1000;
Twin=60;%Stimulation window in s.
tauCaMKII=20;%Self-decay time constant of CaMKII in s.

Ephos_spine=(1.111e-5)*ones(N_spine,length(t_data));
Ephos_PSDB=(9e-4)*ones(N_spine,length(t_data));
Ephos_neck=(1.111e-5)*ones(N_spine,length(t_data));
Ephos_dend=(1.111e-5)*ones(N_dend,length(t_data));

Ephos_spine(StimSpineLoc,1:((Twin*100)+1))=0.111;
Ephos_PSDB(StimSpineLoc,1:((Twin*100)+1))=9;
Ephos_spine(StimSpineLoc,((Twin*100)+1):end)=(1.111e-5)+(0.111-(1.111e-5))*exp(-1*(t_data(((Twin*100)+1):end)-Twin)./tauCaMKII);
Ephos_PSDB(StimSpineLoc,((Twin*100)+1):end)=(9e-4)+(9-(9e-4))*exp(-1*(t_data(((Twin*100)+1):end)-Twin)./tauCaMKII);

%% Spatiotemporal profile of CaN in dendritic branch
tauCaN=100;%Self-decay time constant of Calcineurin in s.
CaNSignalling;
