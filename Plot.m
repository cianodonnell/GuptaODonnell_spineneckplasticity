%% Plot change in stimulated spine head surface area (Exp. data: Patterson et al. 2010)
figure(1);
plot(t_data((100*[0:1800])+1)./60,DeltaPrcnt_SpHdSurfArea_ExpData((100*[0:1800])+1),'k','LineWidth',1);
xlabel('Time (min.)');ylabel('\DeltaSp. Hd. Surf. Area (%)');
title('% change in stimulated spine head surface area');
%% Plot CaMKII profile in stimulated spine head
figure(2);
plot(t_data((100*[0:1800])+1)./60,Ephos_spine(StimSpineLoc,(100*[0:1800])+1),'k','LineWidth',1);
xlabel('Time (min.)');ylabel('k_{phos} (s^{-1})');
title({'CaMKII-dependent AMPAR phos. rate const.', 'in stimulated spine head at 50\mum'});
%% Plot CaN profile in stimulated spine head
figure(3);
subplot(1,2,1);
plot(t_data((100*[0:1800])+1)./60,Edephos_spine(StimSpineLoc,(100*[0:1800])+1),'b','LineWidth',1);
xlim([0,30]);ylim([0,30]);
xlabel('Time (min.)');ylabel('[CaN](\muM)');
title({'Active CaN Conc. in', 'stimulated spine head at 50\mum'});
subplot(1,2,2);
plot(t_data((100*[0:1800])+1)./60,Edephos_dend(StimDendLoc,(100*[0:1800])+1),'k','LineWidth',1);
xlim([0,30]);ylim([0,30]);
xlabel('Time (min.)');ylabel('[CaN](\muM)');
title({'Active CaN Conc. in dendrite','underneath stimulated spine head'});
%% Plot change in total stimulated spine head surface AMPAR fluroresence
load('DeltaPrcnt_TotSpHdSurfFluor_Expdata.mat');
figure(4);
plot(t_rec./60,DeltaPrcnt_TotSpHdSurfFluor_Expdata,'o','MarkerEdgeColor','none','MarkerFaceColor',[250,128,114]./255);
hold on;
Spine_rec_init=Cpsdf_N_save(StimSpineLoc,1)+Cpsdb_N_save(StimSpineLoc,1)+Cesm_N_save(StimSpineLoc,1)+...
    Cpsdf_A_save(StimSpineLoc,1)+Cpsdb_A_save(StimSpineLoc,1)+Cesm_A_save(StimSpineLoc,1);

Spine_rec=Cpsdf_N_save(StimSpineLoc,:)+Cpsdb_N_save(StimSpineLoc,:)+Cesm_N_save(StimSpineLoc,:)+...
    Cpsdf_A_save(StimSpineLoc,:)+Cpsdb_A_save(StimSpineLoc,:)+Cesm_A_save(StimSpineLoc,:);
plot(t_rec./60,((Spine_rec-Spine_rec_init)./Spine_rec_init)*100,'Color',[0,0,0]./255,'LineWidth',1,'LineStyle','-');

xlim([0,30]);ylim([0,250]);
xlabel('Time (min.)');ylabel('\DeltaSp. Hd. Tot. Surf. Fluor.(%)');
title({'% change in total spine head surface AMPAR fluroresence','at the stimulated spine'});
legend('Exp. data: Patterson et al. 2010','Simulated data: Present Model')