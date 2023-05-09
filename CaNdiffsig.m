function dCdt=CaNdiffsig(t,C,t_rec,Comp_indics,NeckLocIndx,Twin,r_neck,l_neck,V_neck,V_spine,V_dend,varspineloc,tau,D_neck,...
    Edephos_spine_base,Edephos_neck_base,Edephos_dend_base,r_neck_time_array,l_neck_time_array,...
    h_spine_vs_neck,V_spine_array,V_neck_to_dend,C_Dephos,M_diag_indices,M_dend_to_dend)

dCdt=zeros(Comp_indics(end),1);
Edephos_spine=C(1:Comp_indics(1));
Edephos_neck=C((Comp_indics(1)+1):Comp_indics(2));
Edephos_dend=C((Comp_indics(2)+1):Comp_indics(3));

Index=find(t_rec<=t,1,'last');
V_spine(varspineloc)=V_spine_array(Index);
r_neck(varspineloc)=r_neck_time_array(Index);
l_neck(varspineloc)=l_neck_time_array(Index);
V_neck(varspineloc)=pi*(r_neck(varspineloc).^2).*l_neck(varspineloc);

h_spine_vs_neck(varspineloc)=((pi*(r_neck(varspineloc).^2)*D_neck)./(0.5*l_neck(varspineloc)));
h_neck_vs_dend=h_spine_vs_neck;

V_spine_to_spine=(-1/tau)+(-1*h_spine_vs_neck./V_spine);
V_neck_to_spine=h_spine_vs_neck./V_spine;

V_spine_to_neck=h_spine_vs_neck./V_neck;
V_neck_to_neck=(-1/tau)+(-1*h_spine_vs_neck./V_neck)+(-1*h_neck_vs_dend./V_neck);
V_dend_to_neck=h_neck_vs_dend./V_neck;

V_neck_to_dend(NeckLocIndx)=h_neck_vs_dend./V_dend;
M_dend_to_dend(M_diag_indices(NeckLocIndx))=M_dend_to_dend(M_diag_indices(NeckLocIndx))+(-1*h_neck_vs_dend'./V_dend);

   %Edephos Dynamics: 
   dCdt(1:Comp_indics(1))=(V_spine_to_spine.*(Edephos_spine-Edephos_spine_base))+(V_neck_to_spine.*(Edephos_neck-Edephos_neck_base));
   dCdt((Comp_indics(1)+1):Comp_indics(2))=(V_spine_to_neck.*(Edephos_spine-Edephos_spine_base))+(V_neck_to_neck.*(Edephos_neck-Edephos_neck_base))+(V_dend_to_neck.*(Edephos_dend(NeckLocIndx)-Edephos_dend_base));
   C_Dephos(NeckLocIndx)=Edephos_neck-Edephos_neck_base;
   dCdt((Comp_indics(2)+1):Comp_indics(3))=(V_neck_to_dend.*C_Dephos)+(M_dend_to_dend*(Edephos_dend-Edephos_dend_base));
   
   if t<=Twin
      dCdt(varspineloc)=0;
   end
   
   M_dend_to_dend(M_diag_indices(NeckLocIndx))=M_dend_to_dend(M_diag_indices(NeckLocIndx))+(h_neck_vs_dend'./V_dend);
   
end