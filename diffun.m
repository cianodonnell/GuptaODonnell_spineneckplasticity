function dCdt=diffun(t,C,Comp_indics,...
    r_neck,l_neck,NeckLocIndx,A_esm,A_neck,A_dend,...
    h_psd_to_esm,h_esm_to_psd,D_neck,...
    Exo_spine,Exo_dend,Endo_spine,Endo_dend,Retrg_spine,Retrg_dend,Antrg_spine,Antrg_dend,CPSD95,Bind_psd_N,Bind_psd_A,Unbind_psd_N,Unbind_psd_A,...
    t_data,var_spineloc,var_dendloc,r_neck_time_array,l_neck_time_array,D_neck_time_array,A_esm_time_array,Exo_spine_time_array,Exo_dend_time_array,...
    Ephos_PSD,Ephos_spine,Ephos_neck,Ephos_dend,Edephos_spine,Edephos_neck,Edephos_dend,kdephos,...
    V_psdb_N_to_psdf_N,V_psdb_A_to_psdf_A,V_psdf_N_to_esm_N,V_psdf_A_to_esm_A,V_esm_N_to_IRspine_N,...
    V_esm_A_to_IRspine_A,V_neck_N_to_dend_N,V_neck_A_to_dend_A,...
    V_dend_N_to_IRdend_N,V_dend_A_to_IRdend_A,C_N,C_A,M_dend_to_dend_N,M_diag_indices,Inp_ext_N,Inp_ext_A)

dCdt=zeros(Comp_indics(end),1);

Cpsdf_N_initial=C(1:Comp_indics(1));
Cpsdb_N_initial=C((Comp_indics(1)+1):Comp_indics(2));
Cesm_N_initial=C((Comp_indics(2)+1):Comp_indics(3));
IRspine_N_initial=C((Comp_indics(3)+1):Comp_indics(4));
Cneck_N_initial=C((Comp_indics(4)+1):Comp_indics(5));
Cdend_N_initial=C((Comp_indics(5)+1):Comp_indics(6));
IRdend_N_initial=C((Comp_indics(6)+1):Comp_indics(7));
Cpsdf_A_initial=C((Comp_indics(7)+1):Comp_indics(8));
Cpsdb_A_initial=C((Comp_indics(8)+1):Comp_indics(9));
Cesm_A_initial=C((Comp_indics(9)+1):Comp_indics(10));
IRspine_A_initial=C((Comp_indics(10)+1):Comp_indics(11));
Cneck_A_initial=C((Comp_indics(11)+1):Comp_indics(12));
Cdend_A_initial=C((Comp_indics(12)+1):Comp_indics(13));
IRdend_A_initial=C((Comp_indics(13)+1):Comp_indics(14));

%       Time-dependent parameter variation:
        Index=find(t_data<=t,1,'last');
        r_neck(var_spineloc)=r_neck_time_array(Index);
        l_neck(var_spineloc)=l_neck_time_array(Index);
        D_neck(var_spineloc)=D_neck_time_array(Index);
        A_neck(var_spineloc)=2*pi*r_neck(var_spineloc).*l_neck(var_spineloc);
        A_esm(var_spineloc)=A_esm_time_array(Index);
        Exo_spine(var_spineloc)=Exo_spine_time_array(Index);
        Exo_dend((var_dendloc-40):var_dendloc)=flip(Exo_dend_time_array(:,Index));
        Exo_dend(var_dendloc:(var_dendloc+40))=Exo_dend_time_array(:,Index);
     
        
        Phos_spine_bound=Ephos_PSD(:,Index);
        Phos_spine=Ephos_spine(:,Index);
        Phos_neck=Ephos_neck(:,Index);
        Phos_dend=Ephos_dend(:,Index);
        Dephos_spine=kdephos*Edephos_spine(:,Index);
        Dephos_neck=kdephos*Edephos_neck(:,Index);
        Dephos_dend=kdephos*Edephos_dend(:,Index);
        
     
        % Time-dependent Force field coefficients: due to possibly one or 
        % more time-dependent parameters' variations and state-dependence
        
        %Since, A_esm is time-variant.
        h_esm_to_psd_var=h_esm_to_psd./A_esm;
        
        %Since neck morphology and septin gating is time-variant
        h_dend_to_neck=((2*pi*r_neck)./(0.5*l_neck)).*D_neck;
        h_neck_to_dend=h_dend_to_neck;
        h_neck_to_esm=h_dend_to_neck;
        h_esm_to_neck=h_neck_to_esm;
        
        h_dend_to_neck=h_dend_to_neck./A_dend;
        h_neck_to_dend=h_neck_to_dend./A_neck;
        h_neck_to_esm=h_neck_to_esm./A_neck;
        h_esm_to_neck=h_esm_to_neck./A_esm;
        
        
        %PSDF Pool:
         %Within Naive:
            Binding_N=Bind_psd_N.*(CPSD95-(Cpsdb_N_initial+Cpsdb_A_initial));
            V_psdf_N_to_psdf_N=(-1*h_psd_to_esm)+(-1*Binding_N)+(-1*Phos_spine);
            V_esm_N_to_psdf_N=h_esm_to_psd_var;
         %Within Active:
            Binding_A=Bind_psd_A.*(CPSD95-(Cpsdb_N_initial+Cpsdb_A_initial));
            V_psdf_A_to_psdf_A=(-1*h_psd_to_esm)+(-1*Binding_A)+(-1*Dephos_spine);
            V_esm_A_to_psdf_A=h_esm_to_psd_var;
         %Between Naive and Active:
            V_psdf_N_to_psdf_A=Phos_spine;
            V_psdf_A_to_psdf_N=Dephos_spine;

        %PSDB Pool:
         %Within Naive:
            V_psdf_N_to_psdb_N=Binding_N;
            V_psdb_N_to_psdb_N=(-1*Unbind_psd_N)+(-1*Phos_spine_bound);
         %Within Active:
            V_psdf_A_to_psdb_A=Binding_A;
            V_psdb_A_to_psdb_A=(-1*Unbind_psd_A)+(-1*Dephos_spine); 
         %Between Naive and Active:
            V_psdb_N_to_psdb_A=Phos_spine_bound;
            V_psdb_A_to_psdb_N=Dephos_spine;

        %ESM Pool:
         %Within Naive:
            V_esm_N_to_esm_N=(-1*h_esm_to_psd_var)+(-1*h_esm_to_neck)+(-1*Endo_spine)+(-1*Phos_spine);
            V_neck_N_to_esm_N=h_neck_to_esm;
            V_IRspine_N_to_esm_N=Exo_spine;
         %Within Active:  
            V_esm_A_to_esm_A=(-1*h_esm_to_psd_var)+(-1*h_esm_to_neck)+(-1*Endo_spine)+(-1*Dephos_spine);
            V_neck_A_to_esm_A=h_neck_to_esm;
            V_IRspine_A_to_esm_A=Exo_spine;
         %Between Naive and Active:
            V_esm_N_to_esm_A=Phos_spine;
            V_esm_A_to_esm_N=Dephos_spine;

        %IRspine Pool:
         %Within Naive: 
            V_IRspine_N_to_IRspine_N=(-1*Exo_spine)+(-1*Retrg_spine)+(-1*Phos_spine);
         %Within Active: 
            V_IRspine_A_to_IRspine_A=(-1*Exo_spine)+(-1*Dephos_spine);
         %Between Naive and Active:
            V_IRspine_N_to_IRspine_A=Phos_spine;
            V_IRspine_A_to_IRspine_N=Dephos_spine;

        %Neck Pool:
         %Within Naive: 
            V_esm_N_to_neck_N=h_esm_to_neck;
            V_neck_N_to_neck_N=(-1*h_neck_to_esm)+(-1*h_neck_to_dend)+(-1*Phos_neck);
            V_dend_N_to_neck_N=h_dend_to_neck;
         %Within Active: 
            V_esm_A_to_neck_A=h_esm_to_neck;
            V_neck_A_to_neck_A=(-1*h_neck_to_esm)+(-1*h_neck_to_dend)+(-1*Dephos_neck);
            V_dend_A_to_neck_A=h_dend_to_neck;
         %Between Naive and Active:
            V_neck_N_to_neck_A=Phos_neck;
            V_neck_A_to_neck_N=Dephos_neck; 
         
            
        %Dend Pool:
         %Within Naive: 
            V_neck_N_to_dend_N(NeckLocIndx)=h_neck_to_dend;
            %M_dend_to_dend_N already defined.Endocytosis and Phos./Dephos. will be implemented seperately in
            %integration, to save futile restoring-computation at M_diag_indices due to these reaction terms.
            M_dend_to_dend_N(M_diag_indices(NeckLocIndx))=M_dend_to_dend_N(M_diag_indices(NeckLocIndx))+(-1*h_dend_to_neck');
            V_IRdend_N_to_dend_N=Exo_dend;
         %Within Active:
            V_neck_A_to_dend_A(NeckLocIndx)=h_neck_to_dend;
            %M_dend_to_dend_A already defined.Endocytosis and Phos./Dephos. will be implemented seperately in
            %integration, to save futile restoring-computation at M_diag_indices due to these reaction terms.
            M_dend_to_dend_A=M_dend_to_dend_N;
            V_IRdend_A_to_dend_A=Exo_dend;
         %Between Naive and Active:
            V_dend_N_to_dend_A=Phos_dend;
            V_dend_A_to_dend_N=Dephos_dend;


        %IRdend Pool:
         %Within Naive:
            V_IRdend_N_to_IRdend_N=(-1*Exo_dend)+(-1*Retrg_dend)+(-1*Phos_dend);
         %Within Active:  
            V_IRdend_A_to_IRdend_A=(-1*Exo_dend)+(-1*Dephos_dend);
         %Between Naive and Active:
            V_IRdend_N_to_IRdend_A=Phos_dend;
            V_IRdend_A_to_IRdend_N=Dephos_dend;
                
             
        %Integration:Naive Receptors(first within-pool terms + lastly cross-pool term)        
        dCdt(1:Comp_indics(1))=(V_psdf_N_to_psdf_N.*Cpsdf_N_initial)+(V_psdb_N_to_psdf_N.*Cpsdb_N_initial)+(V_esm_N_to_psdf_N.*Cesm_N_initial)+(V_psdf_A_to_psdf_N.*Cpsdf_A_initial);
        dCdt((Comp_indics(1)+1):Comp_indics(2))=(V_psdf_N_to_psdb_N.*Cpsdf_N_initial)+(V_psdb_N_to_psdb_N.*Cpsdb_N_initial)+(V_psdb_A_to_psdb_N.*Cpsdb_A_initial);
        dCdt((Comp_indics(2)+1):Comp_indics(3))=(V_psdf_N_to_esm_N.*Cpsdf_N_initial)+(V_esm_N_to_esm_N.*Cesm_N_initial)+(V_IRspine_N_to_esm_N.*IRspine_N_initial)+(V_neck_N_to_esm_N.*Cneck_N_initial)+(V_esm_A_to_esm_N.*Cesm_A_initial);
        dCdt((Comp_indics(3)+1):Comp_indics(4))=(V_esm_N_to_IRspine_N.*Cesm_N_initial)+(V_IRspine_N_to_IRspine_N.*IRspine_N_initial)+Antrg_spine+(V_IRspine_A_to_IRspine_N.*IRspine_A_initial);
        dCdt((Comp_indics(4)+1):Comp_indics(5))=(V_esm_N_to_neck_N.*Cesm_N_initial)+(V_neck_N_to_neck_N.*Cneck_N_initial)+(V_dend_N_to_neck_N.*Cdend_N_initial(NeckLocIndx))+(V_neck_A_to_neck_N.*Cneck_A_initial);
        C_N(NeckLocIndx)=Cneck_N_initial;
        dCdt((Comp_indics(5)+1):Comp_indics(6))=(V_neck_N_to_dend_N.*C_N)+(M_dend_to_dend_N*Cdend_N_initial)+(-1*Endo_dend.*Cdend_N_initial)+(-1*Phos_dend.*Cdend_N_initial)+(V_IRdend_N_to_dend_N.*IRdend_N_initial)+Inp_ext_N+(V_dend_A_to_dend_N.*Cdend_A_initial);
        dCdt((Comp_indics(6)+1):Comp_indics(7))=(V_dend_N_to_IRdend_N.*Cdend_N_initial)+(V_IRdend_N_to_IRdend_N.*IRdend_N_initial)+Antrg_dend+(V_IRdend_A_to_IRdend_N.*IRdend_A_initial);
        
        %Integration:Active Receptors(first within-pool terms + lastly cross-pool term)        
        dCdt((Comp_indics(7)+1):Comp_indics(8))=(V_psdf_A_to_psdf_A.*Cpsdf_A_initial)+(V_psdb_A_to_psdf_A.*Cpsdb_A_initial)+(V_esm_A_to_psdf_A.*Cesm_A_initial)+(V_psdf_N_to_psdf_A.*Cpsdf_N_initial);
        dCdt((Comp_indics(8)+1):Comp_indics(9))=(V_psdf_A_to_psdb_A.*Cpsdf_A_initial)+(V_psdb_A_to_psdb_A.*Cpsdb_A_initial)+(V_psdb_N_to_psdb_A.*Cpsdb_N_initial);
        dCdt((Comp_indics(9)+1):Comp_indics(10))=(V_psdf_A_to_esm_A.*Cpsdf_A_initial)+(V_esm_A_to_esm_A.*Cesm_A_initial)+(V_IRspine_A_to_esm_A.*IRspine_A_initial)+(V_neck_A_to_esm_A.*Cneck_A_initial)+(V_esm_N_to_esm_A.*Cesm_N_initial);
        dCdt((Comp_indics(10)+1):Comp_indics(11))=(V_esm_A_to_IRspine_A.*Cesm_A_initial)+(V_IRspine_A_to_IRspine_A.*IRspine_A_initial)+(V_IRspine_N_to_IRspine_A.*IRspine_N_initial);
        dCdt((Comp_indics(11)+1):Comp_indics(12))=(V_esm_A_to_neck_A.*Cesm_A_initial)+(V_neck_A_to_neck_A.*Cneck_A_initial)+(V_dend_A_to_neck_A.*Cdend_A_initial(NeckLocIndx))+(V_neck_N_to_neck_A.*Cneck_N_initial);
        C_A(NeckLocIndx)=Cneck_A_initial;
        dCdt((Comp_indics(12)+1):Comp_indics(13))=(V_neck_A_to_dend_A.*C_A)+(M_dend_to_dend_A*Cdend_A_initial)+(-1*Endo_dend.*Cdend_A_initial)+(-1*Dephos_dend.*Cdend_A_initial)+(V_IRdend_A_to_dend_A.*IRdend_A_initial)+Inp_ext_A+(V_dend_N_to_dend_A.*Cdend_N_initial);
        dCdt((Comp_indics(13)+1):Comp_indics(14))=(V_dend_A_to_IRdend_A.*Cdend_A_initial)+(V_IRdend_A_to_IRdend_A.*IRdend_A_initial)+(V_IRdend_N_to_IRdend_A.*IRdend_N_initial);
        
        M_dend_to_dend_N(M_diag_indices(NeckLocIndx))=M_dend_to_dend_N(M_diag_indices(NeckLocIndx))+h_dend_to_neck';
end