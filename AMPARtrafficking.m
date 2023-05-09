load('Initialization.mat');
load('Dist_dendExo_FoldsIncr_Expdata.mat');
%% Geometric parameters:
%Geometric dimensions of dendrite and its compartments
r_dend=0.5;%radial thickness of dendrite in um.
l_dend=0.1;%axial length of each dendritic compartment in um.
L_dend=100;%total length of dendritic branch in um.
A_dend=2*pi*r_dend*l_dend;%Area of the dendritic compartment in um^2.
N_dend=L_dend/l_dend;%No. of dendritic compartments.

%Equidistant spine locations
n=10;%choose an integer between 1 and N_dend.
dNeckLoc=n*l_dend;%the equidistance in n multiples of smallest discrete distance l_dend.
NeckLoc=0:dNeckLoc:L_dend;%Neck locations in um from the origin x=0.
NeckLoc(1)=[];%Drop the first 0um location, scheme: a dendritic compartment's distal edge location belongs to that compartment.
NeckLocIndx=NeckLoc./l_dend;%Dendritic compartment indices carrying the spines.
N_spine=length(NeckLoc);%No. of spines finally on the dendritic stretch.

%Geometric dimensions of spine necks
r_neck=0.0564*ones(N_spine,1);%radial thickness of neck in um.
l_neck=0.3*ones(N_spine,1);%neck length in um.
A_neck=2*pi*r_neck.*l_neck;%neck area in um^2.

%Geometric dimensions of spine head:ESM and PSD
A_psd=0.1*ones(N_spine,1);%area of the ESM in um^2.
A_esm=10*A_psd;%area of the PSD in um^2.

%% Model parameters
%Hopping parameters, h (um^2.s^-1)
D_dend=0.1;%Diffusion coefficient um^2.s^-1 in dend membrane.
h_dend_to_dend=(2*pi*r_dend/l_dend)*D_dend;
D_neck=6.7e-3*ones(N_spine,1);%Diffusion coefficient um^2.s^-1 in neck membrane.
D_esm=0.1*ones(N_spine,1);%Diffusion coefficient um^2.s^-1 in spine ESM membrane.
p_esm_vs_psd=1*ones(N_spine,1);%permeability at the boundaries from ESM to PSD compartments.
h_esm_to_psd=(2*pi*D_esm).*p_esm_vs_psd;
h_psd_to_esm=0.1*h_esm_to_psd;

%Desired Resm, Rdend, IRspine and IR dend: for homogeneous parameterization;
%NOTE:These desried values are for unphosphorylated naive receptors.
Resm_desrd=10/(1+gamma);
Rdend_desrd=10/(1+gamma);
IRspine_desrd=100/(1+gamma);
IRdend_desrd=10/(1+gamma);

%Exocytosis and Endocytosis parameters
Exo_dend=1.703430727016970e-04*ones(N_dend,1);%Exocyt. rate for each dendritic comp. (in s^-1) from local IR pool.: OPTIMIZED
Exo_spine=4.198050644475426e-04*ones(N_spine,1);%Exocyt. rate for each spine head (in s^-1) from local IR pool.: OPTIMIZED
Endo_dend=(IRdend_desrd/(Rdend_desrd*A_dend))*Exo_dend;%Endocyt. rate for each dendritic comp. (in s^-1) from unit surface area.: CONSTRAINED.
Endo_spine=(IRspine_desrd/(Resm_desrd*A_esm(1)))*Exo_spine;%Endocyt. rate for each spine head (in s^-1) from unit surface area.: CONSTRAINED.

%Retrograde and Anterograde IR trafficking parameters
Retrg_dend=0.028863480228992*ones(N_dend,1);%Retrograde loss rate for each dendritic comp. (in s^-1) from IR pool.: OPTIMIZED
Retrg_spine=0.071964117255061*ones(N_spine,1);%Retrograde loss rate for each spine head (in s^-1) from IR pool.: OPTIMIZED
Antrg_dend=IRdend_desrd*Retrg_dend;%Anterograde gain rate for each dendritic comp. (in s^-1) to IR pool.: CONSTRAINED.
Antrg_spine=IRspine_desrd*Retrg_spine;%Anterograde gain rate for each spine head (in s^-1) to IR pool.: CONSTRAINED.

%Binding and Unbinding parameters in PSD
CPSD95=100*ones(N_spine,1);%Count of PSD95 slots in the PSD region.
Bind_psd_N=1*ones(N_spine,1);%Binding rate, in s^-1, of free Naive Receptors in the PSD (PSDF) with PSD95.:ARBITRARILY CONSIDERED
Bind_psd_A=1*ones(N_spine,1);%Binding rate, in s^-1, of free Active Receptors in the PSD (PSDF) with PSD95.:ARBITRARILY CONSIDERED
Unbind_psd_N=2.700000003000000e+03*Bind_psd_N;%Unbinding rate, in s^-1, of bound Naive Receptors in the PSD (PSDB) from PSD95.: OPTIMIZED
Unbind_psd_A=33.333333000000000*Bind_psd_A;%Unbinding rate, in s^-1, of bound Active Receptors in the PSD (PSDB) from PSD95.: OPTIMIZED

%Dephosphorylation parameter for the entire dendritic branch
Dephos_spine_base=1e-4;%Baseline dephosp. rate (in s^-1)

%% Time and Dynamical variables
%Timescale in seconds.
Ntme_rec=(Tf*1)+1;%recording variables at every 1s
t_rec=linspace(Ti,Tf,Ntme_rec);

%Variables: Naive Receptors, together with initialization from saved baseline steady-state data
Cpsdf_N=A_psd.*Rpsdf_N_save;
Cpsdb_N=A_psd.*Rpsdb_N_save;
Cesm_N=A_esm.*Resm_N_save;
IRspine_N=IRspine_N_save;
Cneck_N=A_neck.*Rneck_N_save;
Cdend_N=A_dend*Rdend_N_save;
IRdend_N=IRdend_N_save;
%Variables: Active Receptors, initialized for being not initially present in baseline condition
Cpsdf_A=A_psd.*Rpsdf_A_save;
Cpsdb_A=A_psd.*Rpsdb_A_save;
Cesm_A=A_esm.*Resm_A_save;
IRspine_A=IRspine_A_save;
Cneck_A=A_neck.*Rneck_A_save;
Cdend_A=A_dend*Rdend_A_save;
IRdend_A=IRdend_A_save;

%% C_init vector (for ODEsuite) and Compartmental indices list
%C_init:Naive Receptors+Active Receptors 
C_init=[Cpsdf_N;Cpsdb_N;Cesm_N;IRspine_N;Cneck_N;Cdend_N;IRdend_N;...
    Cpsdf_A;Cpsdb_A;Cesm_A;IRspine_A;Cneck_A;Cdend_A;IRdend_A];
%End indices of individual comp.(s) in C_init column
Comp_indics=[N_spine;2*N_spine;3*N_spine;4*N_spine;5*N_spine;(5*N_spine)+N_dend;(5*N_spine)+(2*N_dend);...
    ((5*N_spine)+(2*N_dend))+[N_spine;2*N_spine;3*N_spine;4*N_spine;5*N_spine;(5*N_spine)+N_dend;(5*N_spine)+(2*N_dend)]];
    
%% LTP-associated perturbations at stimulated spine locations
StimDendLoc=10*StimSpineLoc;
A_esm_time_array=A_esm(1)+(A_esm(1)*DeltaPrcnt_SpHdSurfArea_ExpData./100);
Exo_spine_time_array=Exo_spine(1)*ones(1,length(t_data));
Exo_spine_time_array(1:find(t_data==50))=26*Exo_spine(1);
Exo_spine_time_array(find(t_data==50):find(t_data==(50+10)))=3*2*Exo_spine(1);
Exo_dend_time_array=Exo_dend(1)*ones(41,length(t_data));%41 is fixed, given stim. effect upto 4um, equivalent to 41 dend. comp. 
Exo_dend_time_array(:,1:find(t_data==50))=Exo_dend(1)*repmat(Dist_dendExo_FoldsIncr',1,length(1:find(t_data==50)));

%% Load enzyme profiles
kdephos=Dephos_spine_base/Edephos_spine(1,1);%Mapping from Edephos conc. to dephos. rate.

%% Time-independent Force field coefficients
%NOTE:Identical sequence and structure followed in time-dependent coefficients under loop below, so that you can easily move items between the
%time-independent and -dependent sections depending on your exact experimental requirements.

%Since diff. for Ns and As are identical and time-invariant+Areas of compartments except ESM are also time-invariant:
h_dend_to_dend=h_dend_to_dend/A_dend;
h_psd_to_esm=h_psd_to_esm./A_psd;

%Since Bind_psd_N and Bind_psd_A are time-invariant+A_psd is time-invariant:
Bind_psd_N=Bind_psd_N./A_psd;
Bind_psd_A=Bind_psd_A./A_psd;

%PSDF Pool:
 %Within Naive:
    V_psdb_N_to_psdf_N=Unbind_psd_N;
 %Within Active:
    V_psdb_A_to_psdf_A=Unbind_psd_A;
 %Between Naive and Active:

%PSDB Pool:
 %Within Naive:
 %Within Active:
 %Between Naive and Active:
    
%ESM Pool:
 %Within Naive:
    V_psdf_N_to_esm_N=h_psd_to_esm;
 %Within Active:
    V_psdf_A_to_esm_A=h_psd_to_esm;   
 %Between Naive and Active:
    
%IRspine Pool:
 %Within Naive:
    V_esm_N_to_IRspine_N=Endo_spine;
 %Within Active:
    V_esm_A_to_IRspine_A=Endo_spine;
 %Between Naive and Active:
    
%Neck Pool:
 %Within Naive:
 %Within Active:
 %Between Naive and Active:
 
%Dend Pool:
 %Within Naive:
    V_neck_N_to_dend_N=zeros(N_dend,1);
    %M_dend_to_dend_N to be defined below.
 %Within Active:  
    V_neck_A_to_dend_A=zeros(N_dend,1);
    %M_dend_to_dend_A to be defined below.
 %Between Naive and Active:
    
%IRdend Pool:
 %Within Naive:
    V_dend_N_to_IRdend_N=Endo_dend;
 %Within Active:
    V_dend_A_to_IRdend_A=Endo_dend;
 %Between Naive and Active:

%Carrier variables to be used under integration loop:
C_N=zeros(N_dend,1);
C_A=zeros(N_dend,1);

%Defining M_dend_to_dend_N:
M_dend_to_dend_N=zeros(N_dend,N_dend);
for k=1:N_dend
    if k>1 && k<N_dend
        M_dend_to_dend_N(k,k)=-2*h_dend_to_dend;
    else
        M_dend_to_dend_N(k,k)=-1*h_dend_to_dend;
    end
    if k-1>0
        M_dend_to_dend_N(k,k-1)=h_dend_to_dend;
    end
    if k+1<=N_dend
        M_dend_to_dend_N(k,k+1)=h_dend_to_dend;
    end
end
M_diag_indices=(((1:N_dend)-1)*N_dend)+(1:N_dend);

%% Boundary condition at x=0 (origin) of the dendritic stretch
%NOTE:The other end at x=L is perfectly reflective, which means no further flow or leakage from the last i=Nth compartment.
%Fixed-value boundary condition at x=0 i.e. Di=1
Cdend0N=A_dend*Rdend_N_save(1);%Receptor density,.um^-2, of Ns in D0 is fixed.
Cdend0A=A_dend*Rdend_A_save(1);%Receptor density,.um^-2, of As in D0 is fixed.
M_dend_to_dend_N(1,1)=M_dend_to_dend_N(1,1)+(-1*h_dend_to_dend);
M_dend_to_dend_A=M_dend_to_dend_N;
Inp_ext_N=zeros(N_dend,1);
Inp_ext_A=zeros(N_dend,1);
Inp_ext_N(1)=h_dend_to_dend*Cdend0N;
Inp_ext_A(1)=h_dend_to_dend*Cdend0A;

%%
Sol=ode15s(@(t,y)diffun(t,y,Comp_indics,...
    r_neck,l_neck,NeckLocIndx,A_esm,A_neck,A_dend,...
    h_psd_to_esm,h_esm_to_psd,D_neck,...
    Exo_spine,Exo_dend,Endo_spine,Endo_dend,Retrg_spine,Retrg_dend,Antrg_spine,Antrg_dend,CPSD95,Bind_psd_N,Bind_psd_A,Unbind_psd_N,Unbind_psd_A,...
    t_data,StimSpineLoc,StimDendLoc,r_neck_time_array,l_neck_time_array,D_neck_time_array,A_esm_time_array,Exo_spine_time_array,Exo_dend_time_array,...
    Ephos_PSDB,Ephos_spine,Ephos_neck,Ephos_dend,Edephos_spine,Edephos_neck,Edephos_dend,kdephos,...
    V_psdb_N_to_psdf_N,V_psdb_A_to_psdf_A,V_psdf_N_to_esm_N,V_psdf_A_to_esm_A,V_esm_N_to_IRspine_N,...
    V_esm_A_to_IRspine_A,V_neck_N_to_dend_N,V_neck_A_to_dend_A,...
    V_dend_N_to_IRdend_N,V_dend_A_to_IRdend_A,C_N,C_A,M_dend_to_dend_N,M_diag_indices,Inp_ext_N,Inp_ext_A),[Ti,Tf],C_init);

Crec=deval(Sol,t_rec);
%%
Cpsdf_N_save=Crec(1:Comp_indics(1),:);
Cpsdb_N_save=Crec((Comp_indics(1)+1):Comp_indics(2),:);
Cesm_N_save=Crec((Comp_indics(2)+1):Comp_indics(3),:);
IRSpine_N_save=Crec((Comp_indics(3)+1):Comp_indics(4),:);
Cneck_N_save=Crec((Comp_indics(4)+1):Comp_indics(5),:);
Cdend_N_save=Crec((Comp_indics(5)+1):Comp_indics(6),:);
IRDend_N_save=Crec((Comp_indics(6)+1):Comp_indics(7),:);

Cpsdf_A_save=Crec((Comp_indics(7)+1):Comp_indics(8),:);
Cpsdb_A_save=Crec((Comp_indics(8)+1):Comp_indics(9),:);
Cesm_A_save=Crec((Comp_indics(9)+1):Comp_indics(10),:);
IRSpine_A_save=Crec((Comp_indics(10)+1):Comp_indics(11),:);
Cneck_A_save=Crec((Comp_indics(11)+1):Comp_indics(12),:);
Cdend_A_save=Crec((Comp_indics(12)+1):Comp_indics(13),:);
IRDend_A_save=Crec((Comp_indics(13)+1):Comp_indics(14),:);
%%
clearvars -except StimState StimSpineLoc StimDendLoc DeltaPrcnt_SpHdSurfArea_ExpData...
    Ti Tf t_data t_rec r_neck_time_array l_neck_time_array D_neck_time_array Ephos_spine Ephos_PSDB Ephos_neck Ephos_dend...
    Edephos_spine Edephos_neck Edephos_dend...
    Cpsdf_N_save Cpsdb_N_save Cesm_N_save IRSpine_N_save Cneck_N_save Cdend_N_save IRDend_N_save...
    Cpsdf_A_save Cpsdb_A_save Cesm_A_save IRSpine_A_save Cneck_A_save Cdend_A_save IRDend_A_save;

Plot;
