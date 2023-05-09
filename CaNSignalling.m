%% Geometric parameters:
%Geometric dimensions of dendrite and its compartments
r_dend=0.5;%radial thickness of dendrite in um.
l_dend=0.1;%axial length of each dendritic compartment in um.
L_dend=100;%total length of dendritic branch in um.
V_dend=pi*(r_dend^2)*l_dend;%volume of the dendritic compartment in um^3.
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
V_neck=pi*(r_neck.^2).*l_neck;%neck volume in um^3.

%Geometric dimensions of spine head:ESM and PSD
A_psd=0.1*ones(N_spine,1);%area of the ESM in um^2.
A_esm=10*A_psd;%area of the PSD in um^2.
V_spine=(4/3)*pi*sqrt((A_psd+A_esm)./(4*pi)).^3;%volume of the spine head in um^3.

%% Diffusion coefficients
D_neck=0.106;%um^2.s^-1.
D_dend=0.106;

h_spine_vs_neck=((pi*(r_neck.^2)*D_neck)./(0.5*l_neck));
h_dend_vs_dend=(pi*(r_dend^2)*D_dend)/l_dend;

%% Dynamical Variables
Edephos_spine_base=0.1;%uM
Edephos_neck_base=0.1;%uM
Edephos_dend_base=0.1;%uM

Cinit=[Edephos_spine_base*ones(N_spine,1);Edephos_neck_base*ones(N_spine,1);Edephos_dend_base*ones(N_dend,1)];
Comp_indics=[N_spine;2*N_spine;(2*N_spine)+N_dend];

Cinit(StimSpineLoc)=300*Cinit(StimSpineLoc);%Box-Changing Edephos_spine at StimSpineLoc

%% Load Experim. data on spine head surface area
load('DeltaPrcnt_SpHdSurfArea_ExpData.mat');
TotSpHdSurf_timearray=A_psd(1)+(A_esm(1)+(A_esm(1)*DeltaPrcnt_SpHdSurfArea_ExpData./100));
Vspine_timearray=(4/3)*pi*((sqrt(TotSpHdSurf_timearray./(4*pi))).^3);%Assuming spherical geometry

%% Force-Field Coefficients
V_neck_to_dend=zeros(N_dend,1);
C_Dephos=zeros(N_dend,1);
M_dend_to_dend=zeros(N_dend,N_dend);
for k=1:N_dend
    if k>1 && k<N_dend
        M_dend_to_dend(k,k)=(-2*h_dend_vs_dend/V_dend)+(-1/tauCaN);
    else
        M_dend_to_dend(k,k)=(-1*h_dend_vs_dend/V_dend)+(-1/tauCaN);
    end
    if k-1>0
        M_dend_to_dend(k,k-1)=h_dend_vs_dend/V_dend;
    end
    if k+1<=N_dend
        M_dend_to_dend(k,k+1)=h_dend_vs_dend/V_dend;
    end
end
M_diag_indices=(((1:N_dend)-1)*N_dend)+(1:N_dend);
M_dend_to_dend(1,1)=M_dend_to_dend(1,1)+(-1*h_dend_vs_dend/V_dend);

%% Solver-based integration
Sol=ode15s(@(t,y)CaNdiffsig(t,y,t_data,Comp_indics,NeckLocIndx,Twin,r_neck,l_neck,V_neck,V_spine,V_dend,StimSpineLoc,tauCaN,D_neck,...
    Edephos_spine_base,Edephos_neck_base,Edephos_dend_base,r_neck_time_array,l_neck_time_array,...
    h_spine_vs_neck,Vspine_timearray,V_neck_to_dend,C_Dephos,M_diag_indices,M_dend_to_dend),[Ti,Tf],Cinit);

Crec=deval(Sol,t_data);

%% Save Variables
Edephos_spine=Crec(1:Comp_indics(1),:);
Edephos_neck=Crec((Comp_indics(1)+1):Comp_indics(2),:);
Edephos_dend=Crec((Comp_indics(2)+1):Comp_indics(3),:);

clearvars -except StimSpineLoc DeltaPrcnt_SpHdSurfArea_ExpData...
    Ti Tf t_data r_neck_time_array l_neck_time_array D_neck_time_array Ephos_spine Ephos_PSDB Ephos_neck Ephos_dend...
    Edephos_spine Edephos_neck Edephos_dend;

AMPARtrafficking;
