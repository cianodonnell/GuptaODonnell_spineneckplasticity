%% Single Spine Stimulation: set the stimulation conditions

%STATE SPINE LOCATION:of the stimulated spine
StimSpineLoc=[50];%only integer multiples of 1um.

%NECK MORPHOLOGICAL PLASTICITY:of the stimulated spine
%No plasticity=0, 
%Towards widened and shortened neck=1
%Towards narrowed and elongated neck=2
NeckMorphPlast=0;

   %If neck morphological plasticity, state time constant tau_NM in seconds
    tauNM=300;

%NECK SEPTIN PLASTICITY:
%No plasticity:Septin 'ON'=0 , 
%Plasticity:Septin 'OFF'=1
NeckSepPlast=0;

    %If neck septin plasticity, state time constant tau_NS in seconds
    tauNS=1;


Setup;