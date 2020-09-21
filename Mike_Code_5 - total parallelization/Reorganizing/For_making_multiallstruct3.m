%multiallstruct = struct;
%OverStruct = struct;

%current = 0;
%n = current;

%AllStruct = OverStruct.AllStruct
la = length(OverStruct.list)

n = length(multiallstruct)

for i = 1:la

	multiallstruct(n + i).name      = OverStruct.list(i).name     ;
    multiallstruct(n + i).T         = OverStruct.list(i).T        ;
    multiallstruct(n + i).N         = OverStruct.list(i).N        ;
    multiallstruct(n + i).beta      = OverStruct.list(i).beta     ;
    multiallstruct(n + i).sprsvec   = OverStruct.list(i).sprsvec  ;
    multiallstruct(n + i).topology  = OverStruct.list(i).topology ;
    multiallstruct(n + i).Jcontru   = OverStruct.list(i).Jcontru  ;
    multiallstruct(n + i).Jtru      = OverStruct.list(i).Jtru     ;
    multiallstruct(n + i).htru      = OverStruct.list(i).htru     ;
    multiallstruct(n + i).couplings = OverStruct.list(i).couplings;
    multiallstruct(n + i).c         = OverStruct.list(i).c        ;
    %multiallstruct(n + i).S         = OverStruct.list(i).S        ;
    multiallstruct(n + i).BLLH      = OverStruct.list(i).BLLH     ;

end

