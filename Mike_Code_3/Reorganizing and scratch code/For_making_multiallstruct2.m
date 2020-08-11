multiallstruct = struct;

current = 0;
n = current;

%AllStruct = OverStruct.AllStruct
la = length(AllStruct.list)

n = length(multiallstruct)

for i = 1:la

	multiallstruct(n + i).name      = AllStruct.list(i).name     ;
    multiallstruct(n + i).T         = AllStruct.list(i).T        ;
    multiallstruct(n + i).N         = AllStruct.list(i).N        ;
    multiallstruct(n + i).beta      = AllStruct.list(i).beta     ;
    multiallstruct(n + i).sprsvec   = AllStruct.list(i).sprsvec  ;
    multiallstruct(n + i).topology  = AllStruct.list(i).topology ;
    multiallstruct(n + i).Jcontru   = AllStruct.list(i).Jcontru  ;
    multiallstruct(n + i).Jtru      = AllStruct.list(i).Jtru     ;
    multiallstruct(n + i).htru      = AllStruct.list(i).htru     ;
    multiallstruct(n + i).couplings = AllStruct.list(i).couplings;
    multiallstruct(n + i).c         = AllStruct.list(i).c        ;
    multiallstruct(n + i).S         = AllStruct.list(i).S        ;
    multiallstruct(n + i).BLLH      = AllStruct.list(i).BLLH     ;

end

