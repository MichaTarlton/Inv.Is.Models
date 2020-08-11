multiallstruct().Nvec         = []
multiallstruct().Tvec         = []
multiallstruct().betavec      = []
multiallstruct().topologies   = []
multiallstruct().jn           = []
multiallstruct().sparsity     = []
multiallstruct().h_on         = []
multiallstruct().couplings    = []
multiallstruct().c            = []
multiallstruct().time         = []
multiallstruct().AllStruct    = AllStruct


n=2;


multiallstruct(n).AllStruct    = AllStruct;
multiallstruct(n).time         = ['2257-28Jun20'];
multiallstruct(n).betavec      = [0.4];
multiallstruct(n).sparsity     = [0];
multiallstruct(n).Nvec         = [30,60,120];
multiallstruct(n).Tvec         = [1E3];
multiallstruct(n).topologies   = {'sk',5};
multiallstruct(n).jn           = [1];
multiallstruct(n).h_on         = [1];
multiallstruct(n).couplings    = [2];
multiallstruct(n).c            = [2];
n = n + 1


multiallstruct(n).AllStruct    = OverStruct(1).AllStruct 
multiallstruct(n).time         = OverStruct(1).time      
multiallstruct(n).betavec      = OverStruct(1).betavec   
multiallstruct(n).sparsity     = OverStruct(1).sparsity(1)  
multiallstruct(n).Nvec         = OverStruct(1).Nvec      
multiallstruct(n).Tvec         = OverStruct(1).Tvec      
multiallstruct(n).topologies   = OverStruct(1).topologies
multiallstruct(n).jn           = OverStruct(1).jn        
multiallstruct(n).h_on         = OverStruct(1).h_on   
multiallstruct(n).Jstor         = OverStruct(1).Jstor

n = n + 1


Trls_1 - Beta[0.2] - Topo[sk,5] - c[2] - N[30] - T[30]



multiallstruct(n).AllStruct    = OverStruct(1).AllStruct 
multiallstruct(n).time         = OverStruct(1).time      
multiallstruct(n).betavec      = OverStruct(1).betavec   
multiallstruct(n).sparsity     = OverStruct(1).sparsity(1)  
multiallstruct(n).Nvec         = OverStruct(1).Nvec      
multiallstruct(n).Tvec         = OverStruct(1).Tvec      
multiallstruct(n).topologies   = OverStruct(1).topologies
multiallstruct(n).jn           = OverStruct(1).jn        
multiallstruct(n).h_on         = OverStruct(1).h_on   
multiallstruct(n).Jstor        = OverStruct(1).Jstor












%multiallstruct.topdir       = []