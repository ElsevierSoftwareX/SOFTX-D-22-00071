%#########################################################################
%   SR SR kernels, case IR-CPMG or SR-CPMG.
%#########################################################################
function [Kernel_1,Kernel_2] = SR_SR_Kernel
 
 Kernel_1 = inline('1-exp(- Tau * (1./ T1))','Tau','T1');
 Kernel_2 = inline('1-exp(- Tau * (1./ T2))','Tau','T2');
end

