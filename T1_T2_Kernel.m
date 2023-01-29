%#########################################################################
%   T1-T2 kernels, case IR-CPMG or SR-CPMG.
%#########################################################################
function [Kernel_1,Kernel_2] = T1_T2_Kernel(FL_typeKernel)
 if (FL_typeKernel==1)
    Kernel_1 = inline('1-2*exp(- Tau * (1./ T1))','Tau','T1');
  elseif(FL_typeKernel==2)
    Kernel_1 = inline('1-exp(- Tau * (1./ T1))','Tau','T1');
 end
 Kernel_2 = inline('exp(- Tau * (1./ T2))','Tau','T2');
end

