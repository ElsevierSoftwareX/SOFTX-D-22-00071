%###################################################################################################
%PURPOSE :Loads flag parameters from file.
%DATE    :15/07/2017
%VERSION :1.1 [00/00/0000]
%         1.2 [19/10/2017] (vb) Added FLNoContour flag to plot with or without contour.
%         1.3 [03/05/2019] (vb) Added FL_Verbose flag.
%         1.4 [17/07/2020] (vb) Added FL_Debug FL_Debug, FL_Amp_scale and FL_Scale_fact flags.
%                               Now FL_NoContour contains the number of isoline, set 0 for no isoline.
%         1.5 [20/07/2020] (vb) Added new parameters read from filepar.
%         1.6 [03/07/2021] (vb) Now it is used a structured data for flag parameters.
%         1.7 [19/08/2021] (vb) Now it is used "FL_Stat" for statistic purpose.
%         1.8 [19/08/2021] (vb, fz) Added the parFL.FL_Lower_bound flag. Now it is used a local structured variable.
%                                   Set all parametrs to their default value to avoid non defined parameters
%                                   if not defined in the input file.
%NOTES   :
%
%###################################################################################################
%
function [parFLa]= LoadFlags(InputFileName, UseDefault)
 %
 %
 if(UseDefault)
     parFLa.CommentTS='Default Flags';
     parFLa.FL_typeKernel=1;          %1 IR-CPMG; 2 SR-CPMG; 4 T2-T2
     parFLa.FL_InversionTimeLimits=0; %1 automatic, 0 manually selection inversion times
     parFLa.FL_OutputData=0;          %1 create output data file for ILT2D
     parFLa.FL_NoContour=0;           %
	 parFLa.FL_Verbose=0;  
     parFLa.FL_Debug=0;
     parFLa.FL_Amp_scale=1.0E0;
     parFLa.FL_Lower_bound=0.0E0;
     parFLa.FL_Scale_fact=1.0E0;
     parFLa.FL_T1T2Filter=0;
     parFLa.FL_EraseCol=0;
     parFLa.FL_EraseRow=0;
     parFLa.FL_Offset=0;
     parFLa.FL_Stat=0;
 else
     %set all the necessary parametrs. Imput file can be partial.
     parFLa.CommentTS='Default Flags';
     parFLa.FL_typeKernel=1;          
     parFLa.FL_InversionTimeLimits=0; 
     parFLa.FL_OutputData=0;          
     parFLa.FL_NoContour=0;           
	 parFLa.FL_Verbose=0;  
     parFLa.FL_Debug=0;
     parFLa.FL_Amp_scale=1.0E0;
     parFLa.FL_Lower_bound=0.0E0;
     parFLa.FL_Scale_fact=1.0E0;
     parFLa.FL_T1T2Filter=0;
     parFLa.FL_EraseCol=0;
     parFLa.FL_EraseRow=0;
     parFLa.FL_Offset=0;
     parFLa.FL_Stat=0;
     %
     fid = fopen(InputFileName);  %
     parFLa.CommentTS = fgetl(fid);       % a row of comment
     %extract flags values
     while(1)
       stringa=fgetl(fid);
       stringa=strtrim(stringa);
       if(strfind(stringa, 'END')==1) break; end   %stops reading parameters.
       if(strfind(stringa, 'FL_typeKernel         =')==1) 
          parFLa.FL_typeKernel=str2double(strrep(stringa,'FL_typeKernel         =',''));
       end
       if(strfind(stringa, 'FL_InversionTimeLimits=')==1) 
          parFLa.FL_InversionTimeLimits=str2double(strrep(stringa,'FL_InversionTimeLimits=',''));
       end
       if(strfind(stringa, 'FL_OutputData         =')==1) 
          parFLa.FL_OutputData=str2double(strrep(stringa,'FL_OutputData         =',''));
       end
       if(strfind(stringa, 'FL_NoContour          =')==1) 
          parFLa.FL_NoContour=str2double(strrep(stringa,'FL_NoContour          =',''));
       end
       if(strfind(stringa, 'FL_Verbose            =')==1) 
          parFLa.FL_Verbose=str2double(strrep(stringa,'FL_Verbose            =',''));
       end
       if(strfind(stringa, 'FL_Debug              =')==1) 
          parFLa.FL_Debug=str2double(strrep(stringa,'FL_Debug              =',''));
       end
       if(strfind(stringa, 'FL_Amp_scale          =')==1) 
          parFLa.FL_Amp_scale=str2double(strrep(stringa,'FL_Amp_scale          =',''));
       end
       if(strfind(stringa, 'FL_Lower_bound        =')==1)
           parFLa.FL_Lower_bound=str2double(strrep(stringa,'FL_Lower_bound        =',''));
       end
       if(strfind(stringa, 'FL_Scale_fact         =')==1) 
          parFLa.FL_Scale_fact=str2double(strrep(stringa,'FL_Scale_fact         =',''));
       end
       if(strfind(stringa, 'FL_T1T2Filter         =')==1) 
          parFLa.FL_T1T2Filter=str2double(strrep(stringa,'FL_T1T2Filter         =',''));
       end
       if(strfind(stringa, 'FL_EraseCol           =')==1) 
          parFLa.FL_EraseCol=str2double(strrep(stringa,'FL_EraseCol           =',''));
       end
       if(strfind(stringa, 'FL_EraseRow           =')==1) 
          parFLa.FL_EraseRow=str2double(strrep(stringa,'FL_EraseRow           =',''));
       end
       if(strfind(stringa, 'FL_Offset             =')==1) 
          parFLa.FL_Offset=str2double(strrep(stringa,'FL_Offset             =',''));
       end
       if(strfind(stringa, 'FL_Stat               =')==1) 
          parFLa.FL_Stat=str2double(strrep(stringa,'FL_Stat               =',''));
       end
     end
     fclose(fid);
 end
 return;
end

