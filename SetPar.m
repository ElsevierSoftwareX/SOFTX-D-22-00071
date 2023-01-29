%###################################################################################################
%NAME    :SetPar.m
%PURPOSE :Set/load parameters from file.
%DATE    :12/07/2017
%VERSION :1.1 [15/07/2017] manages file of parameters
%IMPUT   :InputFileName = name of the file with parameters
%NOTES   :
%###################################################################################################
function [par]= SetPar(InputFileName, s, B, par, UseDefault)
  if UseDefault
    %Set Defaulf
    % [GRADIENTE PROIETTATO]
    par.gpnr.tol          =0.1;
    par.gpnr.maxiter      =50000;
    %
    % [NEWTON PROIETTATO]
    par.nwtp.maxiter      =1E5;
    par.nwtp.tolrho       =1.0E-8;
    %
    % [CG]
    par.cgn2d.tol         =1.E-3;
    par.cgn2d.maxiter     =3000;
    %
    %[SVD]
    par.svd.svd           =1; % per proiezione-filtro SVD
    par.svd.soglia        =1.0E-6; % SOLO per par.svd=1;
    % 
    %[UPEN]
    par.upen.tol          =0.001;
    par.upen.iter         =100;
    par.upen.tol_res      =1.E-7;
    par.upen.beta00       =1;
    par.upen.beta0        =1.0E-4; 
    par.upen.beta_p       =1.0E-5;
    par.upen.beta_c       =1.0E-1;
    %[TIKHONOV]
    par.tikh.lambda       =5.0E-1;
    %
    %[FISTA]
    par.fista.maxiter     =5.0E5;
    par.fista.tol         =1.0E-7;
    par.fista.weight      =1.0;
    %par.fista.crit        =1;
    %par.fista.Thresh      =1.0E-9;
    %par.fista.tau         =1.0E0;
    %
    %Load from file if file exists.
  else 
    fid = fopen(InputFileName);  %
    %
    CommentTS = fgetl(fid);       % a row of comment
    %extract parameters
    while(1)
      stringa=fgetl(fid);
      stringa=strtrim(stringa);
      if(strfind(stringa, 'END')==1) break; end   %stops reading parameters.
      % [GRADIENTE PROIETTATO]
      if(strfind(stringa, 'par.gpnr.tol          =')==1) 
         par.gpnr.tol=str2double(strrep(stringa,'par.gpnr.tol          =',''));
      end
      if(strfind(stringa, 'par.gpnr.maxiter      =')==1) 
          par.gpnr.maxiter=str2double(strrep(stringa,'par.gpnr.maxiter      =',''));
      end
      % [NEWTON PROIETTATO]
      if(strfind(stringa, 'par.nwtp.maxiter      =')==1) 
          par.nwtp.maxiter=str2double(strrep(stringa,'par.nwtp.maxiter      =',''));
      end
      if(strfind(stringa, 'par.nwtp.tolrho       =')==1) 
          par.nwtp.tolrho=str2double(strrep(stringa,'par.nwtp.tolrho       =',''));
      end
      % [CG]
      if(strfind(stringa, 'par.cgn2d.tol         =')==1) 
          par.cgn2d.tol=str2double(strrep(stringa,'par.cgn2d.tol         =',''));
      end
      if(strfind(stringa, 'par.cgn2d.maxiter     =')==1) 
          par.cgn2d.maxiter=str2double(strrep(stringa,'par.cgn2d.maxiter     =',''));
      end
      %[SVD]
      if(strfind(stringa, 'par.svd.svd           =')==1) 
          par.svd.svd=str2double(strrep(stringa,'par.svd.svd           =',''));
      end
      if(strfind(stringa, 'par.svd.soglia        =')==1) 
          par.svd.soglia=str2double(strrep(stringa,'par.svd.soglia        =',''));
      end
      %[UPEN]
      if(strfind(stringa, 'par.upen.tol          =')==1) 
          par.upen.tol=str2double(strrep(stringa,  'par.upen.tol          =',''));
      end    
      if(strfind(stringa, 'par.upen.iter         =')==1) 
          par.upen.iter=str2double(strrep(stringa, 'par.upen.iter         =',''));
      end
      if(strfind(stringa, 'par.upen.tol_res      =')==1) 
          par.upen.tol_res=str2double(strrep(stringa, 'par.upen.tol_res      =',''));
      end
      if(strfind(stringa, 'par.upen.beta00       =')==1)
         par.upen.beta00=str2double(strrep(stringa,'par.upen.beta00       =',''));
      end
      if(strfind(stringa, 'par.upen.beta0        =')==1)
         par.upen.beta0=str2double(strrep(stringa,'par.upen.beta0        =',''));
      end
      if(strfind(stringa, 'par.upen.beta_p       =')==1)
          par.upen.beta_p=str2double(strrep(stringa,'par.upen.beta_p       =',''));
      end
      if(strfind(stringa, 'par.upen.beta_c       =')==1) 
          par.upen.beta_c=str2double(strrep(stringa,'par.upen.beta_c       =',''));
      end
	  %[TIKHONOV]
      if(strfind(stringa, 'par.tikh.lambda       =')==1) 
          par.tikh.lambda=str2double(strrep(stringa,'par.tikh.lambda       =',''));
      end
      %[FISTA]
      if(strfind(stringa, 'par.fista.maxiter     =')==1) 
          par.fista.maxiter=str2double(strrep(stringa,'par.fista.maxiter     =',''));
      end
      if(strfind(stringa, 'par.fista.tol         =')==1) 
          par.fista.tol=str2double(strrep(stringa,'par.fista.tol         =',''));
      end
      if(strfind(stringa, 'par.fista.weight      =')==1) 
          par.fista.weight=str2double(strrep(stringa,'par.fista.weight      =',''));
      end
      
%        if(strfind(stringa, 'par.fista.crit        =')==1) 
%            %par.fista.crit=str2double(strrep(stringa,'par.fista.crit        =',''));
%            crit=str2double(strrep(stringa,'par.fista.crit        =',''));
%        end
%        if(strfind(stringa, 'par.fista.Thresh      =')==1) 
%            %par.fista.Thresh=str2double(strrep(stringa,'par.fista.Thresh      =',''));
%            Thresh=str2double(strrep(stringa,'par.fista.Thresh      =',''));
%        end
%        if(strfind(stringa, 'par.fista.tau         =')==1) 
%            %par.fista.tau=str2double(strrep(stringa,'par.fista.tau         =',''));
%            tau=str2double(strrep(stringa,'par.fista.tau         =',''));
%        end
     end
    fclose(fid);
  end
  %
  par.upen.B = B;
  return;
%
end

