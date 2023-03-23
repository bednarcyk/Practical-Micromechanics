% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Copyright 2020 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration. No copyright is claimed in the 
% United States under Title 17, U.S. Code. All Other Rights Reserved. BY DOWNLOADING 
% OR USING THIS SOFTWARE, YOU ACKNOWLEDGE THAT YOU HAVE READ THE NASA OPEN SOURCE 
% AGREEMENT V1.3, THAT YOU UNDERSTAND IT, AND THAT YOU AGREE TO BE BOUND BY ITS 
% TERMS. IF YOU DO NOT AGREE TO THE TERMS AND CONDITIONS OF THIS AGREEMENT, DO NOT 
% USE OR DOWNLOAD THE SOFTWARE. THIS SOFTWARE IS PROVIDED AS IS WITHOUT ANY WARRANTY 
% OF ANY KIND. RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST, AND INDEMNIFIES 
% AND HOLDS HARMLESS, THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND 
% SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT. This code was prepared by Drs. 
% B.A. Bednarcyk and S.M. Arnold to complement the book “Practical Micromechanics of 
% Composite Materials” during the course of their government work.
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% Purpose: Driver script for laminate analysis problems using CLT. The problem input 
%          is defined in the function LamProblemDef.m and GetPlyProps.m
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Clear memory and close files
clear;
close all;
fclose('all');
clc;

% -- Add needed function locations to the path
addpath('Functions/CLT');
addpath('Functions/Utilities');
addpath('Functions/WriteResults');
addpath('Functions/Micromechanics');
addpath('Functions/Margins');

%-----------------------------------------------------------------
% 1) Define Laminate Problems
%-----------------------------------------------------------------
[NProblems, OutInfo, Geometry, Loads] = LamProblemDef();

% -- Preallocate LamResults
LamResults = cell(1, NProblems);

%-----------------------------------------------------------------
% 2) Get ply properties
%-----------------------------------------------------------------
plyprops = GetPlyProps(Geometry);

% -- Loop through problems
for NP = 1: NProblems

    % -- Check for missing problem name
    if (ismissing(OutInfo.Name(NP)))
        OutInfo.Name(NP) = string(['Problem ', char(num2str(NP))]);
    end
    
    % -- Check for missing ply angles, thicknesses, materials
    if ~isfield(Geometry{NP},'Orient') || ~isfield(Geometry{NP},'tply') || ...
       ~isfield(Geometry{NP},'plymat')
        error(['Orient, tply, or plymat missing, Problem #', num2str(NP)]);
    end
    
    % -- Echo problem info to command window
    disp(['Problem #',num2str(NP),' - ', char(OutInfo.Name(NP))]);
    disp(['   Per ply angle orientations   [',num2str(Geometry{NP}.Orient), ']']);
    disp(['   Per ply thicknesses          [',num2str(Geometry{NP}.tply), ']']);
    disp(['   Per ply material assignments [',num2str(Geometry{NP}.plymat), ']']);    
    
    % -- Check that the problem's ply materials, orientation, and thickness
    %    have consistent lengths
    if length(Geometry{NP}.tply) ~= length(Geometry{NP}.Orient) || ...
       length(Geometry{NP}.tply) ~= length(Geometry{NP}.plymat) || ...
       length(Geometry{NP}.Orient) ~= length(Geometry{NP}.plymat)
            error('Lengths of tply, Orient, and plymat are not consistent');
    end

    % -- Check that loads are specified for this problem
    if  ~isfield(Loads{NP}, 'Type') || ~isfield(Loads{NP}, 'Value')   
        error(strcat('Problem #', num2str(NP), ' Loads not properly defined'));
    end   
    if ~isfield(Loads{NP}, 'DT') % -- Default to zero DT if not specified
        Loads{NP}.DT = 0;
    end

    % -- Check that the problem's ply materials have been defined
    for k = 1: length(Geometry{NP}.tply)
        if ~isfield(plyprops{Geometry{NP}.plymat(k)}, 'name')
            error(strcat('ply material #', num2str(Geometry{NP}.plymat(k)), ...
                         ' undefined ... check GetPlyProps'));
        end
    end

    % -- Check for slash character in OutInfo.Name
    k = strfind(OutInfo.Name(NP), '/');
    j = strfind(OutInfo.Name(NP), '\');
    if ~isempty(k) || ~isempty(j)
        error('Problem name contains a slash or backslash ... remove');
    end
        
    %-----------------------------------------------------------------
    % 3) Analyze laminate with CLT per problem
    %-----------------------------------------------------------------
    [Geometry{NP}, LamResults{NP}] = CLT(plyprops, Geometry{NP}, Loads{NP});

    %-----------------------------------------------------------------
    % Orthotropic SS Plate Deflection Analysis
    %-----------------------------------------------------------------

    % -- Plate Dimensions
    a_in = 12;     % -- inches
    b_in = 12;
    a = a_in*25.4; % -- mm
    b = b_in*25.4;

    % -- Load
    W_lb = 70;          % -- Applied weight over panel (lb)
    W_N = W_lb*4.44822; % -- Applied force over panel (N)
    p0 = W_N/(a*b);     % -- Applied pressure on panel (N/mm^2)
    
    % -- Number of points in x and y directions (should be odd)
    NumX = 51;
    NumY = 51;
    
    % -- Number of terms in the series (should be odd)
    Nterms = 21;
    
    % -- Effective Torsional Rigidity (H)
    H = LamResults{NP}.D(1,2) + 2*LamResults{NP}.D(3,3);
    
    % -- Bending Stiffnesses
    D11 = LamResults{NP}.D(1,1);
    D22 = LamResults{NP}.D(2,2);
    
    % -- Increments in the x and y coordinates
    x_inc = a/(NumX - 1);
    y_inc = b/(NumY - 1);
    
    % -- Calculate plate deflection Navier solution base on double series
    w = zeros(NumX,NumY);
    err = zeros(NumX,NumY);
    x = 0;
    for i = 1: NumX
        y = 0;
        for j = 1: NumY
            
            xx(i,j) = x;
            yy(i,j) = y;
            
            for m = 1:2:Nterms
                for n = 1:2:Nterms
                    
                    NUM = sin(m*pi*x/a)*sin(n*pi*y/b);
                    
                    DEN = m*n*(D11*m^4/a^4 + 2*H*m^2*n^2/(a^2*b^2) + D22*n^4/b^4); % -- orthotropic
                    %DEN = D11*m*n*((m/a)^2 + (n/b)^2)^2; % -- isotropic
                    
                    Add = (16*p0/pi^6)*NUM/DEN;
                    w(i,j) = w(i,j) + Add;
                    
                    last_add(i,j) = Add;
                end
            end
            
            if w(i,j) ~= 0
                err(i,j) = abs( last_add(i,j)/w(i,j) );
            end
            
            y = y + y_inc;
        end
        x = x + x_inc;
    end
    
    % -- surf plot
    figure;
    % -- w is positive down, so plotting -w
    surf(xx(1:NumX,1),yy(1,1:NumY),-w);   
    colormap(flipud(jet));
    colorbar;
    set(gcf,'color','w');
    hold on;
    
    % -- max deflection
    [M,I] = max(w(:));
    [ind1, ind2] = ind2sub(size(w),I);
    
   
    % -- Plate mass
    mass = 0;
    for i = 1:length(Geometry{NP}.plymat)
        if isfield(plyprops{Geometry{NP}.plymat(i)}, 'rho')
             mass = mass + a * b * Geometry{NP}.tply(i) * plyprops{Geometry{NP}.plymat(i)}.rho;
        else
            mass = nan;
        end
    end

    disp(' ');
    disp(['   Max Error:       ',char(num2str(max(err(:))))])
    disp(['   Max Deflection:  ',char(num2str(M))])
    disp(['   x,y location:   (',char(num2str(xx(ind1,ind2))), ',',char(num2str(yy(ind1,ind2))),')'])
    if ~isnan(mass)
        disp(['   Plate Mass:      ',char(num2str(mass))])
    end
    disp(' ');
   
%     disp(D11);
%     disp(D22);
%     disp(LamResults{NP}.D(1,2));
%     disp(LamResults{NP}.D(3,3));
    
    % -- Display notification of file completion to command window
    disp(['  *** Problem ',char(num2str(NP)),' Completed ***'])
    disp(' ');
    %close all;
    
end