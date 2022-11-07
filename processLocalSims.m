%=============================================================
% Purpose: To create a large structure by processing stochastic simulations and save laxity and force/length data for ligaments and tendons
%
% Output:
% 	StochData structure with fields (time, lig: forces, msl: forces, knee: add laxity and moment)
%
% Other m-files required:
% 	none (nested files at bottom)
%
% Written by: Josh Roth
% Revised by: Matthew Blomquist
%
% Project: bam008-Matthew_kneeInstability
%
% -----------------------------------------
%
% Revision history:
% v1	2022-03-21
%
%=============================================================
function StoData = processLocalSims( Params )

resultsDir = Params.resultsFolder ;

% Temporary variables for this tutorial:
trialNames = { 'flexion' } ;
iTrial = 1 ;


%% Read results
%==============

laxResultFile = 'passive_flexion_states.sto' ;
frcResultsFile = 'passive_flexion_ForceReporter_forces.sto' ;

[ tempData.kine.data, tempData.kine.labels, tempData.kine.header ] = ...
    read_opensim_mot( fullfile( resultsDir, laxResultFile ) ) ; % read in values from lax files

[ tempData.frc.data, tempData.frc.labels, tempData.frc.header ] = ...
    read_opensim_mot( fullfile( resultsDir, frcResultsFile ) ) ; % read in values from force files

StoData.idx = parseOpenSimSto( tempData, Params ) ; % find indices for each trial


% Kinematic data
%================
if Params.kinematicData
    kineDofs = fieldnames( StoData.idx.kine ) ; % fieldnames of kinematic degrees of freedom

    for tempDof = 1 : length( kineDofs )
        tempDofName = kineDofs{ tempDof } ;
        if contains( tempDofName, 'vv' ) || contains( tempDofName, 'ie' )  || contains( tempDofName, 'fe' ) % rotation values
            scaleFactor = 180 / pi() ; % radians to degrees
            StoData.( trialNames{ iTrial } ).kine.( tempDofName )( : , 1 ) = tempData.kine.data( : , StoData.idx.kine.(tempDofName) ) * scaleFactor ;
        elseif contains( tempDofName , 'ap' )
            flexAngles = StoData.( trialNames{ iTrial } ).kine.fe( : , 1 ) ;
            scaleFactor = 1000 ; % m to mm
            StoData.( trialNames{ iTrial } ).kine.( tempDofName )( : , 1 ) = ...
                scaleFactor * ( cosd( flexAngles ) .* tempData.kine.data( : , StoData.idx.kine.ap ) - sind( flexAngles ) .* tempData.kine.data( : , StoData.idx.kine.ap ) ) ;
        elseif contains( tempDofName , 'pd' )
            flexAngles = StoData.( trialNames{ iTrial } ).kine.fe( : , 1 ) ;
            scaleFactor = 1000 ; % m to mm
            StoData.( trialNames{ iTrial } ).kine.( tempDofName )( : , 1 ) = ...
                scaleFactor * ( sind( flexAngles ) .* tempData.kine.data( : , StoData.idx.kine.ap ) + cosd( flexAngles ) .* tempData.kine.data( : , StoData.idx.kine.ap ) ) ;
        elseif contains( tempDofName , 'ml' )
            scaleFactor = 1000 ; % m to mm
            StoData.( trialNames{ iTrial } ).kine.( tempDofName )( : , 1 ) = tempData.kine.data( : , StoData.idx.kine.(tempDofName) ) * scaleFactor ;
        end


    end
end

% Muscle data
%=============

if ~isempty( Params.muscleNames )
    mslNames = fieldnames( StoData.idx.msl ) ;

    for iMsl = 1 : length( mslNames )
        tempMslName = mslNames{ iMsl } ;
        StoData.( trialNames{ iTrial } ).msl.( tempMslName ).frc( : , 1 ) = tempData.frc.data( : , StoData.idx.msl.( tempMslName ).frc ) ; % muscle force
    end
end

% Ligament data
%====================

if ~isempty( Params.ligamentNames )
    ligNames = fieldnames( StoData.idx.lig ) ;

    for iLig = 1 : length( ligNames )

        tempLigName = ligNames{ iLig } ;
        tempLigFrc = zeros( length( tempData.frc.data( : , 1 ) ) , 3 ) ;

        for iStrand = 1 : size( StoData.idx.lig.( ligNames{ iLig } ).frc_total, 1 )
            tempLigFrc( : , 1 ) = tempLigFrc( : , 1 ) + tempData.frc.data( : , StoData.idx.lig.( ligNames{ iLig } ).frc_spring( iStrand ) ) ; % spring force
            tempLigFrc( : , 2 ) = tempLigFrc( : , 2 ) + tempData.frc.data( : , StoData.idx.lig.( ligNames{ iLig } ).frc_damp( iStrand ) ) ; % damp force
            tempLigFrc( : , 3 ) = tempLigFrc( : , 3 ) + tempData.frc.data( : , StoData.idx.lig.( ligNames{ iLig } ).frc_total( iStrand ) ) ; % total force

            StoData.( trialNames{ iTrial } ).lig.( tempLigName ).( [ 'strand' , num2str( iStrand ) ] ).frcSpring( : , 1 ) = ...
                tempData.frc.data( : , StoData.idx.lig.( ligNames{ iLig } ).frc_spring( iStrand ) ) ; % spring force
            StoData.( trialNames{ iTrial } ).lig.( tempLigName ).( [ 'strand' , num2str( iStrand ) ] ).frcDamp( : , 1 ) = ...
                tempData.frc.data( : , StoData.idx.lig.( ligNames{ iLig } ).frc_damp( iStrand ) ) ; % damp force
            StoData.( trialNames{ iTrial } ).lig.( tempLigName ).( [ 'strand' , num2str( iStrand ) ] ).frcTotal( : , 1 ) = ...
                tempData.frc.data( : , StoData.idx.lig.( ligNames{ iLig } ).frc_total( iStrand ) ) ; % total force

            StoData.( trialNames{ iTrial } ).lig.( tempLigName ).( [ 'strand' , num2str( iStrand ) ] ).strain( : , 1 ) = ...
                tempData.frc.data( : , StoData.idx.lig.( ligNames{ iLig } ).strain( iStrand ) ) ;
            StoData.( trialNames{ iTrial } ).lig.( tempLigName ).( [ 'strand' , num2str( iStrand ) ] ).length( : , 1 ) = ...
                tempData.frc.data( : , StoData.idx.lig.( ligNames{ iLig } ).length( iStrand ) ) ;

        end

        StoData.( trialNames{ iTrial } ).lig.( tempLigName ).allStrands.frcSpring( : , 1 ) = tempLigFrc( : , 1 ) ; % add tension from combined strands
        StoData.( trialNames{ iTrial } ).lig.( tempLigName ).allStrands.frcDamp( : , 1 ) = tempLigFrc( : , 2 ) ; % add tension from combined strands
        StoData.( trialNames{ iTrial } ).lig.( tempLigName ).allStrands.frcTotal( : , 1 ) = tempLigFrc( : , 3 ) ; % add tension from combined strands

    end
end

% Contact force data
%===================

if Params.contactForceData
    contactCompartmentNames = Params.contactCompartmentNames ;
    contactLoadNames = Params.contactForces ;
    for iCntComp = 1 : length( contactCompartmentNames )
        for iCntLoad = 1 : length( contactLoadNames )
            StoData.( trialNames{ iTrial } ).cnt.( contactCompartmentNames{ iCntComp } ).( contactLoadNames{ iCntLoad } )( : , 1 ) = ...
                tempData.frc.data( : , StoData.idx.cnt.( contactCompartmentNames{ iCntComp } ).( contactLoadNames{ iCntLoad } ) ) ;
        end
    end
end

end


% end

%% Nested Functions
%==================

function Idx = parseOpenSimSto( StoDataStruct, Params )
%=======================================================================
% Purpose: Parse sto-file outputs from OpenSim-JAM simulations using
% simpleKnee_lab04.osim model (has 1 strand for MCL, LCL, ACL, and PCL)
%
% Inputs:
% 	StoDataStruct: data and labels from OpenSim-JAM (structure)
% 	Params: parameters for simulation (structure)
%
% Outputs:
% 	Idx: index values for data columns (structure)
%
% Example:
% 	Idx = parseOpenSimSto( StoDataStruct, Params )
% 	- Important parameter is Params.simType, which determines which data is
% 	pulled
%
% Other m-files required:
%	None
%
% Written By: Josh Roth
%
% Project: OpenSim-JAM
%
% --------------------------------------------
%
% Revision History:
% v1	2020-11-11	initial release (JDR)
%=======================================================================

%% Parse sto file
%===============

switch Params.model

    case 'lenhart2015'
        kineLabels = StoDataStruct.kine.labels ;
        frcLabels = StoDataStruct.frc.labels ;

        Idx.time = find( strcmp( kineLabels, 'time' )) ;

        %find indices of kinematics
        %--------------------------
        Idx.kine.fe = find( strcmp( kineLabels, '/jointset/knee_r/knee_flex_r/value' ) ) ;
        Idx.kine.vv = find( strcmp( kineLabels, '/jointset/knee_r/knee_add_r/value' ) ) ;
        Idx.kine.ie = find( strcmp( kineLabels, '/jointset/knee_r/knee_rot_r/value' ) ) ;
        Idx.kine.ap = find( strcmp( kineLabels, '/jointset/knee_r/knee_tx_r/value' ) ) ;
        Idx.kine.pd = find( strcmp( kineLabels, '/jointset/knee_r/knee_ty_r/value' ) ) ;
        Idx.kine.ml = find( strcmp( kineLabels, '/jointset/knee_r/knee_tz_r/value' ) ) ;

        %find indices of ligaments for force, length, and strain
        %--------------------------------------------------------
        %         ligNames = { 'dMcl', 'sMcl', 'pmc', 'lcl', 'itb', 'pfl', ...
        %             'aclam', 'aclpl', 'pclpm', 'pclal', ...
        %             'pCap', 'pt', 'lPfl', 'mPfl' } ;

        ligNames = Params.ligamentNames ;

        for iLig = 1 : length( ligNames )
            ligId = ligNames{ iLig } ;

            Idx.lig.( ligNames{ iLig } ).frc_spring = ...
                find( strncmpi( frcLabels, ligId, length( ligId ) ) & endsWith( frcLabels, 'force_spring' ) ) ;
            Idx.lig.( ligNames{ iLig } ).frc_damp = ...
                find( strncmpi( frcLabels, ligId, length( ligId ) ) & endsWith( frcLabels, 'force_damping' ) ) ;
            Idx.lig.( ligNames{ iLig } ).frc_total = ...
                find( strncmpi( frcLabels, ligId, length( ligId ) ) & endsWith( frcLabels, 'force_total' ) ) ;
            Idx.lig.( ligNames{ iLig } ).strain = ...
                find( strncmpi( frcLabels, ligId, length( ligId ) ) & endsWith( frcLabels, 'strain' ) ) ;
            Idx.lig.( ligNames{ iLig } ).length = ...
                find( strncmpi( frcLabels, ligId, length( ligId ) ) & endsWith( frcLabels, 'length' ) ) ;

            clear ligId
        end


        %find indices of muscles
        %-----------------------
        mslNames = {'addbrev_r';'addlong_r';'addmagProx_r';'addmagMid_r';'addmagDist_r';'addmagIsch_r';...
            'bflh_r'; 'bfsh_r';'edl_r';'ehl_r';'fdl_r';'fhl_r';'gaslat_r';'gasmed_r';...
            'gem_r';'glmax1_r';'glmax2_r';'glmax3_r';'glmed1_r';'glmed2_r';'glmed3_r';...
            'glmin1_r';'glmin2_r';'glmin3_r';'grac_r';'iliacus_r';'pect_r';'perbrev_r';...
            'perlong_r';'pertert_r';'piri_r';'psoas_r';'quadfem_r';'recfem_r';'sart_r';'semimem_r';...
            'semiten_r';'soleus_r';'tfl_r';'tibant_r';'tibpost_r';'vasint_r';'vaslat_r';'vasmed_r' };

        numMsl = length( mslNames ) ;

        for iMsl = 1 : numMsl
            Idx.msl.( mslNames{ iMsl } ).frc = find( strcmp( frcLabels, mslNames{ iMsl } ) ) ;
        end

        %find indices of contact forces
        %------------------------------
        contactCompartmentNames = Params.contactCompartmentNames ;
        contactLoadNames = Params.contactForces ;

        for iCntComp = 1 : length( contactCompartmentNames )
            for iCntLoad = 1 : length( contactLoadNames )
                Idx.cnt.( contactCompartmentNames{ iCntComp } ).( contactLoadNames{ iCntLoad } ) = ...
                    find( strncmpi( frcLabels, contactCompartmentNames{ iCntComp }, length( contactCompartmentNames{ iCntComp } ) ) & ...
                    contains( frcLabels, contactLoadNames{ iCntLoad } ) ) ;
            end
        end


end

end


function [data, labels, header] = read_opensim_mot(file)
%%=========================================================================
%READ_OPENSIM_MOT
%--------------------------------------------------------------------------
%Author(s): Colin Smith
%Date: 5/14/2018

%The kneemos_matlab toolkit is a collection of code for developing and
%analyzing musculoskeletal simulations in SIMM and OpenSIM. The developers
%are based at the University of Wisconsin-Madison and ETH Zurich. Please
%see the README.md file for more details. It is your responsibility to
%ensure this code works correctly for your use cases.
%
%Licensed under the Apache License, Version 2.0 (the "License"); you may
%not use this file except in compliance with the License. You may obtain a
%copy of the License at http://www.apache.org/licenses/LICENSE-2.0.
%Unless required by applicable law or agreed to in writing, software
%distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
%WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
%License for the specific language governing permissions and limitations
%under the License.
%
%--------------------------------------------------------------------------
%file : str
%
%
%data : [nFrames x nLabels] matrix
%
%
%labels : {1 x nLabels} cell of strings
%
%
%header : struct
%
%
%==========================================================================

if ~exist('file','var')
    [infile, inpath]=uigetfile('*.mot','Select input file');
    file=[inpath infile];
end

fid=fopen(file,'r');

if fid <0
    mot=[];labels=[];
    disp('File Not Found:\n');
    disp([file '\n']);
    return
end

disp(['Loading file...' file] );


%read the file name line
header.filename=strtrim(fgetl(fid));


% Read Header
line = fgetl(fid);
while ~strncmpi(line,'endheader',length('endheader'))


    if (line == -1)
        disp('ERROR: Reached EOF before "endheader"')
        return
    end
    line_space = strrep(line,'=',' ');
    split_line = strsplit(line_space);

    if (length(split_line)==2)
        var = split_line{1};
        value = split_line{2};

        if strcmpi(var,'version')
            header.version = str2double(value);
        elseif strcmpi(var,'nRows') || strcmpi(var,'datarows')
            nr = str2double(value);
            header.nRows = nr;
        elseif strcmpi(var,'nColumns') || strcmpi(var,'datacolumns')
            nc = str2double(value);
            header.nColumns = nc;
        elseif strcmpi(var,'indegrees')
            header.indegrees = strtrim(value);
        end
    end

    line = fgetl(fid);
end


%Load Column Headers
line=fgetl(fid);

labels=cell(nc,1);

j=1;
jl=length(line);
for i=1:nc
    name=sscanf(line(j:jl),'%s',1);
    ii = findstr(line(j:jl), name);
    j=j+ii(1)+length(name);
    labels(i,1)=cellstr(name);
end

% Now load the data
data=zeros(nr,nc);
i=0;
while ((feof(fid)==0)&&(i<nr))
    i=i+1;
    line=fgetl(fid);
    try
        data(i,:)=sscanf(line,'%f');
    catch
        line = strrep( line, 'nan(ind)', '0' ) ;
        data(i,:)=sscanf(line,'%f');
        disp( ' ' )
        disp( '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' )
        disp( '!Warning! Data contained nan-values that were set to 0!' )
        disp( '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' )
        disp( ' ' )
    end
end
if (i<nr)
    disp(['Number of rows (',num2str(i),') is less than that specified in header (',num2str(nr),')']);
    data=data(1:i,:);
end
fclose(fid);



if (nargout>1)
    % return all data in a single array
    mot=data;
elseif (nargout==1)
    % return all information in a single structure
    mot.data=data;
    mot.hdr=labels;
end


    function [t,q]=load_exp(file);
        global NQ NM;
        rawdata=load_motionfile(file);
        [t,data]=extractcolumns(rawdata,1);
        [q,data]=extractcolumns(data,NQ);
    end


    function [x,outdata]=extractcolumns(data,nc);
        x=data(:,1:nc);
        [m,n]=size(data);
        outdata=data(:,(nc+1):n);
    end

end
