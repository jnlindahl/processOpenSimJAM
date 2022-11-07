%% Forward Simulation Example Wrapper
%
% Purpose: To customize inputs before using functions to generate files

%% Specify creating files or post-processing
%-------------------------------------------

% Specify whether you want to create files to run simulations, or post
% process results that have already been run
%   'createFiles' to create files, 'postProcess' to post process
purpose = 'postProcess' ;

switch purpose
    
    case 'postProcess'
        % Specify parameters (change these to desired parameters and file locations)
        %====================

        % Directory paths:
        %-----------------
        
        % Specify folder where results files are contained
        Params.resultsFolder = 'C:\Users\mbb201\Documents\MATLAB\Research\openSimTutorial\opensim-jam\opensim-jam-release\examples\passive_flexion\results\forsim' ;

        Params.model = 'lenhart2015' ;
        
        % Data to pull:
        %--------------
        
        % 1 to extract kinematic/kinetic data, 0 to not extract data
        Params.kinematicData = 1 ; % all 6 DOF
        Params.contactForceData = 1 ; % contact forces
        
        % Muscles to extract
        %   Leave as {} if you don't want to extract muscle data
        Params.muscleNames = { 'addbrev_r';'addlong_r';'addmagProx_r';'addmagMid_r';'addmagDist_r';'addmagIsch_r';...
            'bflh_r'; 'bfsh_r';'edl_r';'ehl_r';'fdl_r';'fhl_r';'gaslat_r';'gasmed_r';...
            'gem_r';'glmax1_r';'glmax2_r';'glmax3_r';'glmed1_r';'glmed2_r';'glmed3_r';...
            'glmin1_r';'glmin2_r';'glmin3_r';'grac_r';'iliacus_r';'pect_r';'perbrev_r';...
            'perlong_r';'pertert_r';'piri_r';'psoas_r';'quadfem_r';'recfem_r';'sart_r';'semimem_r';...
            'semiten_r';'soleus_r';'tfl_r';'tibant_r';'tibpost_r';'vasint_r';'vaslat_r';'vasmed_r' } ;
        
        % Ligaments to extract
        %   Leave as {} if you don't want to extract ligament data
        Params.ligamentNames = { 'MCLs', 'MCLd' , 'MCLp', 'ACLpl' , 'ACLam' , 'LCL', 'ITB', 'PFL', 'pCAP', ...
            'PCLpm', 'PCLal', 'PT', 'lPFL', 'mPFL' } ;

        % Compartments to extract contact forces from
        Params.contactCompartmentNames = { 'tf_contact', 'pf_contact' } ;
        
        % Contact force values to extract
        Params.contactForces = { 'mean_pressure' , 'max_pressure' , ...
            'contact_force_x' , 'contact_force_y' , 'contact_force_z' , ...
            'contact_moment_x' , 'contact_moment_y' , 'contact_moment_z' } ;
        
        %% Run post processing (no need to change anything in this portion)
        %======================
        
        % Run processLocalSims.m
        simData = processLocalSims( Params ) ;
        
        
end


%% ADD CODE HERE TO EXTRACT DATA OF INTEREST

% Switch this between 1 and 0. If you want to plot the results for example
% 1 (below), then set this to 1. If you don't want to plot the results for
% example 1, then set this to 0.
plotExample1 = 0 ;

if plotExample1
    % Example: If you wanted to plot LCL tension vs flexion angle
    figure() ; hold on ; grid on ;

    flexionAngles = simData.flexion.kine.fe ; % 251 x 1 double of flexion angles
    LCLtension = simData.flexion.lig.LCL.allStrands.frcTotal ; % 251 x 1 double of LCL tensions

    plot( flexionAngles , LCLtension )

    xlabel( 'Flexion Angle (^o)' ) ; ylabel( 'LCL Tension (N)' ) ;
end

% Switch this between 1 and 0. If you want to plot the results for example
% 1 (below), then set this to 1. If you don't want to plot the results for
% example 1, then set this to 0.
plotExample2 = 0 ;

if plotExample2
    % Example: If you wanted to plot tension of each MCLs strand vs flexion angle
    figure() ; hold on ; grid on ;

    flexionAngles = simData.flexion.kine.fe ; % 251 x 1 double of flexion angles
   
    % There are 6 strands for the MCLs. The best way to plot would be to
    % loop through them all using a for loop, but if you aren't as
    % comfortable with MATLAB, then this way should make more sense:
    MCLs1_tension = simData.flexion.lig.MCLs.strand1.frcTotal ; % 251 x 1 double
    MCLs2_tension = simData.flexion.lig.MCLs.strand2.frcTotal ; % 251 x 1 double
    MCLs3_tension = simData.flexion.lig.MCLs.strand3.frcTotal ; % 251 x 1 double
    MCLs4_tension = simData.flexion.lig.MCLs.strand4.frcTotal ; % 251 x 1 double
    MCLs5_tension = simData.flexion.lig.MCLs.strand5.frcTotal ; % 251 x 1 double
    MCLs6_tension = simData.flexion.lig.MCLs.strand6.frcTotal ; % 251 x 1 double

    % Plot each strand's tension vs flexion angle. They will all stay on
    % the same plot since we added "hold on" on line 92.
    plot( flexionAngles , MCLs1_tension )
    plot( flexionAngles , MCLs2_tension )
    plot( flexionAngles , MCLs3_tension )
    plot( flexionAngles , MCLs4_tension )
    plot( flexionAngles , MCLs5_tension )
    plot( flexionAngles , MCLs6_tension )

    xlabel( 'Flexion Angle (^o)' ) ; ylabel( 'MCLs Tension (N)' ) ;
    legend( { 'Strand 1' , 'Strand 2' , 'Strand 3' , 'Strand 4' , 'Strand 5' , 'Strand 6' } )
end

% Make various plots here to get more comfortable with the simData
% structure and plotting in MATLAB! Some things you could try (PLEASE reach
% out to me if you are unable to get any of these things, I'm happy to
% help!)

% 1) Tibiofemoral superior-inferior (S-I) contact force vs flexion angle
% (S-I is the "y" direction in OpenSim)
% 
% 2) Anterior-Posterior, Compression-Distraction, and Medial-Lateral
% kinematics on the same plot vs flexion angle
%
% 3) LCL length vs tension for each strand (there are 4 strands for the
% LCL)
%
% 4) Come up with something on your own and let me know and we can double
% check that we have the same results!