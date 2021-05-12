%---------- File Parsing Script  ----------
%This script runs the NAV and OBS files parsers, and saves the parsed data
%as a .mat file. The need for a separate script is due to the very long time
%that it takes to parse the NAV and OBS files.


%% 0-Clear the workspace

clc
clear
format longg
tic() %DEBUG

%% 1-Set the file paths

%Input the name of the output file
outPath = "MGUE00ARG_R_20202990000_TEMP.mat";

%Input RINEX observations file
obs_file_path = "data/RINEX examples/MGUE00ARG_R_20202990300_01H_01S_GO.rnx";

%Input RINEX navigation file
nav_file_path = "data/RINEX examples/MGUE00ARG_R_20202982359_25H_GN_processed.nav";


%% 2-Open the input files

obs_fp=fopen(obs_file_path,'r');
if (obs_fp == -1)
    disp("Error reading Observation file")
    return;
end

nav_fp=fopen(nav_file_path,'r');
if (nav_fp == -1)
    disp("Error reading Navigation file")
    return;
end


%% 3-Parse the files

disp("Parsing RINEX observation data")
[gps_observations,sbas_observations] = parseRINEX_obs(obs_fp);

%disp("Parsing RINEX navigation data")
[alpha,beta,a,gps_navigation,sbas_navigation] = parseRINEX_nav(nav_fp);

%% 4-Sort the Navigation Messages by time of transmission
%It is assumed that the input file is sorted by time and not PRN

gps_navigation = sortRINEX_nav(gps_navigation);


%% 5-Save the variables to a .mat file
save(outPath,'gps_observations','alpha','beta','a','gps_navigation');


%% Finished sound
toc()
audio = load('train');
sound(audio.y,audio.Fs)