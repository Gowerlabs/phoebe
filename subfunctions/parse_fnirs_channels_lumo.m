function [measurement_matrix,type] = parse_fnirs_channels_lumo(inlet)
% PARSE_FNIRS_CHANNELS_LUMO Parse fNIRS channels from a LUMO LSL stream
% 
% [measurement_matrix,type] = PARSE_FNIRS_CHANNELS_LUMO(inlet)
%
% PARSE_FNIRS_CHANNELS_LUMO is used to parse fNIRS channel enumeration from
% a LUMO LSL stream. 
% 
% The code checks if the LUMO LSL stream contains valid channel enumeration 
% structure, and will return an error if invalid data is found.
% 
% Parameters
%   inlet: an open LUMO LSL inlet
%
% Returns
%   meausurement_matrix:    [#channels x 4] contains the S-D pairings (cols 1-2)
%                           and the LSL vector indices (cols 3-4) for either wavelength 1-2 
%
%   type: 0 for raw data, 1 for hemoglobin
% 
% ND, Gowerlabs, October 2024

% LUMO transmits raw data
type = 0; 

% Get information about the stream
stream_info = inlet.info();

% Get individual stream attributes
stream_mac= stream_info.type();  %Type NIRS
stream_n_channels = stream_info.channel_count(); %Get total channel count

% Validate stream attributes
if ~strcmp(stream_mac, 'NIRS') 
    error('Chosen LSL stream is not a fNIRS LUMO stream.');
end

%% Scan LSL metadata to parse information about fNIRS channels only 
ch = stream_info.desc().child('channels').child('channel');

% Fill measurement list channel-by-channel
% Columns: [s d wl lsl_vector_index) 
nirs_channel = 1;
meas_list = zeros(stream_n_channels,4); %Prepare for efficiency

for k = 1:(stream_n_channels)
    if strcmp(ch.child_value('type'),'Intensity')   

        % extract source information
        source = ch.child_value('source');
        rg1= regexp(source,'N(?<node>\d+)/(?<idx>\w)','names');
        src_node = str2double(rg1.node);
    
        if (src_node < 1 || src_node > 54)
            error('Chosen LSL stream has inconsistent source enumeration.');
        end

        switch(rg1.idx)
            case 'A'
                source_idx = 1;
            case 'B'
                source_idx = 2;
            case 'C'
                source_idx = 3;
            otherwise
                error('Chosen LSL stream has inconsistent source enumeration.');
        end
        
        % extract detector information
        detector = ch.child_value('detector');
        rg2= regexp(detector,'N(?<node>\d+)/(?<idx>\d)','names');
        det_node = str2double(rg2.node);
        det_idx = str2double(rg2.idx);

        if (det_node < 1 || det_node > 54) || (det_idx < 1 || det_idx > 4)
            error('Chosen LSL stream has inconsistent detector enumeration.');
        end

        % extract wavelength information 
        wavelen = ch.child_value('wavelen');
        if strcmp(wavelen,'735.000000')
           wl = 0;
        elseif strcmp(wavelen,'850.000000')
           wl = 1;
        else
           error('Chosen LSL stream has inconsistent wavelength information.');
        end

        meas_list(nirs_channel,1) = (src_node - 1)*3 + source_idx;
        meas_list(nirs_channel,2) = (det_node - 1)*4 + det_idx;
        meas_list(nirs_channel,3) = wl;
        meas_list(nirs_channel,4) = k;   % This is the index(position) of channel in the LSL vector being streamed
        nirs_channel = nirs_channel +1;   % Increase nirs channel counter 
    end
    ch = ch.next_sibling(); % Move to next LSL channel
end

%% Convert to measurement_matrix: [s d wl1_lsl_index wl2_lsl_index]

% Extract list of all wavelengths streamed by this device
wl_list = unique(meas_list(:,3));
%If wl_list is only two wavelengths, we are good. 
if size(wl_list, 1) ~= 2
    error('Chosen LSL stream has inconsistent number of wavelengths.');
end

measurement_matrix = zeros(size(meas_list,1)/length(wl_list),4); %Prepare for efficiency
channel_list = unique(meas_list(:,1:2),'rows');   % Pull unique list of s-d pairings

% For all occurrences of each s-d pairings, copy LSL index into column 3 and 4 of new matrix 
for i = 1: size(channel_list,1)
   [r1,~] = find(ismember(meas_list(:,1:3),[channel_list(i,:) wl_list(1)],'rows'));   % Find row where s-d pairing at wl1 is
   [r2,~] = find(ismember(meas_list(:,1:3),[channel_list(i,:) wl_list(2)],'rows'));   % Find row where s-d pairing at wl2 is
   measurement_matrix(i,1:2) = channel_list(i,:);   % Place unique s-d pairings
   measurement_matrix(i,3) = meas_list(r1,4);       % Place lsl_index of wl1
   measurement_matrix(i,4) = meas_list(r2,4);       % Place lsl_index of wl2
end