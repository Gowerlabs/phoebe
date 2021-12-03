function measurement_matrix = parse_fnirs_channels_nirstar(inlet)

% get the full stream info (including custom meta-data) and dissect it
inf = inlet.info();
tags = inf.desc();
tags.first_child().name();
tags.first_child().next_sibling().name();
tags.first_child().next_sibling().next_sibling().name();

%Get information about the stream
stream_info = inlet.info();

% Get individual stream attributes
stream_name= stream_info.name();
stream_mac= stream_info.type();
stream_n_channels = stream_info.channel_count();

%% Scan LSL metadata to parse information about fNIRS channels only 

ch = stream_info.desc().child('channels').child('channel');

% Fill measurement list channel-by-channel
% Columns: [s d wl lsl_vector_index) 
nirs_channel = 1;
for k = 1:(stream_n_channels)  
    if strcmp(ch.child_value('type'),'nirs_raw')||strcmp(ch.child_value('type'),'nirs_hb')    %If the channel is NIRS, let's parse the metadata
        label = ch.child_value('label');
        wl = ch.child_value('wavelength');
        rg = regexp(label,'(?<s>\d+)-(?<d>\d+)','names'); %extract the source number s and detector number d
        meas_list(nirs_channel,1) = str2double(rg.s);
        meas_list(nirs_channel,2) = str2double(rg.d);
        meas_list(nirs_channel,3) = str2double(wl);
        meas_list(nirs_channel,4) = k;   % This is the index(position) of channel in the LSL vector being streamed
        nirs_channel = nirs_channel +1;   % Increase nirs channel counter 
    end
    ch = ch.next_sibling(); % Move to next LSL channel
end

%% Convert to measurement_matrix: [s d wl1_lsl_index wl2_lsl_index]

% Extract list of all wavelengths streamed by this device
wl_list = unique(meas_list(:,3));
%If wl_list is only two wavelengths, we are good. If more than two, we must select the two preferred wavelengths

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

