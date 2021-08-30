%% 
% Code to add GGA messages to EK60 files, based on RMC messages that are
% already in the file.

% Trip 78 - DONE
dataDir = 'E:\Aqualyd\SIO_ORH\Data\WW ES80 2018-2021\Trip 78\Hull\ES60';
outDir = 'E:\Aqualyd\SIO_ORH\Data\WW ES80 2018-2021\Trip 78\Hull\ES60-GGA';

% Trip 67
dataDir = 'E:\Aqualyd\SIO_ORH\Data\2005-17 selected\2016 Acoustics\Trip 67';
outDir = 'E:\Aqualyd\SIO_ORH\Data\2005-17 selected\2016 Acoustics\Trip 67-GLL';

dataDir = 'E:\Aqualyd\SIO_ORH\Data\2005-17 selected\2016 Acoustics\Trip 67\survey';
outDir = 'E:\Aqualyd\SIO_ORH\Data\2005-17 selected\2016 Acoustics\Trip 67-GLL\survey';

dataDir = 'E:\Aqualyd\SIO_ORH\Data\2005-17 selected\2016 Acoustics\Trip 67\open ocean survey';
outDir = 'E:\Aqualyd\SIO_ORH\Data\2005-17 selected\2016 Acoustics\Trip 67-GLL\open ocean survey';

% Trip 77 - DONE
dataDir = 'E:\Aqualyd\SIO_ORH\Data\WW ES80 2018-2021\Trip 77\All Trip';
outDir = 'E:\Aqualyd\SIO_ORH\Data\WW ES80 2018-2021\Trip 77\All Trip-GLL';
% D20180624-T234738.raw failed
% D20180704-T010518.raw failed

dataDir = 'E:\Aqualyd\SIO_ORH\Data\2005-17 selected\2016 Acoustics\Trip 67\calibration willwtch 2016 10 5';
outDir = 'E:\Aqualyd\SIO_ORH\Data\2005-17 selected\2016 Acoustics\Trip 67-GLL\calibration willwtch 2016 10 5';


d = dir(fullfile(dataDir, '*.raw'));

numFiles = length(d);

for i = 1:numFiles
    disp(['Doing ' d(i).name ' (' num2str(i) ' of ' num2str(numFiles) ')'])
    dfile = fullfile(d(i).folder, d(i).name);
    ofile = fullfile(outDir, d(i).name);
    
%   % read in the raw file
%   [header, data] = readEKRaw(dfile, 'VesselSpeed', 1, 'VTGSource', 'GPRMC', ...
%       'GPSSource', 'GPRMC', 'RawNMEA', true, 'RequireGPSChecksum', false);
   
%   % write out the raw file
%    writeEKRaw(ofile, header, data, 'GPS', true, 'VesselSpeed', true);

    headerlength = cHeader.length(); % bytes

    fid = fopen(dfile,'r');
    nfid = fopen(ofile, 'w');

    if (fid == -1)
        warning(['Could not open file ' dfile]);
    else
        while(1)
            dglength = fread(fid, 1, 'int32'); % the datagram header
            if feof(fid)
                break
            end
            header = cHeader;
            
            % read the datagram
            header = header.read(fid);
            dgData = fread(fid, dglength-headerlength);
            fread(fid, 1, 'int32'); % the trailing datagram marker

            % if RMC datagram, make up a GLL datagram
            if strcmp(header.type, 'NME0')
                %
                nmeadata = char(dgData');
                if strncmp(nmeadata, '$GPRMC', 6)
                   %format = '%2d %2d %f %c %2d %f %c %3d %f %c %f %f %6d %f';
                   %[out{1:14}] = strread(nmeadata, format, 1, 'delimiter', ',');
                   out = split(nmeadata, ',');
                   % make a GLL from the RMC
                   gll = ['GPGLL,' out{4} ',' out{5} ',' out{6} ',' out{7} ',' out{2}];
                   % write the gll out as a new datagram, using the same
                   % header as for the RMC message
                   dgLength = length(gll) + header.length();
                   fwrite(nfid, dgLength, 'int32');
                   header.write(nfid)
                   fwrite(nfid, gll);
                   fwrite(nfid, dgLength, 'int32');
                end
            end
            
            
            % write out the datagram we just read in
            dgLength = length(dgData) + header.length();
            fwrite(nfid, dgLength, 'int32');
            header.write(nfid)
            fwrite(nfid, dgData);
            fwrite(nfid, dgLength, 'int32');

        end
        fclose(fid);
        fclose(nfid);
    end
end
