%% 
% Code to add GGA messages to EK60 files, based on RMC messages that are
% already in the file.

% Trip 78 - DONE
dataDir = 'E:\Aqualyd\SIO_ORH\Data\WW ES80 2018-2021\Trip 78\Hull\ES60';
outDir = 'E:\Aqualyd\SIO_ORH\Data\WW ES80 2018-2021\Trip 78\Hull\ES60-GGA';

% Trip 67 - DONE
dataDir = 'E:\Aqualyd\SIO_ORH\Data\2005-17 selected\2016 Acoustics\Trip 67';
outDir = 'E:\Aqualyd\SIO_ORH\Data\2005-17 selected\2016 Acoustics\Trip 67-GLL';

dataDir = 'E:\Aqualyd\SIO_ORH\Data\2005-17 selected\2016 Acoustics\Trip 67\survey';
outDir = 'E:\Aqualyd\SIO_ORH\Data\2005-17 selected\2016 Acoustics\Trip 67-GLL\survey';

dataDir = 'E:\Aqualyd\SIO_ORH\Data\2005-17 selected\2016 Acoustics\Trip 67\open ocean survey';
outDir = 'E:\Aqualyd\SIO_ORH\Data\2005-17 selected\2016 Acoustics\Trip 67-GLL\open ocean survey';

dataDir = 'E:\Aqualyd\SIO_ORH\Data\2005-17 selected\2016 Acoustics\Trip 67\calibration willwtch 2016 10 5';
outDir = 'E:\Aqualyd\SIO_ORH\Data\2005-17 selected\2016 Acoustics\Trip 67-GLL\calibration willwtch 2016 10 5';

% Trip 77 - DONE
dataDir = 'E:\Aqualyd\SIO_ORH\Data\WW ES80 2018-2021\Trip 77\All Trip';
outDir = 'E:\Aqualyd\SIO_ORH\Data\WW ES80 2018-2021\Trip 77\All Trip-GLL';
% D20180624-T234738.raw failed
% D20180704-T010518.raw failed



mainDir = 'E:\Aqualyd\SIO_ORH\Data\2005-17 selected';
mainDir = 'E:\Aqualyd\SIO_ORH\Data\Survey_data2';

dirs = dir(mainDir);

for j = 1:length(dirs)
    if dirs(j).isdir == 1 && dirs(j).name ~= "." && dirs(j).name ~= ".."
        dataDir = fullfile(dirs(j).folder, dirs(j).name);
        d = dir(fullfile(dataDir, '**/*.raw'));
        
        numFiles = length(d);
        disp(['Processing top level directory: ' dirs(j).name])
        % output directory will be under dataDir, but with '-NMEA' added
        for i = 1:numFiles
            disp(['Doing ' fullfile(d(i).folder, d(i).name) ' (' num2str(i) ' of ' num2str(numFiles) ')'])

            outDir = replace(d(i).folder, dataDir, [dataDir '-NMEA']);
            if ~isfolder(outDir)
                mkdir(outDir)
            end

            dfile = fullfile(d(i).folder, d(i).name);
            ofile = fullfile(outDir, d(i).name);

            headerlength = cHeader.length(); % bytes

            fid = fopen(dfile,'r');
            nfid = fopen(ofile, 'w');

            if (fid == -1)
                warning(['Could not open file ' dfile]);
            else
                try
                    while(1) % read in each datagram
                        dglength = fread(fid, 1, 'int32'); % the datagram header
                        if feof(fid)
                            break
                        end
                        header = cHeader;

                        % read the datagram
                        header = header.read(fid);
                        dgData = fread(fid, dglength-headerlength);
                        fread(fid, 1, 'int32'); % the trailing datagram marker

                        % if RMC datagram, make up a GLL and VTG datagram
                        if strcmp(header.type, 'NME0')
                            %
                            nmeadata = char(dgData');
                            if strncmp(nmeadata, '$GPRMC', 6)
                               out = split(nmeadata, ',');
                               if length(out) >= 8
                                   % make a GLL from the RMC
                                   gll = ['GPGLL,' out{4} ',' out{5} ',' out{6} ',' out{7} ',' out{2}];
                                   % write the gll out as a new datagram, using the same
                                   % header as for the RMC message
                                   dgLength = length(gll) + header.length();
                                   fwrite(nfid, dgLength, 'int32');
                                   header.write(nfid)
                                   fwrite(nfid, gll);
                                   fwrite(nfid, dgLength, 'int32');

                                   % make a VTG from the RMC
                                   vtg = ['GPVTG,' out{9} ',T,,M,' out{8} ',N,,K,'];
                                   % write the vtg out as a new datagram, using the same
                                   % header as for the RMC message
                                   dgLength = length(vtg) + header.length();
                                   fwrite(nfid, dgLength, 'int32');
                                   header.write(nfid)
                                   fwrite(nfid, vtg);
                                   fwrite(nfid, dgLength, 'int32');
                               end
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
                catch ME
                    % probably because the file ended unexpectedly, so just go to
                    % the next file.
                    disp(['Skipping to next file due to error: "' ME.message '" at line ' num2str(ME.stack.line)])
                    fclose(fid);
                    fclose(nfid);
                end
            end
        end
    end
end