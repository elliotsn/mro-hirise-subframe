% 
% Function to read a portion of a HiRISE map-projected JP2. Assumes the
% label file is stored in the same directory as the HiRISE map-projected JP2.
%
% lims1 and lims2 are the min and max bounds to read in the JP2 of
% whichever coordinate is intended to be used.
%
% The string mode is used to specify which system is used:
%   System                              modestr
%   lon, lat                            lonlat
%   equidistant cylindrical x,y         xy
%   sample, line                        pix
%
% xvec,yvec      - equidistant cylindrical coords
% latvec,lonvec  - lat-lon vectors for image
%
% Elliot Sefton-Nash  
%
%  06/10/2017 Original
%  28/02/2018 Update to include modes of operation, provide either pixel
%             numbers, lat-lon or xy coords as boundaries of subframe.
function [lonvec, latvec, xvec, yvec, im, flbl, lbl, linelims, samplelims]  = ...
    readHiriseJp2Subframe(fjp2, lims1, lims2, modestr)

% Read the label and get tile info
flbl = [fjp2(1:end-3),'LBL'];
% Test jp2 and lbl exist
if (exist(fjp2,'file') ~= 2) || (exist(flbl,'file') ~= 2)
    error([fjp2, ' and it''s label do not exist.']);
end
lbl = read_pds_lbl(flbl);

% Mapping parameters
m = lbl.image_map_projection;
% Core
c = lbl.uncompressed_file.image;

% If we already have line and sample numbers we can skip to reading the
% subframe.

% Otherwise we need to calculate either the xy and pixel coords,
% or just the pixel coords.

% Mode is a logical vector that is true in element n when:
% n   modestr     condition
% 1   'lonlat'    Given lon-lat.
% 2   'xy'        Given xy in equidistant cylindrical map-projection of HiRISE image.
% 3   'pix'       Given pixel numbers of image subframe
% 4   'all'       Regardless of content of lims1 and lims2, read whole image.
% 5   ''          For empty string, do the same as 4, read whole image.
mode = [strcmpi(modestr, 'lonlat') strcmpi(modestr, 'xy')...
        strcmpi(modestr, 'pix')    strcmpi(modestr, 'all')...
        isempty(modestr)];

% Ensure ascending limits
lims1 = sort(lims1, 'ascend');
lims2 = sort(lims2, 'ascend');
    
% mode tells us how to interpret lims1 and lims2 in order to read HiRISE
% subframe.
if any(mode)
    
    % HiRISE images use the Mars sphere shape model, all radii are
    % 3396.190 km.
    R    = str1sttok2double(m.a_axis_radius)*1e3;
    lonp = str1sttok2double(m.center_longitude);
    latp = str1sttok2double(m.center_latitude);
    map_scale = str1sttok2double(m.map_scale);
    % Sample projection offset should be consistent with center lat and
    % center lon.
    s0 = str1sttok2double(m.sample_projection_offset);
    l0 = str1sttok2double(m.line_projection_offset);
    
    % Size of image
    minsample = str1sttok2double(m.sample_first_pixel);
    maxsample = str1sttok2double(m.sample_last_pixel);
    minline = str1sttok2double(m.line_first_pixel);
    maxline = str1sttok2double(m.line_last_pixel);
    
    % If lonlat or xy we must calculate line and sample numbers from xy
    if mode(2) || mode(1)
        
        % If lonlat we must first calculate xy.
        if mode(1)            
            % Assumed planetocentric coordinates, because HiRISE is.
            % Calculate xy of lon-lat limits.
            xlims = R.*(deg2rad( lims1 - lonp )) .* cosd(latp);
            ylims = R.*deg2rad( lims2 );        
        else
            % If mode == 2, xy then we already have coords.
            xlims = lims1;
            ylims = lims2;
        end
        
        % Now we have x and y lims in map coords, get line and sample
        % numbers.
        samplelims = round(   xlims./map_scale+s0+1  );
        linelims   = round( -(ylims./map_scale-l0+1) );
        
    elseif mode(4) || mode(5)
        
        % If modestr is empty or contains 'all' then we read whole image, 
        % so sample numbers are:
        samplelims = [ minsample maxsample ];
        linelims   = [ minline   maxline   ];
    end
    
    % If outside image bounds error and exit
    if any([samplelims(:); linelims(:)] < 1)
        warning('Subframe out of image bounds');
        lonvec=NaN; latvec=NaN; xvec=NaN; yvec=NaN; im=NaN; linelims=NaN; samplelims=NaN;
    else
        % Line numbers decrease with increasing y-coord, so reverse order of
        % limits.

        % Read subframe
        im = double(imread(fjp2, 'PixelRegion', {[linelims(2) linelims(1)], samplelims}));
        im(im == 0) = NaN;

        % Convert to science values
        offset = str2double(lbl.uncompressed_file.image.offset);
        scale = str2double(lbl.uncompressed_file.image.scaling_factor);
        im = (double(im) .* scale) + offset;

        % HiRISE pixels are aligned with latitude and longitude, and
        % map-projection axes, so we can return xvec, yvec,  latvec and lonvec
        % i.e., vectors. Not matrices.
        xvec = ((samplelims(1):samplelims(2)) - s0).*map_scale;
        % There is a mistake in the HiRISE SIS
        %yvec = (-l0 -(linelims(1):-1:linelims(2)) ).*map_scale;
        % This appears to work instead, l0 need not have its sign changed.
        yvec = (l0-(linelims(2):linelims(1))).*map_scale;

        % Now calculate latlon from  xvec yvec.
        %latvec = 180/pi*(yvec./R);
        %lonvec = lonp + xvec./(R*cosd(latp));
        [latvec, lonvec] = equirec2latlon(xvec, yvec, 0, 0, R, latp, lonp);
    end
else
	error('Unrecognised mode. Valid modes are: lonlat, xy, pix'); 
end