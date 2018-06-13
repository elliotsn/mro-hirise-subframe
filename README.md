# mro-hirise-subframe
% 
% readHiriseJp2Subframe.m
%
% Function to read a portion of a HiRISE map-projected JP2. Assumes the
% label file is stored in the same directory as the HiRISE JP2.
%
% Usage:
%
%  [lonvec, latvec, xvec, yvec, im, flbl, lbl, linelims, samplelims]  = ...
%    readHiriseJp2Subframe(fjp2, lims1, lims2, modestr);
%
% Inputs: 
%  
%   fjp2     Path to HiRISE map-projected JP2. Label with the same filename
%            stem should exist in the same directory.
%   lims1/2  1x2 vector, min and max bounds to read in the JP2 of
%            whichever coordinate is intended to be used. 
%   modestr  String identifier of the mode in which to operate.
%            
%            Coord. system                       modestr
%            -------------------------------------------
%            lon, lat                            'lonlat'
%            equidistant cylindrical x,y         'xy'
%            sample, line                        'pix'
%
% Outputs:
%
%   lonvec, latvec  Vectors of longitude and latitude for each pixel centre
%                   (in same geographic coordinate system as HiRISE image).
%   xvec,yvec       Vectors of equidistant cylindrical coordinates for each
%                   pixel centre.
%   im              Subframe pixel array.
%   flbl            Filepath to the label.
%   lbl             Structure containing the contents of the label.
%   line/samplelims Min. and max. pixel numbers for the subframe.
%
% Elliot Sefton-Nash
%
% 06/10/2017  Original
% 28/02/2018  Update to include modes of operation, provide either pixel
%             numbers, lat-lon or xy coords as boundaries of subframe.
