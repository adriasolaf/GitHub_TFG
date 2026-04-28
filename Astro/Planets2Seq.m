function [ seq ] = Planets2Seq ( planets )
%PLANETS2SEQ Planet strings to Planet sequence
%   Returns a list with the sequence from a given planet name list
%
% Inputs:
%   planets: planet strings cellarray {string}
%
% Outputs:
%   seq: planet sequence [string]
%
% Example:
%   seq = Planets2Seq ( {'Earth','Venus','Mars','Jupiter'} );
%
% References:
%	[-]
%
%David de la Torre Sangra
%August 2016

% Preallocate planet names array
seq = '';

% Iterate sequence
for i=1:length(planets)
    
    % Get planet string
    if strcmpi(planets{i},'Mercury'), s = 'm';
    elseif strcmpi(planets{i},'Venus'), s = 'V';
    elseif strcmpi(planets{i},'Earth'), s = 'E';
    elseif strcmpi(planets{i},'Moon'), s = 'M';
    elseif strcmpi(planets{i},'Mars'), s = 'M';
    elseif strcmpi(planets{i},'Jupiter'), s = 'J';
    elseif strcmpi(planets{i},'Saturn'), s = 'S';
    elseif strcmpi(planets{i},'Uranus'), s = 'U';
    elseif strcmpi(planets{i},'Neptune'), s = 'N';
    elseif strcmpi(planets{i},'Pluto'), s = 'P';
    elseif strcmpi(planets{i},'Didymos'), s = 'D';
    else % Planet not in list
        error('Planet ID "%s" not implemented',planets{i});
    end
    
    % Save planet ID into planet sequence array
    seq(i) = s;
    
end

end

