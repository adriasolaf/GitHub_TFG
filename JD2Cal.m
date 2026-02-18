function [ year, month, day, hour, minute, second ] = JD2Cal ( jd, REF )
%JD2CAL Julian day to Calendar date
%   Compute the Calendar date form Julian date using the algorithm from [1]
%
% Inputs:
%   jd: Julian day
%   REF: JD reference ('JD','J2000','RJD','MJD','TJD','DJD'). Default JD
%   CAL: Calendar input ('Gregorian','Julian'). Default Gregorian
%
% Outputs:
%   Calendar: year, month, day, hour, minute, second
%
% Example:
%   [ Y, M, D, h, m, s ] = JD2Cal ( 6351.34, 'J2000' ); % J2000
%
% References:
%   [1] https://www.mathworks.com/matlabcentral/fileexchange/
%       39191-low-precision-ephemeris/content/gdate.m
%
%David de la Torre Sangra
%May 2016

% Get default inputs
if nargin < 2 || isempty(REF), REF = 'JD'; end

% Adjust for optional reference frame
switch REF
    case 'JD' % J2000
        % Do nothing
    case 'J2000' % J2000
        jd = jd + 2451545.0;
    case 'RJD' % Reduced JD
        jd = jd + 2400000.0;
    case 'MJD' % Modified JD, introduced by SAO in 1957
        jd = jd + 2400000.5;
    case 'TJD' % Truncated JD, introduced by NASA in 1979
        jd = jd + 2440000.5;
    case 'DJD' % Dublin JD, introduced by the IAU in 1955
        jd = jd + 2415020.0;
    otherwise % Reference not implemented
        warning(['Reference "',REF,'" not impemented. ',...
            'Options are: JD J2000 RJD MJD TJD DJD. ',...
            'Assuming Julian Day JD']);
end

% Convert JD to Calendar (Gregorian)
z = fix(jd + .5);
fday = jd + .5 - z;
if (fday < 0)
    fday = fday + 1;
    z = z - 1;
end
if (z < 2299161)
    a = z;
else
    alpha = floor((z - 1867216.25) / 36524.25);
    a = z + 1 + alpha - floor(alpha / 4);
end
b = a + 1524;
c = fix((b - 122.1) / 365.25);
d = fix(365.25 * c);
e = fix((b - d) / 30.6001);
day = b - d - fix(30.6001 * e) + fday;
if (e < 14)
    month = e - 1;
else
    month = e - 13;
end
if (month > 2)
    year = c - 4716;
else
    year = c - 4715;
end

% Get fractional part of the day (H:M:S)
fracday = day - fix(day);
day = fix(day);
dechours = fracday * 24;
hour = fix(dechours);
decminutes = (dechours - hour) * 60;
minute = fix(decminutes);
second = (decminutes - minute) * 60;

end

