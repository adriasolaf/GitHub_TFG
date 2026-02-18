function [ jd2k, r, v, vd, va, rpga, dvga ] = MGA2_PGA2 ( ...
    planets, jd2k0, tofs )
%MGA2_PGA2 Multi-gravity assist trajectory
%   Performs a multi-gravity assist trajectory by means of:
%   - Powered Gravity Assist, 2D
%   - Lambert arcs
%   The algorithm relies on the patched conics method.
%
% Inputs:
%   planets: encounter planets sequence
%   jd2k0: sequence initial date [days from J2000]
%   tofs: times of flight between sequence objects [days]
%
% Outputs:
%   jd2k: time at each encounter [days from J2000]
%   r: heliocentric position of planet at each encounter [km]
%   v: heliocentric velocity of planet at each encounter [km/s]
%   vd: heliocentric departure velocity of s/c at each transfer [km/s]
%   va: heliocentric arrival velocity of s/c at each transfer [km/s]
%   dvga: Gravity assist DeltaV manoeure at each encounter [km/s]
%   rpga: Gravity assist periapsis radius at each encounter [km]
%
% Note: indices are considered as follows:
%   |
%   |-- Sequence: EMJ (Earth->Mars->Jupiter)
%   |
%   |-> Encounter i=1 (departure from Earth)
%   |   |-> jd2k(1), r(1,:), v(1,:)
%   |
%   |-> Transfer i=1 (transfer arc from Earth to Mars)
%   |   |-> vd(1,:), va(1,:)
%   |
%   |-> Encounter i=2 (GA with Mars)
%   |   |-> jd2k(2), r(2,:), v(2,:), dvga(1)
%   |
%   |-> Transfer i=2 (transfer arc from Mars to Jupiter)
%   |   |-> vd(2,:), va(2,:)
%   |
%   |-> Encounter i=3 (arrival to Jupiter)
%       |-> jd2k(3), r(3,:), v(3,:)
%
% Example:
%   <https://solarsystem.nasa.gov/multimedia/gallery/Voyager_Path.jpg>
%   seq = 'EJSUN'; t0 = -8169; tof = [688,778,1613,1309]; % Voyager 1
%   [ t, r, v, vd, va ] = MGA_PGA2 ( seq, t0, tof ); % Compute MGA
%   Voyager-1 c3: 105.5 (10.27km/s), Voyager-2 c3: 102.4 (10.12km/s)
%
% References:
%	[-] n/a
%
% See also:
%	Lambert, GA_PGA2
%
%David de la Torre Sangra
%January 2015

% Constants
mu = GetBodyProps('Sun'); % Standard gravitational parameter (Sun) [km3/s2]
days2secs = 86400; % Days to seconds

% Lambert configuration
mr = 0; % No multi-revolutions
lp = 0; % Short-period solutions

% Auxiliary magnitudes
lplanets = length(planets); % Number of planets
ltransfers = lplanets - 1; % Number of transfers

% Preallocate arrays
jd2k = zeros(lplanets,1);
r = zeros(lplanets,3);
v = zeros(lplanets,3);
vd = zeros(ltransfers,3);
va = zeros(ltransfers,3);
dvga = zeros(ltransfers-1,1);
rpga = zeros(ltransfers-1,1);

% Dates of encounters with planets
jd2k(1) = jd2k0; % Departure date
for i=2:lplanets
    jd2k(i) = jd2k(i-1) + tofs(i-1); % Encounter i
end

% Planet state vectors at each respective encounter
for i=1:lplanets
    [ r(i,:), v(i,:) ] = GetBodyICF ( planets{i}, jd2k(i), mu, 1 );
end

% Sequence of Lambert transfer arcs
for i=1:ltransfers % Iterate sequence of transfers

    % Compute transfer angle, accounting for prograde motion
    dnu = DeltaNu3 ( r(i,:), r(i+1,:), 1 ); % Transfer angle
    if dnu > pi, lw = 1; else, lw = 0; end % Long/short way tranfer

    % Compute Lambert arc from the current planet up to the next planet
    [ vd(i,:), va(i,:) ] = Lambert ( ...
        r(i,:), r(i+1,:), tofs(i)*days2secs, mu, lw, mr, lp );
end

% Sequence of gravity assists
for i=1:(ltransfers-1) % Iterate sequence of transfers

    % Compute planetocentric s/c velocity vectors at encounter
    % Point-mass approximation is used to favour computational performance
    vinfi = va(i+0,:) - v(i+1,:); % Hyperbolic excess velocity at arrival
    vinfo = vd(i+1,:) - v(i+1,:); % Hyperbolic excess velocity at departure
    vinfin = norm(vinfi); % Norm of vinfi
    vinfon = norm(vinfo); % Norm of vinfo

    % Turning angle [rad]
    delta = acos(dot(vinfi,vinfo) / (vinfin * vinfon));

    % Standard gravitational parameter of the booster planet [km3/s2]
    [ mu_planet, vmr ] = GetBodyProps ( planets{i+1} );

    % Powered gravity assist manoeuvre, 2D: PGA2_Rp tentative
    [ dvga(i), rpga(i) ] = GA_PGA2_Rp ( ...
        vinfin, vinfon, delta, mu_planet );

    % Powered gravity assist manoeuvre, 2D: PGA2_Vinfo if low rp
    vmr_safety = 1.05 * vmr; % Safety margin for low periapsis
    if rpga(i) < vmr_safety % GA periapsis too low
        rpga(i) = vmr_safety; % GA at lowest admisible periapsis
        [ dvga(i) ] = GA_PGA2_Vinfo ( ...
            vinfin, vinfon, delta, vmr_safety, mu_planet );
    end

end

end

