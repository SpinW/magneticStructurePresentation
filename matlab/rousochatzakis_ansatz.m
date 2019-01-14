function [S, k, n, name, pname, limit] = rousochatzakis_ansatz(absS, x)

if nargout <= 3
	k = [1/2 1/2 1];
	n = [0 0 1];

	cth = cos(x(1));
	sth = sin(x(1));
	cthp = cos(x(2));
	sthp = sin(x(2));
	ex = [1; 0; 0];
	ey = [0; 1; 0];
	ez = [0; 0; 1];
	ep = (ex + ey) / sqrt(2);
	em = (ex - ey) / sqrt(2);

	% The following assumes the following set-up:
	%   bfo.genlattice('lat_const', [8 8 6], 'angled', [90 90 90], 'spgr', 'P 4 b m');
	%   bfo.addatom('r', [0 0 0.5], 'S', 2.5, 'label', 'MFe3', 'color', 'gold');
	%   bfo.addatom('r', [1/3 1/3+0.5 0.5], 'S', 2.5, 'label', 'MFe3', 'color', 'black');
	S = zeros(3, 6);
	S(:,1) = absS(1) * (cth*ex - sth*ez);    % Sf
	S(:,2) = absS(2) * (cth*ey - sth*ez);    % Sa
	S(:,3) = absS(3) * (cthp*em + sthp*ez);  % Se
	S(:,4) = absS(3) * (cthp*em + sthp*ez);  % Sc
	S(:,5) = absS(5) * (-cthp*ep + sthp*ez); % Sd
	S(:,6) = absS(6) * (-cthp*ep + sthp*ez); % Sb
    % Note that the p_i phase factor in eq 2-5 of the paper is not
    % needed here because that is handled by SpinW when we put in
    % that k = [1/2 1/2 1]. We only need to define the moments in
    % the first unit cell.
else
    % provide the limits for the parameters
    name  = 'Rousochatzakis Ansatz for Cairo Lattice';
    % parameter names
	pname = {'Theta', 'Theta_prime'}
    % limits on input parameters
	limit = [0 0; pi pi];
    % garbage
    S = []; k = []; n = [];
end

