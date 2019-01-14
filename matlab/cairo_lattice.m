% SpinW Magnetic Structure tutorial 2, Horace/SpinW Workshop 2019
%
% The aim of this tutorial is to calculate the phase diagram of a Cairo 
% pentagon lattice.
% The theory is described in this paper:
%   Quantum magnetism on the Cairo pentagonal lattice, 
%   I. Rousochatzakis, A.M. Laeuchli, R. Moessner,
%   Phys. Rev. B 85, 104415 (2012), 
%   arXiv:1201.3079 https://arxiv.org/abs/1201.3079
% We will only carry out the first (classical) part of the calculations 
% described.
%
% The Cairo pentagon lattice is a geometrically frustrated magnetic
% lattice comprised of edge-shared irregular pentagons.
% The pentagons each have one side longer or shorter than the other
% four. Ions linked by the unequal sides have three-fold symmetry,
% whilst the other ions have 4-fold symmetry. There are two types
% of nearest neighbour bonds:
%   J33 between sites of 3-fold symmetry, along the unequal side
%   J43 between a 3-fold and a 4-fold site.
%
% We first set up the crystal structure, which is based on a
% simplified version of the crystal structure of Bi2Fe4O9 which is
% the only known physical realisation of the Cairo lattice.

cairo = spinw;
cairo.genlattice('lat_const', [8 8 6], 'angled', [90 90 90], 'spgr', 'P 4 b m');

% Define the 4-fold sites:
cairo.addatom('r', [0 0 0.5], 'S', 2.5, 'label', 'MFe3', 'color', 'gold');

% Now define the 3-fold sites. The position p determines whether the
% unequal bond is shorter (p > 0.363) or longer than the other bonds
p = 0.33;
cairo.addatom('r', [p+0.5 p 0.5], 'S', 2.5, 'label', 'MFe3', 'color', 'black');

% Generate the couplings and define the interactions
cairo.gencoupling()
% Because the J33 may be longer or shorter than J43 we have to 
% figure which is the case, since SpinW indexes bonds by length
if size(cairo.table('bond',1),1) > size(cairo.table('bond',2),1)
    j43_bond = 1;
    j33_bond = 2;
else
    j43_bond = 2;
    j33_bond = 1;
end
cairo.addmatrix('label', 'J33', 'value', 1, 'color', 'green')
cairo.addmatrix('label', 'J43', 'value', 1, 'color', 'white')
cairo.addcoupling('mat', 'J33', 'bond', j33_bond)
cairo.addcoupling('mat', 'J43', 'bond', j43_bond)

% Plots the resulting lattice - compare to Fig 1 of the paper
plot(cairo, 'range', [2 2 1])

%%
% The frustration in this lattice comes from the competition
% between J33 and J43. When J33 is large, the system forms
% a set of dimers, whilst when J43 is large it forms a
% connected network. Rousochatzakis et al. define a frustration
% index x = J43/J33 which we will also use:

x = 1;

J33 = 1;
J43 = x * J33;

% Substitute in the new values for J33 and J43
cairo.matparser('mat', {'J33', 'J43'}, 'param', [J33 J43]);

% According to the paper, the magnetic structure depends on x
% We will now use optmagstr to try to compute it
%
% There are actually 6 magnetic ions in the unit cell.
% The "orthogonal" structure has k = [1/2 1/2 1] whilst
% the "ferrimagnetic" stucture has k = [1 1 1]
%
% We first use the built-in gm_planar function.
ang_min = [0 0 0 0 0 0];
ang_max = [pi pi pi pi pi pi];
nm = [0 0];   % The normal is the a-b plane
% First with a propagation vector corresponding to the
% orthogonal structure
k_orth = [0.5 0.5 1];
res_orth = cairo.optmagstr('func', @gm_planar, ...
    'xmin', [ang_min k_orth nm], 'xmax', [ang_max k_orth nm])
% Then with that corresponding to the ferrimagnetic one
k_feri = [1 1 1];
res_feri = cairo.optmagstr('func', @gm_planar, ...
    'xmin', [ang_min k_feri nm], 'xmax', [ang_max k_feri nm])
% Print out the energy
opt_gs_energy = [res_orth.e res_feri.e]
% Plot the structure with the lowest energy
if res_orth.e < res_feri.e
    plot(res_orth.obj, 'range', [2 2 1])
else
    plot(res_feri.obj, 'range', [2 2 1])
end

% Now run this section (ctrl+enter) several times
% - what do you notice about the structure, and gs energy?
% Does the optmised structure match what you expect from
% the theory paper?

%%
% We could try to help things along by doing an optmagsteep
% gradient descent afterwards

x = 1;

J33 = 1;
J43 = x * J33;
cairo.matparser('mat', {'J33', 'J43'}, 'param', [J33 J43]);
ang_min = [0 0 0 0 0 0];
ang_max = [pi pi pi pi pi pi];
nm = [0 0];
k_orth = [0.5 0.5 1];
res_orth = cairo.optmagstr('func', @gm_planar, ...
    'xmin', [ang_min k_orth nm], 'xmax', [ang_max k_orth nm])

% Do a steepest descent here to try to nudge the spins into
% the optimum configuration
res_orth = cairo.optmagsteep()

k_feri = [1 1 1];
res_feri = cairo.optmagstr('func', @gm_planar, ...
    'xmin', [ang_min k_feri nm], 'xmax', [ang_max k_feri nm])
res_feri = cairo.optmagsteep()

% Print out the energy 
% - the energy is now a vector of length the number of
%   gradient steps
opt_gs_energy = [res_orth.e(end) res_feri.e(end)]
% Plot the structure with the lowest energy
if res_orth.e(end) < res_feri.e(end)
    plot(res_orth.obj, 'range', [2 2 1])
else
    plot(res_feri.obj, 'range', [2 2 1])
end

%%
% Hopefully you should be able to see that in the optimised
% structure for small x, the two 4-fold sites are orthogonal
% to each other, and the two sets of 3-fold sites, linked by
% differently oriented J33 bonds are also orthogonal to
% each other, whilst the phase angle between the 3- and 
% 4-fold sites are free, and has no effect on the ground 
% state energy. The phase angle between the spins and
% the lattice also does not affect the gs energy.
% In the theory paper, the authors fix this phase angle to
% be zero, and the 3-fold / 4-fold phase angle to be 45 deg
% in the ansatz defined in eqs 2-5.

% One other thing you hopefully may have noticed is that
% the "mixed phase" found by the authors is not stabilised
% in the code - because we force the spins to be coplanar
% We now try a different built-in function which handles
% non-coplanar spin structures (at the cost of more free
% parameters which makes finding a minimum harder).

x = 1.8;

J33 = 1;
J43 = x * J33;
cairo.matparser('mat', {'J33', 'J43'}, 'param', [J33 J43]);
a3_min = [0 0 0 0 0 0 0 0 0 0 0 0];
a3_max = [pi pi pi pi pi pi pi pi pi pi pi pi];
nm = [0 0];
k_orth = [0.5 0.5 1];
res_3D = cairo.optmagstr('func', @gm_spherical3d, ...
    'xmin', [a3_min k_orth nm], 'xmax', [a3_max k_orth nm])
res_3D = cairo.optmagsteep()
% Run the other optimisations to compare
ang_min = [0 0 0 0 0 0];
ang_max = [pi pi pi pi pi pi];
nm = [0 0];
k_orth = [0.5 0.5 1];
res_orth = cairo.optmagstr('func', @gm_planar, ...
    'xmin', [ang_min k_orth nm], 'xmax', [ang_max k_orth nm])
res_orth = cairo.optmagsteep()
k_feri = [1 1 1];
res_feri = cairo.optmagstr('func', @gm_planar, ...
    'xmin', [ang_min k_feri nm], 'xmax', [ang_max k_feri nm])
res_feri = cairo.optmagsteep()

opt_gs_energy = [res_orth.e(end) res_feri.e(end) res_3D.e(end)]

plot(res_3D.obj, 'range', [2 2 1])

%%
% The above unfortunately does not work well to find the
% ground state because there are too many parameters
% (all the different phase angles) to minimise.
% Luckily for the Cairo lattice case, we have a theoretical
% ansazt for the magnetic structure, which is programmed
% into the function rousochatzakis_ansatz. We can run that
% and compare it with the two planar cases
% It requires only two parameters, theta and theta'
% which are defined in the paper.
% Also look at the file rousochatzakis_ansatz.m and check that 
% you understand the code and how it relates eqs 2-5 in the paper

x = 2.1;

J33 = 1;
J43 = x * J33;
cairo.matparser('mat', {'J33', 'J43'}, 'param', [J33 J43]);

res_ans = cairo.optmagstr('func', @rousochatzakis_ansatz, ...
    'xmin', [0 0], 'xmax', [pi pi])

% Run the other optimisations to compare
ang_min = [0 0 0 0 0 0];
ang_max = [pi pi pi pi pi pi];
nm = [0 0];
k_orth = [0.5 0.5 1];
res_orth = cairo.optmagstr('func', @gm_planar, ...
    'xmin', [ang_min k_orth nm], 'xmax', [ang_max k_orth nm])
res_orth = cairo.optmagsteep()
k_feri = [1 1 1];
res_feri = cairo.optmagstr('func', @gm_planar, ...
    'xmin', [ang_min k_feri nm], 'xmax', [ang_max k_feri nm])
res_feri = cairo.optmagsteep()

opt_gs_energy = [res_orth.e(end) res_feri.e(end) res_ans.e]

% Plots the structure with the lowest energy
objs = [res_orth.obj res_feri.obj res_ans.obj];
[~, idm] = min(opt_gs_energy);
plot(objs(idm), 'range', [2 2 1])

%%
% Finally, we can try to reconstruct the phase diagram
% in figure 2 of the paper

J33 = 1;
x = 0.1:0.1:3
for ix = 1:numel(x)
    J43 = x(ix) * J33;
    cairo.matparser('mat', {'J33', 'J43'}, 'param', [J33 J43]);
    res = cairo.optmagstr('func', @rousochatzakis_ansatz, ...
        'xmin', [0 0], 'xmax', [pi pi])
    theta(ix) = res.x(1);
    theta_p(ix) = res.x(2);
end

figure; hold all;
% Plots the SpinW calculated angles
plot(x, sin(theta), 'ok');
plot(x, sin(theta_p), 'sr');
% Plots the theoretical angles from the paper
sth_theo = sqrt(2 - 4./(x.^2));
plot(x, sth_theo, '-k');
plot(x, (x./2) .* sth_theo, '-r');

% Hopefully you should get perfect agreement!
% (Note that the theoretical expression is only valid for x<2)
