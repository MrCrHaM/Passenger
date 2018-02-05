We used the real populations of the northeastern counties to generate the experimental data.

The struct variable “cons” the population and shape files for the 129 counties in northeastern.

Vector “S” contains the 23 counties that we chose as abnormal counties.

“yy_1″ is the observed case number for each county as shown in our experiment.

“pop_vec” is the population of each county.

You are free to use the population data to create data of your own like this:

l_0 = 5e-5; lam_0 = pop_vec * l_0; l_1 = 4 * l_0; lam_1 = lam_0; lam_1(S) = pop_vec(S) * l_1; yy_0 = poissrnd(lam_0); % null hypothesis yy_1 = poissrnd(lam_1);

The following code can be used to draw map shaped figures:

fall = flipud(autumn(numel(cons))); minpop = min([cons.rate]); maxpop = max([cons.rate]); densityColors = makesymbolspec(‘Polygon’, {‘rate’, … [minpop maxpop], ‘FaceColor’, fall}); figure geoshow(cons, ‘DisplayType’, ‘polygon’,’SymbolSpec’, densityColors);

Here we are drawing the case incidence rate of each county where cons(i).rate = yy_1(i)/pop_vec(i); You can add additional fields to the cons struct variable and visualize the population or ground truth or your results.
