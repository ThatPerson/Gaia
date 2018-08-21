/*----------------------------------------------------------------------*
 * File:    gaia.c                                                      *
 *                                                                      *
 * Purpose: Simulate a 2-D cellular automata model "DEAsy-World".       *
 *                                                                      *
 *                                                                      *
 * Author:  Ben Tatman                                                  *
 *          University of Cambridge                                     *
 *          ben@tatmans.co.uk                                           *
 *                                                                      *
 * License: Copyright 2017-2018 Ben Tatman                              *
 * Permission is hereby granted, free of charge, to any person          *
 * obtaining a copy of this software and associated documentation files *
 * (the "Software"), to deal in the Software without restriction,       *
 * including without limitation the rights to use, copy, modify, merge, *
 * publish, distribute, sublicense, and/or sell copies of the Software, *
 * and to permit persons to whom the Software is furnished to do so,    *
 * subject to the following conditions:                                 *
 *                                                                      *
 * The above copyright notice and this permission notice shall be       *
 * included in all copies or substantial portions of the Software.      *
 *                                                                      *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,      *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF   *
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND                *
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS  *
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN   *
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN    *
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE     *
 * SOFTWARE.                                                            *
 *                                                                      *
 * History: 21-Aug-2018, version 1.0                                    *
 *----------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "parse.c"

/*----------------------------------------------------------------------*
 * Constants                                                            *
 *----------------------------------------------------------------------*/
#define LANDSCAPE_X 50						// LANDSCAPE_X and Y define the dimensions of the map
#define LANDSCAPE_Y 50
#define CARRYING_CAP 10000					// Carrying capacity implements an artificial maximum carrying capacity.
#define MAX_PROGENY 10						// The number of progeny a daisy can have can mutate - this limits that.
#define PI 3.141592654

/*----------------------------------------------------------------------*
 * Global variables                                                     *
 *----------------------------------------------------------------------*/
int initial_mutation_rate = 0.05 * 1000; 	// This is the mutation rate (eg % of time it will mutate) multiplied by 1000. Multiplication occurs because later on this is compared to rng(1,1000).
float initial_colour = 0.5;					// Now we define some initial characteristics; colour,
float initial_t_opt = 25;					// optimum temperature,
int initial_dispersal =	3;					// dispersal, (this is the radius around which baby daisies can be placed
int initial_progeny = 2;					// progeny, (this is the number of baby daisies produced)
int initial_pop = 10;						// population, (the number of daisies placed at the start.
int cheat = 0;								// In earlier tests we experimented with adding cheats - that is to say, a daisy at [20, 20] can mate with one at [40, 40]. We hypothesised that this might lead to the development of cheats (eg a daisy which is dark and so warming and likes warming environments, and so leads to a runaway effect). This is left in here.
int edea = 0;								// How are the alleles expressed? edea = 1 produces differential allelic expression, =0 produces random expression.
int peak_verb = 0;							// If this is set to 1 then it tries to output verbose data when it reaches a peak in the oscillation.

//float mutation_deviation[7] = {0.05, 0.25, 0, 0.0, 0, 0.25, 0.05};
//float mutation_deviation[7] = {0.01, 0.001, 0, 0, 0, 0.001, 0.01};
float mutation_deviation[7] = {0.05, 0.1, 0, 0, 0, 0.1, 0.05};
											// These define the steps in mutations -
											// 			[0] - colour mutation in allele A
											//			[1] - temperature mutation in allele A
											//			[2] - dispersal mutation (integer!)
											//			[3] - progeny mutation (integer!)
											//			[4] - mutation rate mutation
											//			[5] - temperature mutation in allele B
											//			[6] - colour mutation in allele B

float goldschmidt_mm[7] = {0.5, 4, 0, 0, 0, 4, 0.5};
											// These are macromutations. Earlier in modeling we were experimenting with adding lower frequency large mutations.
int goldschmidt_freq = 0; 					// The frequency of macromutations. This is times per 10,000.
int age_of_death = 50;						// The age at which daisies will die naturally.
float mean_radiation_intensity = 1;
float radiation_intensity = 1;				// The radiation intensity experienced.
int filled_positions_x[CARRYING_CAP * 1000];
int filled_positions_y[CARRYING_CAP * 1000];
int filled_length = 0;
float global_temperature = 25;				// The average temperature over the whole grid.
int resources_for_reproducing = 30;			// How many resources are required for a daisy to reproduce.

int diploid = 1; 							// 1 if diploid, 0 if haploid.
int sexual = 1; 							// 1 if sexual, 0 if asexual. Can't have a sexual haploid.
int sim_length = 60;						// How long the simulation runs for.
char output_folder[500] = ".";				// Defines where the model can find the output/ folder
int verbose = 0;							// Activates verbose mode whereby it will output a map and temperature_map (AFTER RUN()) every verbose_s steps
int verbose_s = 20;							// How regularly the map and temperature_map are expressed.
int goldschmidt = 0;						// If this is 1 then the model will try to enable macromutations (see above).
int cheat_freq = 0;							// (see int cheat) How frequently cheating events occur.
FILE *vertical, *horizontal;
float oscillation_wavelength = 0;			// This is generally not used in favour of defining the equation below.
int runprint = 1;							// This outputs "Running (TIME) RAD" every time.
float sex_freq = 0; 						// This allows for asexual diploids to have low frequency sex.

struct Daisy {
	int pos_x, pos_y, dispersal, progeny, age, generation, living, mutation_rate, cheat, switchs, current;
	float colour[2], t_opt[2], local_te, cumulated_resources;
};											// The general daisy structure. switchs and current are relating to which temperature allele is being expressed.

struct Daisy daisies[CARRYING_CAP*1000];	// This holds the daisies.

struct Daisy * daisy_map[LANDSCAPE_X][LANDSCAPE_Y];
											// This is an array of pointers to the daisy present at position [x][y].

int num_daisies = 0;						// This keeps track of the operating length of daisies[]
int num_alive = 0;							// At certain points in the simulation daisies[] may contain dead daisies, which will still be counted by num_daisies. This variable counts the number of those which are alive.
int time_q;
int period = 200;

float temperature_map[LANDSCAPE_X][LANDSCAPE_Y];
											// The equivalent of daisy_map, only this holds the temperatures over the map.

char filename[500];

float sigmaconstant;						// 5.67*10^-8, this is defined below.
float solar_intensity = 1366; // Wm^-2

float dominance = 0.5;						// This affects how t_opt is expressed by changing the relative proportions of each temperature allele.

/*----------------------------------------------------------------------*
 * Function: rng                                                        *
 * Purpose:  Generates a random number between n and m                  *
 * Params:   n = lower bound                                            *
 *           m = upper bound                                            *
 * Returns:  Random number from n to m inclusive                        *
 *----------------------------------------------------------------------*/
int rng(int n, int m) {
	return n + rand()%(m-n+1); 				// rand()%(m-n+1) gives a random number from 0 to m-n, so adding n gives a random number from n to m.
}

/*----------------------------------------------------------------------*
 * Function: t_opt                                                      *
 * Purpose:  Calculates the optimum temperature for a given daisy       *
 * Params:   d = pointer to a Daisy                                     *
 *           temperature = local temperature                            *
 * Returns:  Optimum temperature of daisy                               *
 *----------------------------------------------------------------------*/
float t_opt(struct Daisy * d, float temperature) {
	if (edea == 1) {
		float a = (dominance * d->t_opt[0]) + ((1-dominance) * d->t_opt[1]);
											// We calculate the possible growth temperatures
		float b = (dominance * d->t_opt[1]) + ((1-dominance) * d->t_opt[0]);
											// Then we calculate the absolute difference between these growth temperatures and the local temperature.
		float delta_a = fabs(a-temperature);
		float delta_b = fabs(b-temperature);
		if (delta_a < delta_b) {			// If A is closer to the local temperature than B then we want to express A.
			if (d->current == 1)
					d->switchs = 1;
			if (dominance > 0.5)
				d->current = 1;
			else
				d->current = 0;

			return a;
		}
		if (d->current == 0)
			d->switchs = 1;
		if (dominance > 0.5)
			d->current = 0;
		else
			d->current = 1;
		return b;							// Else we want to express B.
	} else {
		return d->t_opt[rng(0,1)];			// If DAE is disabled then we express the temperature alleles at random
	}
}

/*----------------------------------------------------------------------*
 * Function: colour                                                     *
 * Purpose:  Calculates the expressed colour of a daisy                 *
 * Params:   d = pointer to a Daisy                                     *
 * Returns:  Current colour of daisy                                    *
 *----------------------------------------------------------------------*/
float colour(struct Daisy *d) {
	return dominance * d->colour[0] + (1-dominance) * d->colour[1]; // Colour is not expressed with differential allelic expression, but it still takes into account the dominance.
}

/*----------------------------------------------------------------------*
 * Function: radiation_factor                                           *
 * Purpose:  Calculates the radiation multiplier at a position on the   *
 *           map.                                                       *
 * Params:   n = number varying from 0 - 1 from top to bottom of map    *
 * Returns:  Multiplier                                                 *
 *----------------------------------------------------------------------*/
float radiation_factor(float n) {
	return 0.8 + 0.4*n;						// 0.8 at pole, 1.2 at equator.
}

/*----------------------------------------------------------------------*
 * Function: nutrient_gradient                                          *
 * Purpose:  Calculates the nutrient multiplier at a position on the    *
 *           map.                                                       *
 * Params:   n = number varying from 0 - 1 from left to right of map    *
 * Returns:  Multiplier                                                 *
 *----------------------------------------------------------------------*/
float nutrient_gradient(float n) {
	return 1;								// This is similar to radiation_factor only it operates east to west.
}

/* Runs each time step and varies the overall radiation intensity */
void vary_ri(void) {
	radiation_intensity += 0;				// This has generally been superseeded by just having the equation in main().
}



/*----------------------------------------------------------------------*
 * Function: check_pos                                                  *
 * Purpose:  Checks a position to see if a daisy is present             *
 * Params:   x = the x coordinate of input position                     *
 *           y = the x coordinate of input position                     *
 * Returns:  0 if the position is vacant, and                           *
             1 if the position is filled                                *
 *----------------------------------------------------------------------*/
int check_pos(int x, int y) {
	if (daisy_map[x][y] == NULL)
		return 0;
	return 1;
}

/*----------------------------------------------------------------------*
 * Function: mate                                                       *
 * Purpose:  Produces sexual progeny of two daisies.                    *
 * Params:   p = array of pointers to parental daisies.                 *
 *           progenitors = pointer to array of progeny                  *
 *           v = unused                                                 *
 * Returns:  Number of progeny produced                                 *
 *----------------------------------------------------------------------*/
int mate(struct Daisy *p[2], struct Daisy *progenitors, int v) {
	assert(p);
	assert(progenitors);
	int i, randoms[18], pq, current = 0;
	for (i = 0; i < 18; i++) {
		randoms[i] = rng(0,1);				// This is just an array of 0s and 1s at random, which are used to determine from which parent each allele is from.
	}
	float total_resources = p[0]->cumulated_resources + p[1]->cumulated_resources;
											// The progeny have access to both parents resources.
											// Then we remove these resources from the parents.


	for (i = 0; i < p[randoms[0]]->progeny; i++) {
		int y_delta = rng(0, 2*p[randoms[1]]->dispersal) - p[randoms[1]]->dispersal;
		int x_delta = pow(-1, rng(0, 1)) * rng(0, (int) sqrt(pow(p[randoms[1]]->dispersal, 2) - pow(y_delta, 2)));
											// y delta is a number from -dispersal to +dispersal. x delta then uses Pythagoras' equation a^2 + b^2 = c^2 to determine the relating x value (eg for dispersal = 3, when y delta = 0 x delta = 3, and visa versa.
		float deltas[7] = {0, 0, 0, 0, 0, 0, 0};
		for (pq = 0; pq < 7; pq++) {
			if (rng(0, 1000) < p[randoms[2]]->mutation_rate)
				deltas[pq] = pow(-1, rng(0, 1)) * mutation_deviation[pq];
		}									// deltas[] is an array laid out the same way as mutation_deviation which refers to the actual change.

		if (goldschmidt == 1) {
			// float goldschmidt_mm[7] = {0.5, 4, 0, 0, 0, 4, 0.5};
			// int goldschmidt_freq = 1; // times in 10000
			for (pq = 0; pq < 7; pq++) {
				if (rng(0, 10000) <= goldschmidt_freq)
					deltas[pq] += pow(-1, rng(0, 1)) * goldschmidt_mm[pq]; // Macromutations as mentioned before.
			}
		}
											// What follows is the definition of the new alleles according to their current alleles (randomly chosen from each parent's set of alleles) and the modifier from deltas[].
		float new_ca = p[0]->colour[randoms[12]] + deltas[0];
		float new_cb = p[1]->colour[randoms[11]] + deltas[6];
		if (new_ca > 1)
			new_ca = 1;
		if (new_cb > 1)
			new_cb = 1;

		float new_tea = p[0]->t_opt[randoms[4]] + deltas[1];
		float new_teb = p[1]->t_opt[randoms[5]] + deltas[5];
		int new_d = (int) p[randoms[6]]->dispersal + deltas[2];
		int new_p = (int) p[randoms[7]]->progeny + deltas[3];
		if (new_p > MAX_PROGENY)
			new_p = MAX_PROGENY;
		int new_m = (int) p[randoms[8]]->mutation_rate + deltas[4];
		int new_x = p[randoms[9]]->pos_x + x_delta;
		int new_y = p[randoms[10]]->pos_y + y_delta;

		if (check_pos(new_x, new_y) != 1 && new_x < LANDSCAPE_X && new_x >= 0 && new_y >= 0 && new_y < LANDSCAPE_Y) {
											// Is the position vacant? If so, is it in bounds?

			progenitors[current].pos_x = new_x;
			progenitors[current].pos_y = new_y;
			progenitors[current].t_opt[0] = (new_tea>0)?new_tea:0;
			progenitors[current].t_opt[1] = (new_teb>0)?new_teb:0;
			progenitors[current].dispersal = new_d;
			progenitors[current].progeny = (new_p > 0)?new_p:0;
			progenitors[current].age = 0;
			progenitors[current].local_te = 0;
			progenitors[current].generation = p[randoms[11]]->generation + 1;
			progenitors[current].mutation_rate = new_m;
			progenitors[current].living = 1;
			progenitors[current].switchs = 0;
			progenitors[current].cheat = (p[randoms[9]]->cheat == 1)?2:p[randoms[9]]->cheat;
			progenitors[current].cumulated_resources = total_resources / (p[randoms[0]]->progeny + 1);											// The total resources are divided between all progeny.
			progenitors[current].colour[0] = ((new_ca > 0)?new_ca:0);
			progenitors[current].colour[1] = ((new_cb > 0)?new_cb:0);
			filled_positions_x[filled_length] = new_x;
			filled_positions_y[filled_length] = new_y;
			daisy_map[new_x][new_y] = &progenitors[current];
			if (v == 1) {
				printf("\n\n\nP {(%f %f) (%f %f)} (%f %f) x {(%f %f) (%f %f)} (%f %f)\n", p[0]->t_opt[0], p[0]->t_opt[1], p[0]->colour[0], p[0]->colour[1], t_opt(p[0], temperature_map[p[0]->pos_x][p[0]->pos_y]), colour(p[0]), p[1]->t_opt[0], p[1]->t_opt[1], p[1]->colour[0], p[1]->colour[1], t_opt(p[1], temperature_map[p[1]->pos_x][p[1]->pos_y]), colour(p[1]));
				printf(" F1 {(%f %f) (%f %f)} -> %f, %f DELTAS %f %f\n", progenitors[current].t_opt[0], progenitors[current].t_opt[1], progenitors[current].colour[0], progenitors[current].colour[1], t_opt(&progenitors[current], temperature_map[new_x][new_y]), colour(&progenitors[current]), deltas[1], deltas[5]);
			}
			filled_length++;
			current = current + 1;
		}


	}
	return current;
}

/*----------------------------------------------------------------------*
 * Function: s_reproduce                                                *
 * Purpose:  Locates a mate and reproduces sexually.                    *
 * Params:   d = pointer to a daisy                                     *
 *           progenitors = pointer to array of progeny                  *
 * Returns:  Number of progeny produced                                 *
 *----------------------------------------------------------------------*/
int s_reproduce(struct Daisy*d, struct Daisy*progenitors) {
	assert(d);
	assert(progenitors);
	int min_x, max_x, min_y, max_y, current = 0;	// These define the range of potential mates.
	min_x = d->pos_x-1;
	max_x = d->pos_x+1;
	min_y = d->pos_y-1;
	max_y = d->pos_y+1;
	if (min_x < 0)
		min_x = 0;
	if (max_x > LANDSCAPE_X-1)
		max_x = LANDSCAPE_X-1;
	if (min_y < 0)
		min_y = 0;
	if (max_y > LANDSCAPE_Y-1)
		max_y = LANDSCAPE_Y-1;
	int x, y;
	struct Daisy * p[2];
	p[0] = d;
	int found = 0;
	//printf("%d %d\n", d->age, p[0]->age);
	for (x = min_x; x <= max_x; x++) {
		for (y = min_y; y <= max_y; y++) {
			if (x != d->pos_x && y != d->pos_y && daisy_map[x][y] != NULL) {
				if (daisy_map[x][y]->living == 1 && daisy_map[x][y]->cumulated_resources > resources_for_reproducing) {
					p[1] = daisy_map[x][y];
					found = 1;
				}
			}
		}
	}
	if (found == 0) {
		//printf("No partners!");
		return 0;
	}
		//printf("%d %d\n", p[1]->age, p[0]->age);

	// The parents are now p[0] and p[1]. Fill in progeny as before, only crossing over the genes.

	current = mate(p, progenitors, 0);

	return current;
}

/*----------------------------------------------------------------------*
 * Function: reproduce                                                  *
 * Purpose:  Reproduces clonally                                        *
 * Params:   d = pointer to a daisy                                     *
 *           progenitors = pointer to array of progeny                  *
 * Returns:  Number of progeny produced                                 *
 *----------------------------------------------------------------------*/
int reproduce(struct Daisy * d, struct Daisy * progenitors) {
	assert(progenitors);

	int i, p, pq, current = 0 ;
	for (i = 0; i < d->progeny; i++) {
		int y_delta = rng(0, 2*d->dispersal) - d->dispersal;
		int x_delta = pow(-1, rng(0, 1)) * rng(0, (int) sqrt(pow(d->dispersal, 2) - pow(y_delta, 2)));
		float deltas[7] = {0, 0, 0, 0, 0, 0, 0};
		for (p = 0; p < 7; p++) {
			if (rng(0, 1000) < d->mutation_rate)
				deltas[p] = pow(-1, rng(0, 1)) * mutation_deviation[p];

		}

		if (goldschmidt == 1) {
			// float goldschmidt_mm[7] = {0.5, 4, 0, 0, 0, 4, 0.5};
			// int goldschmidt_freq = 1; // times in 10000
			for (pq = 0; pq < 7; pq++) {
				if (rng(0, 10000) <= goldschmidt_freq)
					deltas[pq] += pow(-1, rng(0, 1)) * goldschmidt_mm[pq];
			}
		}
											// Same as before only for one parent.
		float new_ca = d->colour[0] + deltas[0];
		float new_cb = d->colour[1] + deltas[6];
		if (new_ca > 1)
			new_ca = 1;
		if (new_cb > 1)
			new_cb = 1;
		float new_tea = d->t_opt[0] + deltas[1];
		float new_teb = d->t_opt[1] + deltas[5];

		if (diploid == 0) {
			new_tea = new_teb;
			new_ca = new_cb;
		}

		int new_d = (int) d->dispersal + deltas[2];
		int new_p = (int) d->progeny + deltas[3];
		if (new_p > MAX_PROGENY) {
			new_p = MAX_PROGENY;
		}
		int new_m = (int) d->mutation_rate + deltas[4];
		int new_x = d->pos_x + x_delta;
		int new_y = d->pos_y + y_delta;

		if (check_pos(new_x, new_y) != 1 && new_x < LANDSCAPE_X && new_x > 0 && new_y > 0 && new_y < LANDSCAPE_Y) {
			progenitors[current].pos_x = new_x;
			progenitors[current].pos_y = new_y;
			progenitors[current].t_opt[0] = (new_tea>0)?new_tea:0;
			progenitors[current].t_opt[1] = (new_teb>0)?new_teb:0;
			progenitors[current].dispersal = new_d;
			progenitors[current].progeny = (new_p > 0)?new_p:0;
			progenitors[current].age = 0;
			progenitors[current].local_te = 0;
			progenitors[current].generation = d->generation + 1;
			progenitors[current].mutation_rate = new_m;
			progenitors[current].living = 1;
			progenitors[current].switchs = 0;
			progenitors[current].cheat = (d->cheat == 1)?2:d->cheat;
			progenitors[current].cumulated_resources = d->cumulated_resources / (d->progeny + 1);
			progenitors[current].colour[0] = ((new_ca > 0)?new_ca:0);
			progenitors[current].colour[1] = ((new_cb > 0)?new_cb:0);
			filled_positions_x[filled_length] = new_x;
			filled_positions_y[filled_length] = new_y;
			daisy_map[new_x][new_y] = &progenitors[current];
			filled_length++;
			current = current + 1;
		}
	}
	return current;
}

/*----------------------------------------------------------------------*
 * Function: update_t_map                                               *
 * Purpose:  Smoothes temperature landscape and calculated non occupied *
 *           space temperatures.                                        *
 * Params:   none                                                       *
 * Returns:  none                                                       *
 *----------------------------------------------------------------------*/
void update_t_map(void) {
	int x, y, xi, yi, min_x, max_x, min_y, max_y;
	float b_temperature_map[LANDSCAPE_X][LANDSCAPE_Y];
	float albedo_temp;
	for (x = 0; x < LANDSCAPE_X; x++) {
		for (y = 0; y < LANDSCAPE_Y; y++) {
			if (check_pos(x, y) == 0) {			// Then temperature uses base albedo of 0.5
				albedo_temp = pow((solar_intensity * radiation_factor((float) y/LANDSCAPE_Y) * radiation_intensity * (1-0.5))/(4*sigmaconstant), 1/4.);
				albedo_temp = albedo_temp + 25 - 234;
				temperature_map[x][y] = 0.7 * albedo_temp + 0.3*(temperature_map[x][y]);
			}

			b_temperature_map[x][y] = temperature_map[x][y];
												// Used so we don't use new values in calculations.
		}
	}

	// If a position doesn't have anything growing on it we need to set its temp.

	for (x = 0; x < LANDSCAPE_X; x++) {
		for (y = 0; y < LANDSCAPE_Y; y++) {
			min_x = x - 1;
			max_x = x + 2;
			min_y = y - 1;
			max_y = y + 2;
			if (min_x < 0)
				min_x = 0;
			if (max_x > LANDSCAPE_X - 1)
				max_x = LANDSCAPE_X - 1;
			if (min_y < 0)
				min_y = 0;
			if (max_y > LANDSCAPE_Y - 1)
				max_y = LANDSCAPE_Y - 1;
			float total = 0;
			for (xi = min_x; xi < max_x; xi++) {
				for (yi = min_y; yi < max_y; yi++) {
					total += b_temperature_map[xi][yi];
				}
			}
			temperature_map[x][y] = total / ((max_y - min_y) * (max_x - min_x));
		}
	}
	return;
}

/*----------------------------------------------------------------------*
 * Function: grow                                                       *
 * Purpose:  Grows the daisy, accumulates resources                     *
 * Params:   d = pointer to a daisy                                     *
 * Returns:  0                                                          *
 *----------------------------------------------------------------------*/
int grow(struct Daisy * daisy) {
	assert(daisy);

	if (daisy->living == 0) {
		return -1;
	}

	daisy->age++;
	if (daisy->age > age_of_death || daisy->cumulated_resources < 0 || check_pos(daisy->pos_x, daisy->pos_y) != 1) {						// We kill off old daisies.
		if (daisy->living == 1) {
			daisy->living = 0;
			daisy_map[daisy->pos_x][daisy->pos_y] = NULL;
			num_alive--;					// We haven't yet cleaned up the daisy array (and can't in this context) so we can't remove from num_daisies.
		}
		return -1;
	}

	float albedo_temp = pow((solar_intensity * radiation_factor((float) daisy->pos_y/LANDSCAPE_Y) * radiation_intensity * (1-colour(daisy)))/(4*sigmaconstant), 1/4.);
	// Stefan-Boltzmann Equation to calculate the daisy temperature */
	albedo_temp = albedo_temp + 25 - 234;
	daisy->local_te = 0.7*albedo_temp + 0.3*(temperature_map[daisy->pos_x][daisy->pos_y]);
	/* 70% of the local temperature is due to the daisy, 30% is due to the thermal insulation of the ground */
	temperature_map[daisy->pos_x][daisy->pos_y] = daisy->local_te;

	float delta_resources = (5 - pow(daisy->local_te - t_opt(daisy, daisy->local_te), 2)) * nutrient_gradient((float) daisy->pos_x/LANDSCAPE_X);
	daisy->cumulated_resources += delta_resources;

	return 0;
}

/*----------------------------------------------------------------------*
 * Function: run                                                        *
 * Purpose:  Grows all daisies, sorts out temperature, and outputs data *
 * Params:   n = current timestep                                       *
 * Returns:  none                                                       *
 *----------------------------------------------------------------------*/
void run(int n) {
	vary_ri();
	float planetary_albedo = 0;
	int i, p;
										// We prepare the output data for vertical.csv
	for (i = 0; i < num_daisies; i++) {
		planetary_albedo += (daisies[i].living == 1)?colour(&daisies[i]):0;
	}

	planetary_albedo += 0.5 * (LANDSCAPE_Y * LANDSCAPE_X - num_alive);
	planetary_albedo = planetary_albedo / (LANDSCAPE_X * LANDSCAPE_Y); // Average albedo
	float total = 0;
	int x, y;
	for (x = 0; x < LANDSCAPE_X; x++) {
		for (y = 0; y < LANDSCAPE_Y; y++) {
			total+=temperature_map[x][y];
		}
	}
	global_temperature = total/(LANDSCAPE_X*LANDSCAPE_Y); // Average temperature
	float sd_global_temp = 0;
	for (x = 0; x < LANDSCAPE_X; x++) {
		for (y = 0; y < LANDSCAPE_Y; y++) {
			sd_global_temp += pow((temperature_map[x][y] - global_temperature), 2);
		}
	}
	sd_global_temp = sqrt(sd_global_temp / (LANDSCAPE_X * LANDSCAPE_Y));

	/* Output data */
	float n_t_opta[2] = {0, 0};
	float n_t_optb[2] = {0, 0};
	float n_colour[2] = {0, 0};
	int n_progeny[2] = {0, 0};
	int n_dispersal[2] = {0, 0};
	int n_mutation_rate[2] = {0, 0};
	int n_count[2] = {0, 0};
	int min_y = LANDSCAPE_Y;
	int max_y = 0;
	int white = 0, black = 0, grey = 0, num_cheat = 0, switching = 0;
	float average_t_opt = 0;
	for (i = 0; i < num_daisies; i++) {
		if (daisies[i].living == 1) {
			int c = (int) round((float) daisies[i].pos_y / LANDSCAPE_Y);


			if (colour(&daisies[i]) > 0.55) {
				white++;
				c = 1;
			} else if (colour(&daisies[i]) < 0.45) {
				black++;
				c = 0;
			} else {
				grey++;
				c = 0;
			}
			n_count[c]++;
			n_t_opta[c]+=daisies[i].t_opt[0];
			n_t_optb[c]+=daisies[i].t_opt[1];
			average_t_opt += t_opt(&daisies[i], temperature_map[daisies[i].pos_x][daisies[i].pos_y]);
			n_colour[c]+=colour(&daisies[i]);
			n_progeny[c]+=daisies[i].progeny;
			if (daisies[i].switchs == 1)
				switching++;
			if (daisies[i].cheat == 1)
				num_cheat++;
			n_dispersal[c]+=daisies[i].dispersal;
			n_mutation_rate[c]+=daisies[i].mutation_rate;

			if (min_y > daisies[i].pos_y)
				min_y = daisies[i].pos_y;
			if (max_y < daisies[i].pos_y)
				max_y = daisies[i].pos_y;



		}
	}

	float divergence = 0;
	float sd_divergence = 0;
	float q = 0;
	int num_pop = 0;
	for (i = 0; i < num_daisies; i++) {
		if (daisies[i].living == 1) {
			q = fabs(daisies[i].t_opt[0] - daisies[i].t_opt[1]);
			divergence += q;
			num_pop++;
		}

	}

	divergence = divergence/num_pop;

	for (i = 0; i < num_daisies; i++) {
		if (daisies[i].living == 1) {
			sd_divergence += pow(fabs(daisies[i].t_opt[0] - daisies[i].t_opt[1]) - divergence, 2);
		}
	}
	sd_divergence = sqrt(sd_divergence/num_pop);

	/* Reproduce */
	for (i = 0; i < num_daisies; i++) {
		if (daisies[i].living == 1) {
			grow(&daisies[i]);
		}
	}
								// And now we run the simulation...
	struct Daisy * sp[2];
	struct Daisy new_daisies[MAX_PROGENY];
	average_t_opt = average_t_opt/num_pop;

	for (i = 0; i < num_daisies; i++) {
		if (daisies[i].living == 1 && daisies[i].cumulated_resources > resources_for_reproducing) {
			int length = 0;
			if (sexual == 1) {
				if (cheat == 1 && rng(1, 10000) < cheat_freq) {
					sp[0] = &daisies[i];
					sp[1] = &daisies[rng(0, num_daisies)];
					printf("Cheating...\n");
					length = mate(sp, new_daisies, 1);
					for (p = 0; p < length; p++) {
						new_daisies[p].cheat = 1;
					}
				} else {
					length = s_reproduce(&daisies[i], new_daisies);
				}
			} else {
				if (rng(0, 10000) < sex_freq * 10000) {
					// Sexual
					length = s_reproduce(&daisies[i], new_daisies);
				} else {
					length = reproduce(&daisies[i], new_daisies);
				}
			}
			for (p = 0; p < length; p++) {
				daisies[num_daisies] = new_daisies[p];
				daisy_map[daisies[num_daisies].pos_x][daisies[num_daisies].pos_y] = &daisies[num_daisies];
				num_daisies++;
				num_alive++;
			}
		}
	}




	if (n_count[0] == 0)
		n_count[0] = 1000000000;
	if (n_count[1] == 0)
		n_count[1] = 1000000000;
	/* Output to file */
	fprintf(vertical, "%d, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %d, %d, %d, %d, %d, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %d, %d, %d, %d, %d, %0.2f, %0.2f, %0.2f, %0.2f\n",
		n, global_temperature, (float) n_t_opta[0]/n_count[0], (float) n_t_optb[0]/n_count[0], (float)n_t_opta[1]/n_count[1], (float)n_t_optb[1]/n_count[1], (float)n_colour[0]/n_count[0], (float)n_colour[1]/n_count[1], (n_count[0] != 1000000000)?n_count[0]:0, (n_count[1] != 1000000000)?n_count[1]:0, min_y, max_y, num_pop, (float) n_progeny[0]/n_count[0], (float) n_progeny[1]/n_count[1], (float) n_dispersal[0]/n_count[0], (float) n_dispersal[1]/n_count[1], (float) n_mutation_rate[0]/n_count[0], (float) n_mutation_rate[1]/n_count[1], radiation_intensity, white, black, grey, num_cheat, switching, sd_global_temp, divergence, sd_divergence, average_t_opt);
												// The data from above is output.

	/* Cull the number of daisies down to the carrying capacity */
	int number_to_cull = num_alive - CARRYING_CAP;
	for (i = 0; i < number_to_cull; i++) {
		int q = rng(0, num_daisies);
		if (daisies[q].living == 1) {
			daisies[q].living = 0;
			daisy_map[daisies[q].pos_x][daisies[q].pos_y] = NULL;
			num_alive--;
		}
	}
	return;
}


/*----------------------------------------------------------------------*
 * Function: setup                                                      *
 * Purpose:  Sets up initial map.                                       *
 * Params:   none                                                       *
 * Returns:  none                                                       *
 *----------------------------------------------------------------------*/
void setup(void) {
	int i;
	int x, y;
	for (x = 0; x < LANDSCAPE_X; x++) {
		for (y = 0; y < LANDSCAPE_Y; y++) {
			temperature_map[x][y] = 25;

		}
	}
										// We use predefined initial values.
	for (i = 0; i < initial_pop; i++) {
		int x = rng(0, LANDSCAPE_X);
		int y = rng((LANDSCAPE_Y/2)-1, (LANDSCAPE_Y/2));
		if (check_pos(x, y) == 0) {
			daisies[num_daisies].pos_x = x;
			daisies[num_daisies].pos_y = y; // check pos 1 if something there
			filled_positions_x[filled_length] = x;
			filled_positions_y[filled_length] = y;
			filled_length++;
			daisies[num_daisies].colour[0] = initial_colour;
			daisies[num_daisies].colour[1] = initial_colour;
			daisies[num_daisies].t_opt[0] = initial_t_opt;
			daisies[num_daisies].t_opt[1] = initial_t_opt;
			daisies[num_daisies].dispersal = initial_dispersal;
			daisies[num_daisies].progeny = initial_progeny;
			daisies[num_daisies].age = 0;
			daisies[num_daisies].local_te = 0;
			daisies[num_daisies].generation = 0;
			daisies[num_daisies].cheat = 0;
			daisies[num_daisies].mutation_rate = initial_mutation_rate;
			daisies[num_daisies].living = 1;
			daisies[num_daisies].cumulated_resources = rng(0, 5);
			daisies[num_daisies].switchs = 0;
			daisy_map[x][y] = &daisies[num_daisies];
			num_daisies++;
			num_alive++;
		} else {
			i--;
		}
	}
}

/*----------------------------------------------------------------------*
 * Function: output_map                                                 *
 * Purpose:  Outputs pd*csv map                                         *
 * Params:   n = current timestep                                       *
 * Returns:  none                                                       *
 *----------------------------------------------------------------------*/
void output_map(int n) {
	char filename[500];
	sprintf(filename, "output/pd%05d.csv", n);

	FILE * locations;
	locations = fopen(filename, "w+");
	int i;
	for (i = 0; i < num_daisies; i++) {
		if (daisies[i].living == 1) {
			fprintf(locations, "%d, %d, %0.2f, %0.2f, %d, %0.2f, %0.2f, %d, %0.2f, %0.2f, %d\n", daisies[i].pos_x, daisies[i].pos_y, colour(&daisies[i]), daisies[i].local_te, daisies[i].progeny, t_opt(&daisies[i], daisies[i].local_te), temperature_map[daisies[i].pos_x][daisies[i].pos_y], daisies[i].cheat, daisies[i].t_opt[0], daisies[i].t_opt[1], daisies[i].switchs);
		}
	}
	fclose(locations);
	return;
}

/*----------------------------------------------------------------------*
 * Function: output_tmap                                                *
 * Purpose:  Outputs td*csv temperature map                             *
 * Params:   n = current timestep                                       *
 * Returns:  none                                                       *
 *----------------------------------------------------------------------*/
void output_tmap(int n) {
	char filename[500];
	sprintf(filename, "output/td%05d.csv", n);

	FILE * locations;
	locations = fopen(filename, "w+");
	int x, y;
	for (x = 0; x < LANDSCAPE_X; x++) {
		for (y = 0; y < LANDSCAPE_Y; y++) {
			fprintf(locations, "%d, %d, %f\n", x, y, temperature_map[x][y]);
		}
	}
	fclose(locations);
	return;
}

/*----------------------------------------------------------------------*
 * Function: main                                                       *
 * Purpose:  Entry point to program.                                    *
 * Params:   argc = number of arguments                                 *
 *           argv -> array of arguments                                 *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
int main(int argc, char* argv[]) {
	srand(time(NULL));
	sigmaconstant = 5.67*pow(10, -8);
	if (-1 == parse_settings(argc, argv))
		return 1;
	setup();

	/* Open the files required for constant output */
	vertical = fopen("output/vertical.csv",   "w+");
	horizontal = fopen("output/horizontal.csv", "w+");
	if (vertical == NULL || horizontal == NULL) {
		printf("File failed to open.\n");
	}

	fprintf(vertical, "n, global_temperature, (float) n_t_opta[0]/n_count[0], (float) n_t_optb[0]/n_count[0], (float)n_t_opta[1]/n_count[1], (float)n_t_optb[1]/n_count[1], (float)n_colour[0]/n_count[0], (float)n_colour[1]/n_count[1], (n_count[0] != 1000000000)?n_count[0]:0, (n_count[1] != 1000000000)?n_count[1]:0, min_y, max_y, num_alive, (float) n_progeny[0]/n_count[0], (float) n_progeny[1]/n_count[1], (float) n_dispersal[0]/n_count[0], (float) n_dispersal[1]/n_count[1], (float) n_mutation_rate[0]/n_count[0], (float) n_mutation_rate[1]/n_count[1], radiation_intensity, white, black, grey, cheat\n");
	int i, p, q, l;
	int pl = 0;
	int last = 0;
	float radiation_intensity_next, radiation_intensity_last;

	int going_down = 0, going_up = 0;

	for (i = 0; i < sim_length; i++) {

										// In some sexual models the daisies can be far apart from one another and so can't reproduce.
		if (i == 3 && num_alive == 0 && sexual == 1) {
			// None of the initial progeny are in the right place. Return -1 and start over.
			return -1;
		}
		if (runprint == 1)
			printf("Running %d (%0.8f)...\n", i, radiation_intensity);
		//radiation_intensity += 0.05*((PI * i / (100*pow((2000-(i/100)), 2))) + (PI / (2000-(i/100)))) * cos((PI * i)/(2000-(i/100)));


					// ######################
					// The following code calculates the radiation intensity. This equation needs to be varied to change the intensity series.
					// ######################

		float psdsd = (float)60000/200000;
		float qwe = 400 - (pl/36);
		//radiation_intensity = 1 + (psdsd * sin((3.14*(float) pl)/qwe));
		radiation_intensity_last = radiation_intensity;
		radiation_intensity = 1 + (psdsd * sin((3.14*(float) pl)/qwe));
		if (i > 500)
			pl++;

		if (peak_verb == 1) {
			radiation_intensity_next = 1 + (psdsd * sin((3.14*(float) (pl+1))/(400 - ((pl+1)/36))));
			going_down = 0;
			going_up = 0;

			if (radiation_intensity > radiation_intensity_next) {
				going_down = 1;
			}
			if (radiation_intensity > radiation_intensity_last) {
				going_up = 1;

			}
			if (going_down == 1 && going_up == 1 && (i - last) > 200) {
				output_map(i);
				output_tmap(i);
				last = i;
			}

		}


		if (oscillation_wavelength != 0)
			radiation_intensity += 0.00005 * sin(i * 3.14 / oscillation_wavelength);



							// Finally we get to running. We update the temperature map and then run().
		time_q++;
		update_t_map();
		run(i);
		if (i%verbose_s == 0 && verbose == 1) {
			output_map(i);
			output_tmap(i);
		}



		// Shrink array by removing dead daisies.
		l = num_daisies;
		int r_num_alive = 0;
		for (p = 0; p < l; p++) {
			if (daisies[p].switchs == 1)
				daisies[p].switchs = 0;
			if (daisies[p].living == 0) {
				for (q = p; q < num_daisies; q++) {
					daisy_map[daisies[q].pos_x][daisies[q].pos_y] = NULL;
					daisies[q] = daisies[q+1];
					daisy_map[daisies[q].pos_x][daisies[q].pos_y] = &daisies[q];
					filled_positions_x[q] = filled_positions_x[q+1];
					filled_positions_y[q] = filled_positions_y[q+1];
				}
				num_daisies--;
				filled_length--;
			} else {
				r_num_alive++;
			}
		}

		/* Clean daisy map */
		int x, y;
		for (x = 0; x < LANDSCAPE_X; x++) {
			for (y = 0; y < LANDSCAPE_Y; y++) {
				daisy_map[x][y] = NULL;
			}
		}

		for (x = 0; x < num_daisies+1; x++) {
			daisy_map[daisies[x].pos_x][daisies[x].pos_y] = &daisies[x];
		}

		num_alive = r_num_alive;
	}

	if (vertical != NULL)
		fclose(vertical);
	if (horizontal != NULL)
		fclose(horizontal);

	return 1;
}
