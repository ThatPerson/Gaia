/*----------------------------------------------------------------------*
 * File:	gaia.c													  *
 *																	  *
 * Purpose: Simulate a 2-D cellular automata model "DEAsy-World".	   *
 *																	  *
 *																	  *
 * Author:  Ben Tatman												  *
 *		  University of Cambridge									 *
 *		  ben@tatmans.co.uk										   *
 *																	  *
 * License: Copyright 2017-2018 Ben Tatman							  *
 * Permission is hereby granted, free of charge, to any person		  *
 * obtaining a copy of this software and associated documentation files *
 * (the "Software"), to deal in the Software without restriction,	   *
 * including without limitation the rights to use, copy, modify, merge, *
 * publish, distribute, sublicense, and/or sell copies of the Software, *
 * and to permit persons to whom the Software is furnished to do so,	*
 * subject to the following conditions:								 *
 *																	  *
 * The above copyright notice and this permission notice shall be	   *
 * included in all copies or substantial portions of the Software.	  *
 *																	  *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,	  *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF   *
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND				*
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS  *
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN   *
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN	*
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE	 *
 * SOFTWARE.															*
 *																	  *
 * History: 20-Dec-2017, version 1.0									*
 *			* Jul-2018, version 1.1										*
 *----------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "parse.c"

/*----------------------------------------------------------------------*
 * Constants															*
 *----------------------------------------------------------------------*/

#define LANDSCAPE_X 50
#define LANDSCAPE_Y 50
#define CARRYING_CAP 10000
#define MAX_PROGENY 10
#define PI 3.141592654

/*----------------------------------------------------------------------*
 * Global variables													 *
 *----------------------------------------------------------------------*/

int initial_mutation_rate = 0.05 * 1000; // Times in 1000
float initial_colour = 0.5;
float initial_t_opt = 25;
int initial_dispersal =	3;
int initial_progeny = 2;
int initial_pop = 10;
int cheat = 0;
int edea = 0;
int peak_verb = 0;

float mutation_deviation[7] = {0.05, 0.1, 0, 0, 0, 0.1, 0.05};
float temperature_mutation = 0.1; // 1 or 5 in mutdev.
float colour_mutation = 0.05; // Equiv to 0 or 6 in mutation_deviation
float goldschmidt_mm[7] = {0.5, 4, 0, 0, 0, 4, 0.5};
int goldschmidt_freq = 0; // times in 10000
int age_of_death = 50;
float mean_radiation_intensity = 1;
float radiation_intensity = 1;
int filled_positions_x[CARRYING_CAP * 1000];
int filled_positions_y[CARRYING_CAP * 1000];
int filled_length = 0;
float global_temperature = 25;
int resources_for_reproducing = 30;
int polyploid = 0; // Enabling this allows the ploidy level of diatoms to mutate, which affects the resources for reproducing and the t_opt/colour calculations. These are _obligate asexual DAE_.

/* So for polyploids.
 *  Colour becomes average of all colour alleles.
 *  t_opt has 100% contribution from nearest, 50% from next nearest, 25% from next nearest, working its way down.
 *  At very low frequency polyploid can evolve.
 */


int diploid = 1; // 1 if diploid, 0 if haploid.
int sexual = 1; // 1 if sexual, 0 if asexual. Can't have a sexual haploid.
int sim_length = 60;
char output_folder[500] = ".";
int verbose = 0;
int verbose_s = 20;
int goldschmidt = 0;
int cheat_freq = 0;
FILE *vertical, *horizontal;
float oscillation_wavelength = 800;
int runprint = 1;
float sex_freq = 0; // If we're an asexual diploid then we can have sex occasionally.

struct Daisy {
	int pos_x, pos_y, dispersal, progeny, age, generation, living, mutation_rate, cheat, switchs, current, ploidy;
	float colour[40], t_opt[40], local_te, cumulated_resources;
	char sheltered_load[40][11];
};

struct RecombinationPair {
	char a[2][11];
};

struct Daisy daisies[CARRYING_CAP*1000];

struct Daisy * daisy_map[LANDSCAPE_X][LANDSCAPE_Y];

int num_daisies = 0;
int num_alive = 0;
int time_q;
int period = 200;

float temperature_map[LANDSCAPE_X][LANDSCAPE_Y];

char filename[500];

float sigmaconstant;
float solar_intensity = 1366; // Wm^-2

float dominance = 0.5;
int radia = 0; // determines how the environment changes
int slm = 0; // sheltered load mutations. Times in 75000
float slmmod = 0.8;
int wthreeinc = 6000;
int recombenefit = 10;
/*----------------------------------------------------------------------*
 * Function: rng														*
 * Purpose:  Generates a random number between n and m				  *
 * Params:   n = lower bound											*
 *		   m = upper bound											*
 * Returns:  Random number from n to m inclusive						*
 *----------------------------------------------------------------------*/
int rng(int n, int m) {
	return n + rand()%(m-n+1);
}

/*----------------------------------------------------------------------*
 * Function: t_opt													  *
 * Purpose:  Calculates the optimum temperature for a given daisy	   *
 * Params:   d = pointer to a Daisy									 *
 *		   temperature = local temperature							*
 * Returns:  Optimum temperature of daisy							   *
 *----------------------------------------------------------------------*/
float t_opt(struct Daisy * d, float temperature) {
	/*int as, bs;
	if (edea == 1) {
		if (polyploid == 1) {
			// New calculation for t_opt goes here!. TODO

			// t_opts are in d->t_opt[0 -> ploidy]. Find smallest difference, and apply calculation below.

			as = curr_low;
			bs = curr_sec;

		} else {
			as = 0;
			bs = 1;
		}
		float a = (dominance * d->t_opt[as]) + ((1-dominance) * d->t_opt[bs]);
		float b = (dominance * d->t_opt[bs]) + ((1-dominance) * d->t_opt[as]);
		float delta_a = fabs(a-temperature);
		float delta_b = fabs(b-temperature);
		if (delta_a < delta_b) {
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
		return b;

	} else {
		return d->t_opt[rng(0,1)];
		float delta_a = fabs(d->t_opt[0] - temperature);
		float delta_b = fabs(d->t_opt[1] - temperature);
		if (delta_a > delta_b) {
			if (d->current == 0)
				d->switchs = 1;
			d->current = 1;
			return d->t_opt[1];
		} else {
			if (d->current == 1)
				d->switchs = 1;
			d->current = 0;
			return d->t_opt[0];
		}
	}*/

	if (edea == 1) {
		if (polyploid == 1) {


			int i;
			float curr_l = 100000;
			int returns;
			for (i = 0; i < d->ploidy; i++) {

				if (fabs(d->t_opt[i] - temperature) < curr_l) {
					returns = i;
					curr_l = fabs(d->t_opt[i] - temperature);
				}
			}
			return returns;
			return 25;

		} else {
			float a = (dominance * d->t_opt[0]) + ((1-dominance) * d->t_opt[1]);
			float b = (dominance * d->t_opt[1]) + ((1-dominance) * d->t_opt[0]);
			float delta_a = fabs(a-temperature);
			float delta_b = fabs(b-temperature);
			if (delta_a < delta_b) {
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
			return b;
		}
	} else {
		return d->t_opt[rng(0,1)];
		float delta_a = fabs(d->t_opt[0] - temperature);
		float delta_b = fabs(d->t_opt[1] - temperature);
		if (delta_a > delta_b) {
			if (d->current == 0)
				d->switchs = 1;
			d->current = 1;
			return d->t_opt[1];
		} else {
			if (d->current == 1)
				d->switchs = 1;
			d->current = 0;
			return d->t_opt[0];
		}
	}
}

/*----------------------------------------------------------------------*
 * Function: colour													 *
 * Purpose:  Calculates the expressed colour of a daisy				 *
 * Params:   d = pointer to a Daisy									 *
 * Returns:  Current colour of daisy									*
 *----------------------------------------------------------------------*/
float colour(struct Daisy *d) {
	if (polyploid == 1) {
		int i;
		float temp = 0;
		for (i = 0; i < d->ploidy; i++) {
			temp += d->colour[i];
		}
		return temp / (d->ploidy);
	} else {
		return dominance * d->colour[0] + (1-dominance) * d->colour[1];
	}
}

/*----------------------------------------------------------------------*
 * Function: radiation_factor										   *
 * Purpose:  Calculates the radiation multiplier at a position on the   *
 *		   map.													   *
 * Params:   n = number varying from 0 - 1 from top to bottom of map	*
 * Returns:  Multiplier												 *
 *----------------------------------------------------------------------*/
float radiation_factor(float n) {
	return 0.8 + 0.4*n;
}

/*----------------------------------------------------------------------*
 * Function: nutrient_gradient										  *
 * Purpose:  Calculates the nutrient multiplier at a position on the	*
 *		   map.													   *
 * Params:   n = number varying from 0 - 1 from left to right of map	*
 * Returns:  Multiplier												 *
 *----------------------------------------------------------------------*/
float nutrient_gradient(float n) {
	return 1;
}

/* Runs each time step and varies the overall radiation intensity */
void vary_ri(void) {
	radiation_intensity += 0;
}

/*----------------------------------------------------------------------*
 * Function: check_pos												  *
 * Purpose:  Checks a position to see if a daisy is present			 *
 * Params:   x = the x coordinate of input position					 *
 *		   y = the x coordinate of input position					 *
 * Returns:  0 if the position is vacant, and						   *
			 1 if the position is filled								*
 *----------------------------------------------------------------------*/
int check_pos(int x, int y) {
	if (daisy_map[x][y] == NULL)
		return 0;
	return 1;
}

/*----------------------------------------------------------------------*
 * Function: recombine													   *
 * Purpose:  Recombines.					*
 * Params:   n - Recombination Pair										 *
 * Returns:  new pair								 *
 *----------------------------------------------------------------------*/
struct RecombinationPair recombine(struct RecombinationPair n) {
	// Recombine n.a and n.b
	struct RecombinationPair m;
	int bp = rng(0, 9);
	int st = rng(0, 1);
	int ts = (st == 1)?0:1;
	//printf("%d %d, %d\n", st, ts, bp);

	int i;
	for (i = 0; i < 10; i++) {
		if (i == bp) {
			st = (st == 1)?0:1;
			ts = (ts == 1)?0:1;
		}
		m.a[st][i] = n.a[0][i];
		m.a[ts][i] = n.a[1][i];
	}
	m.a[0][10] = n.a[0][10];
	m.a[1][10] = n.a[1][10];

	return m;
}


/*----------------------------------------------------------------------*
 * Function: sl_score													   *
 * Purpose:  Calculates sheltered load score..					*
 * Params:   sheltered_load												 *
 * Returns:  score								 *
 *----------------------------------------------------------------------*/
int sl_score(char sheltered_load[11]) {
	int i;
	int score = 0;
	for (i = 0; i < 10; i++) {
		if (sheltered_load[i] != '0')
			score++;
	}
	return score;
}

/*----------------------------------------------------------------------*
 * Function: sl_mutate													   *
 * Purpose:  Mutates sheltered load					*
 * Params:   sheltered_load												 *
 * Returns:  sheltered_load								 *
 *----------------------------------------------------------------------*/

char * sl_mutate(char sheltered_load[11]) {
	int i = rng(0, 9);
	sheltered_load[i] = (rng(1, 50) <= 49)?'-':'0';
	return sheltered_load;
}


/*----------------------------------------------------------------------*
 * Function: mate													   *
 * Purpose:  Produces sexual progeny of two daisies.					*
 * Params:   p = array of pointers to parental daisies.				 *
 *		   progenitors = pointer to array of progeny				  *
 *		   v = unused												 *
 * Returns:  Number of progeny produced								 *
 *----------------------------------------------------------------------*/
int mate(struct Daisy *p[2], struct Daisy *progenitors, int v) {
	assert(p);
	assert(progenitors);
	int i, randoms[18], pq, current = 0;
	for (i = 0; i < 18; i++) {
		randoms[i] = rng(0,1);
	}
	float total_resources = p[0]->cumulated_resources + p[1]->cumulated_resources;
	for (i = 0; i < p[randoms[0]]->progeny; i++) {
		int y_delta = rng(0, 2*p[randoms[1]]->dispersal) - p[randoms[1]]->dispersal;
		int x_delta = pow(-1, rng(0, 1)) * rng(0, (int) sqrt(pow(p[randoms[1]]->dispersal, 2) - pow(y_delta, 2)));
		float deltas[8] = {0, 0, 0, 0, 0, 0, 0, 0};
		for (pq = 0; pq < 8; pq++) {
			if (rng(0, 1000) < p[randoms[2]]->mutation_rate)
				deltas[pq] = pow(-1, rng(0, 1)) * mutation_deviation[pq];
		}
		if (goldschmidt == 1) {
			// float goldschmidt_mm[7] = {0.5, 4, 0, 0, 0, 4, 0.5};
			// int goldschmidt_freq = 1; // times in 10000
			for (pq = 0; pq < 7; pq++) {
				if (rng(0, 10000) <= goldschmidt_freq)
					deltas[pq] += pow(-1, rng(0, 1)) * goldschmidt_mm[pq];
			}
		}
		float new_ca = p[0]->colour[randoms[12]] + deltas[0];
		float new_cb = p[1]->colour[randoms[11]] + deltas[6];
		if (new_ca > 1)
			new_ca = 1;
		if (new_cb > 1)
			new_cb = 1;

		float new_tea = p[0]->t_opt[randoms[4]] + deltas[1];
		float new_teb = p[1]->t_opt[randoms[5]] + deltas[5];

		char new_sla[11];
		char new_slb[11];
		//printf("%d %d Recombining. %d %d.. %s ;; %s :: %s ;; %s\n", p[0]->ploidy, p[1]->ploidy, randoms[4], (randoms[4]==1)?0:1, p[0]->sheltered_load[randoms[4]], p[1]->sheltered_load[randoms[4]], p[0]->sheltered_load[(randoms[4]==1)?0:1], p[1]->sheltered_load[(randoms[4]==1)?0:1]);
		struct RecombinationPair a;
		struct RecombinationPair b;
		strcpy(a.a[0], p[0]->sheltered_load[randoms[4]]);
		strcpy(a.a[1], p[1]->sheltered_load[randoms[4]]);
		strcpy(b.a[0], p[0]->sheltered_load[(randoms[4] == 1)?0:1]);
		strcpy(b.a[1], p[1]->sheltered_load[(randoms[4] == 1)?0:1]);
		//printf("A: %s ::: %s\n", a.a[0], a.a[1]);
		//printf("B: %s ::: %s\n", b.a[0], b.a[1]);
		//printf("REC\n");
		if (rng(1, (25000/slm)/recombenefit) == 35) {
		//	printf("RECOMBINING\n");
			a = recombine(a);
			b = recombine(b);
			//exit(0);
		}
		//printf("A: %s ::: %s\n", a.a[0], a.a[1]);
		//printf("B: %s ::: %s\n", b.a[0], b.a[1]);
		strcpy(new_sla, a.a[rng(0,1)]);
		strcpy(new_slb, b.a[rng(0,1)]);
		//printf("Recombined - %s, %s\n", new_sla, new_slb);

		int new_d = (int) p[randoms[6]]->dispersal + deltas[2];
		int new_p = (int) p[randoms[7]]->progeny + deltas[3];
		if (new_p > MAX_PROGENY)
			new_p = MAX_PROGENY;
		int new_m = (int) p[randoms[8]]->mutation_rate + deltas[4];
		int new_x = p[randoms[9]]->pos_x + x_delta;
		int new_y = p[randoms[10]]->pos_y + y_delta;



		/* Will only do sexual reproduction if they are diploid. This code won't work for polyploids!!!!  */

		int q = 0;
		if (check_pos(new_x, new_y) != 1 && new_x < LANDSCAPE_X && new_x >= 0 && new_y >= 0 && new_y < LANDSCAPE_Y) {
			if (polyploid == 1) { /* Ploidy will always be 2 but this allows it to double/halve*/
				int sc = rng(0, 100);
				if (sc == 5) {
					progenitors[current].ploidy = (p[randoms[10]]->ploidy >= 2)?(p[randoms[10]]->ploidy)/2:p[randoms[10]]->ploidy;
				} else if (sc == 12) {
					progenitors[current].ploidy = (p[randoms[10]]->ploidy <= 8)?2*(p[randoms[10]]->ploidy):p[randoms[10]]->ploidy;

				} else {
					progenitors[current].ploidy = p[randoms[10]]->ploidy;
				}
			} else {
				progenitors[current].ploidy = p[randoms[10]]->ploidy; // = 1 or 2 depending on haploid or diploid.
			}
			progenitors[current].pos_x = new_x;
			progenitors[current].pos_y = new_y;
			progenitors[current].t_opt[0] = (new_tea>0)?new_tea:0;
			progenitors[current].t_opt[1] = (new_teb>0)?new_teb:0;
			progenitors[current].colour[0] = ((new_ca > 0)?new_ca:0);
			progenitors[current].colour[1] = ((new_cb > 0)?new_cb:0);
			strcpy(progenitors[current].sheltered_load[0], new_sla);
			strcpy(progenitors[current].sheltered_load[1], new_slb);
			/* Biological relevance? */
			if (progenitors[current].ploidy > 2) {
				for (q = 0; q < floor(progenitors[current].ploidy); q++) {
					progenitors[current].t_opt[q*2] = (new_tea>0)?new_tea:0;
					progenitors[current].t_opt[(2*q)+1] = (new_teb>0)?new_teb:0;
					progenitors[current].colour[q*2] = ((new_ca > 0)?new_ca:0);
					progenitors[current].colour[(2*q) + 1] = ((new_cb > 0)?new_cb:0);
					strcpy(progenitors[current].sheltered_load[q*2], new_sla);
					strcpy(progenitors[current].sheltered_load[(2*q) + 1], new_slb);
				}
				//printf("%s %s %s %s\n", progenitors[current].sheltered_load[0], progenitors[current].sheltered_load[1], progenitors[current].sheltered_load[2], progenitors[current].sheltered_load[3]);
			}

			//printf("OUTPUTS %s :::: %s\n", progenitors[current].sheltered_load[0], progenitors[current].sheltered_load[1]);
			progenitors[current].dispersal = new_d;
			progenitors[current].progeny = (new_p > 0)?new_p:0;
			progenitors[current].age = 0;
			progenitors[current].local_te = 0;
			progenitors[current].generation = p[randoms[11]]->generation + 1;
			progenitors[current].mutation_rate = new_m;
			progenitors[current].living = 1;
			progenitors[current].switchs = 0;
			progenitors[current].cheat = (p[randoms[9]]->cheat == 1)?2:p[randoms[9]]->cheat;
			progenitors[current].cumulated_resources = total_resources / (p[randoms[0]]->progeny + 1);

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
 * Function: s_reproduce												*
 * Purpose:  Locates a mate and reproduces sexually.					*
 * Params:   d = pointer to a daisy									 *
 *		   progenitors = pointer to array of progeny				  *
 * Returns:  Number of progeny produced								 *
 *----------------------------------------------------------------------*/
int s_reproduce(struct Daisy*d, struct Daisy*progenitors) {
	assert(d);
	assert(progenitors);
	int min_x, max_x, min_y, max_y, current = 0;
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
					if (polyploid == 1) {
						if (daisy_map[x][y]->ploidy == 2) {
							p[1] = daisy_map[x][y];
							found = 1;
						}
					} else {
						p[1] = daisy_map[x][y];
						found = 1;
					}
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
	//printf("Not a problem in s-reproduce\n");
	current = mate(p, progenitors, 0);

	return current;
}

/*----------------------------------------------------------------------*
 * Function: reproduce												  *
 * Purpose:  Reproduces clonally										*
 * Params:   d = pointer to a daisy									 *
 *		   progenitors = pointer to array of progeny				  *
 * Returns:  Number of progeny produced								 *
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



		float new_colours[20];
		float new_t_opt[20];
		char new_sl[20][11];
		int lk;
		for (lk = 0; lk < d->ploidy; lk++) {
			new_colours[lk] = d->colour[lk] + (pow(-1, rng(0, 1)) * colour_mutation);
			if (new_colours[lk] > 1)
				new_colours[lk] = 1;
			if (new_colours[lk] < 0)
				new_colours[lk] = 0;
			new_t_opt[lk] = d->t_opt[lk] + (pow(-1, rng(0, 1)) * temperature_mutation);
			strcpy(new_sl[lk], d->sheltered_load[lk]); // This is where we would recombine but as the alleles [0] and [0] from both parents are the same parent there is nothing to recombine. At the moment recombination is only between the same n - so we treat t_opt[0], t_opt[1] etc as being on different chromosomes, which recombine with t_opt[0] and t_opt[1] respectively (eg no [0]->[1] recombination).
		}
		if (polyploid == 1 && d->ploidy < 10) {
			for (lk = d->ploidy; lk < 2*(d->ploidy); lk++) {

				new_colours[lk] = d->colour[lk - (d->ploidy)] + (pow(-1, rng(0, 1)) * colour_mutation);
				if (new_colours[lk] > 1)
					new_colours[lk] = 1;
				if (new_colours[lk] < 0)
					new_colours[lk] = 0;
				new_t_opt[lk] = d->t_opt[lk - (d->ploidy)] + (pow(-1, rng(0, 1)) * temperature_mutation);
				strcpy(new_sl[lk], d->sheltered_load[lk - (d->ploidy)]);
			}

		}

		if (diploid == 0) {
			new_colours[0] = new_colours[1];
			new_t_opt[0] = new_t_opt[1];
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

			/* Evolve ploidy. 1 in 20 increase ploidy level by multiplication, 1 in 20 decrease by division. Can't go above 16 or below 1.*/
			if (polyploid == 1) {
				int sc = rng(0, 100);
				if (sc == 5) {
					progenitors[current].ploidy = (d->ploidy >= 2)?(d->ploidy)/2:d->ploidy;
				} else if (sc == 12) {
					progenitors[current].ploidy = (d->ploidy <= 8)?2*(d->ploidy):d->ploidy;

				} else {
					progenitors[current].ploidy = d->ploidy;
				}
			} else {
				progenitors[current].ploidy = d->ploidy; // = 1 or 2 depending on haploid or diploid.
			}



			for (lk = 0; lk < progenitors[current].ploidy; lk++) {
				progenitors[current].t_opt[lk] = (new_t_opt[lk]>0)?new_t_opt[lk]:0;
				progenitors[current].colour[lk] = new_colours[lk];
				strcpy(progenitors[current].sheltered_load[lk], new_sl[lk]);
			}

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
 * Function: update_t_map											   *
 * Purpose:  Smoothes temperature landscape and calculated non occupied *
 *		   space temperatures.										*
 * Params:   none													   *
 * Returns:  none													   *
 *----------------------------------------------------------------------*/
void update_t_map(void) {
	int x, y, xi, yi, min_x, max_x, min_y, max_y;
	float b_temperature_map[LANDSCAPE_X][LANDSCAPE_Y];
	float albedo_temp;
	for (x = 0; x < LANDSCAPE_X; x++) {
		for (y = 0; y < LANDSCAPE_Y; y++) {
			if (check_pos(x, y) == 0) {
				albedo_temp = pow((solar_intensity * radiation_factor((float) y/LANDSCAPE_Y) * radiation_intensity * (1-0.5))/(4*sigmaconstant), 1/4.);
				albedo_temp = albedo_temp + 25 - 234;
				temperature_map[x][y] = 0.7 * albedo_temp + 0.3*(temperature_map[x][y]);
			}

			b_temperature_map[x][y] = temperature_map[x][y];
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
 * Function: grow													   *
 * Purpose:  Grows the daisy, accumulates resources					 *
 * Params:   d = pointer to a daisy									 *
 * Returns:  0														  *
 *----------------------------------------------------------------------*/
int grow(struct Daisy * daisy) {
	assert(daisy);

	if (daisy->living == 0) {
		return -1;
	}

	daisy->age++;
	if (daisy->age > age_of_death || daisy->cumulated_resources < 0 || check_pos(daisy->pos_x, daisy->pos_y) != 1) {
		if (daisy->living == 1) {
			daisy->living = 0;
			daisy_map[daisy->pos_x][daisy->pos_y] = NULL;
			num_alive--;
		}
		return -1;
	}

	float albedo_temp = pow((solar_intensity * radiation_factor((float) daisy->pos_y/LANDSCAPE_Y) * radiation_intensity * (1-colour(daisy)))/(4*sigmaconstant), 1/4.);
	// Stefan-Boltzmann Equation to calculate the daisy temperature */
	albedo_temp = albedo_temp + 25 - 234;
	daisy->local_te = 0.7*albedo_temp + 0.3*(temperature_map[daisy->pos_x][daisy->pos_y]);
	/* 70% of the local temperature is due to the daisy, 30% is due to the thermal insulation of the ground */
	temperature_map[daisy->pos_x][daisy->pos_y] = daisy->local_te;
	float opt_t = (polyploid == 1)?daisy->t_opt[(int) t_opt(daisy, daisy->local_te)]:t_opt(daisy, daisy->local_te);

	/* spontaneous mutations in sheltered load, only on non expressed alleles. */
	int i;

	float delta_resources;

	if (slm != 0 && polyploid == 1) {
		int ot = (int) t_opt(daisy, daisy->local_te);

		for (i = 0; i < daisy->ploidy; i++) {
			if (rng(1, 25000/slm) == 35) {
				// Then we mutate this one.
				//daisy->sheltered_load[i] += (rng(0, 5) == 1)?-1:1;
				//printf("ORIGINAL %s\n", daisy->sheltered_load[i]);
				strcpy(daisy->sheltered_load[i], sl_mutate(daisy->sheltered_load[i]));

				//printf("Mutated sheltered load, now %f\n", daisy->sheltered_load[i]);
			}
		}
		if (sl_score(daisy->sheltered_load[ot]) > 4) {
			//printf("%d\n", sl_score(daisy->sheltered_load[ot]));
			//printf("%s\n", daisy->sheltered_load[ot]);
		}

		delta_resources = (5 - pow(daisy->local_te - opt_t, 2) - (slmmod *sl_score(daisy->sheltered_load[ot]))) * nutrient_gradient((float) daisy->pos_x/LANDSCAPE_X);
	} else {
		delta_resources = (5 - pow(daisy->local_te - opt_t, 2)) * nutrient_gradient((float) daisy->pos_x/LANDSCAPE_X);
	}



	daisy->cumulated_resources += delta_resources;

	return 0;
}

/*----------------------------------------------------------------------*
 * Function: run														*
 * Purpose:  Grows all daisies, sorts out temperature, and outputs data *
 * Params:   n = current timestep									   *
 * Returns:  none													   *
 *----------------------------------------------------------------------*/
void run(int n) {

	vary_ri();

	float planetary_albedo = 0;
	int i, p;
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
	float average_ploidy = 0, ploidy_sd = 0;
	int n_c = 0;
	for (i = 0; i < num_daisies; i++) {
		if (daisies[i].living == 1) {
			n_c++;
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

			grow(&daisies[i]);

		}
	}




	/* Reproduce */
	struct Daisy * sp[2];
	struct Daisy new_daisies[MAX_PROGENY];
	average_t_opt = average_t_opt/n_c;

	float divergence = 0;
	float sd_divergence = 0;
	float average_sl = 0;
	float local_sl = 0;
	float q = 0;
	int l;
	n_c = 0;
	for (i = 0; i < num_daisies; i++) {
		if (daisies[i].living == 1){
			n_c++;
			q = fabs(daisies[i].t_opt[0] - daisies[i].t_opt[1]);
			divergence += q;

			average_ploidy += daisies[i].ploidy;
			for (l = 0; l < daisies[i].ploidy; l++) {
				local_sl += fabs(sl_score(daisies[i].sheltered_load[l]));
			}
			local_sl = local_sl / daisies[i].ploidy;
			average_sl += local_sl;
			local_sl = 0;
		}
	}
	average_ploidy = average_ploidy / n_c;
	divergence = divergence/n_c;
	average_sl = average_sl / n_c;
	int numberof[5] = {0, 0, 0, 0, 0};



	for (i = 0; i < num_daisies; i++) {

		if (daisies[i].living == 1) {
			numberof[(int) (log(daisies[i].ploidy) / log(2))]++;
			sd_divergence += pow(fabs(daisies[i].t_opt[0] - daisies[i].t_opt[1]) - divergence, 2);
			ploidy_sd += pow(fabs(daisies[i].ploidy - average_ploidy), 2);
		}
	}
	sd_divergence = sqrt(sd_divergence/n_c);
	ploidy_sd = sqrt(ploidy_sd/n_c);


	float rfr = resources_for_reproducing;


	for (i = 0; i < num_daisies; i++) {
		rfr = resources_for_reproducing;
		if (polyploid == 1)
			rfr += 2*daisies[i].ploidy; // Cost of polyploidy.
		if (daisies[i].living == 1 && daisies[i].cumulated_resources > rfr && daisies[i].age > 2) {
			int length = 0;
			if (polyploid == 0) {
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
			} else {
				//printf("%d\n", daisies[i].ploidy);
				if (daisies[i].ploidy == 2) {
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
	fprintf(vertical, "%d, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %d, %d, %d, %d, %d, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %d, %d, %d, %d, %d, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %d, %d, %d, %d, %d\n",
		n, global_temperature, (float) n_t_opta[0]/n_count[0], (float) n_t_optb[0]/n_count[0], (float)n_t_opta[1]/n_count[1], (float)n_t_optb[1]/n_count[1], (float)n_colour[0]/n_count[0], (float)n_colour[1]/n_count[1], (n_count[0] != 1000000000)?n_count[0]:0, (n_count[1] != 1000000000)?n_count[1]:0, min_y, max_y, num_alive, (float) n_progeny[0]/n_count[0], (float) n_progeny[1]/n_count[1], (float) n_dispersal[0]/n_count[0], (float) n_dispersal[1]/n_count[1], (float) n_mutation_rate[0]/n_count[0], (float) n_mutation_rate[1]/n_count[1], radiation_intensity, white, black, grey, num_cheat, switching, sd_global_temp, divergence, sd_divergence, average_t_opt, average_ploidy, ploidy_sd, average_sl, numberof[0], numberof[1], numberof[2], numberof[3], numberof[4]);


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
 * Function: setup													  *
 * Purpose:  Sets up initial map.									   *
 * Params:   none													   *
 * Returns:  none													   *
 *----------------------------------------------------------------------*/
void setup(void) {
	int i;
	int x, y;
	for (x = 0; x < LANDSCAPE_X; x++) {
		for (y = 0; y < LANDSCAPE_Y; y++) {
			temperature_map[x][y] = 25;

		}
	}
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
			strcpy(daisies[num_daisies].sheltered_load[0], "0000000000");
			strcpy(daisies[num_daisies].sheltered_load[1], "0000000000");
			daisies[num_daisies].dispersal = initial_dispersal;
			daisies[num_daisies].ploidy = diploid + 1;
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
 * Function: output_map												 *
 * Purpose:  Outputs pd*csv map										 *
 * Params:   n = current timestep									   *
 * Returns:  none													   *
 *----------------------------------------------------------------------*/
void output_map(int n) {
	char filename[500];
	sprintf(filename, "output/pd%05d.csv", n);

	FILE * locations;
	locations = fopen(filename, "w+");
	int i;
	for (i = 0; i < num_daisies; i++) {
		if (daisies[i].living == 1) {
			fprintf(locations, "%d, %d, %0.2f, %0.2f, %d, %0.2f, %0.2f, %d, %0.2f, %0.2f, %d, %d\n", daisies[i].pos_x, daisies[i].pos_y, colour(&daisies[i]), daisies[i].local_te, daisies[i].progeny, t_opt(&daisies[i], daisies[i].local_te), temperature_map[daisies[i].pos_x][daisies[i].pos_y], daisies[i].cheat, daisies[i].t_opt[0], daisies[i].t_opt[1], daisies[i].switchs, daisies[i].ploidy);
		}
	}
	fclose(locations);
	return;
}

/*----------------------------------------------------------------------*
 * Function: output_tmap												*
 * Purpose:  Outputs td*csv temperature map							 *
 * Params:   n = current timestep									   *
 * Returns:  none													   *
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
 * Function: main													   *
 * Purpose:  Entry point to program.									*
 * Params:   argc = number of arguments								 *
 *		   argv -> array of arguments								 *
 * Returns:  None.													  *
 *----------------------------------------------------------------------*/
int main(int argc, char* argv[]) {
	srand(time(NULL));
	if (-1 == parse_settings(argc, argv))
		return 1;
	sigmaconstant = 5.67*pow(10, -8);
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
	float zl = 0;
	for (i = 0; i < sim_length; i++) {
		if (i == 3 && num_alive == 0 && sexual == 1) {
			// None of the initial progeny are in the right place. Return -1 and start over.
			return -1;
		}
		if (runprint == 1)
			printf("Running %d (%0.8f)...\n", i, radiation_intensity);
		//radiation_intensity += 0.05*((PI * i / (100*pow((2000-(i/100)), 2))) + (PI / (2000-(i/100)))) * cos((PI * i)/(2000-(i/100)));

		float psdsd = (float)60000/200000;


		radiation_intensity_last = radiation_intensity;
		if (radia == 0) {
			float qwe = oscillation_wavelength - (pl/36);
			//radiation_intensity = 1 + (psdsd * sin((3.14*(float) pl)/qwe));
			radiation_intensity = 1 + (psdsd * sin((3.14*(float) pl)/qwe)) + zl;
			if (i > 500)
				pl++;
			if (i > 1000)
				zl += 0.000; // Global Warming
		} else if (radia == 1) { // Global Warming
			float qwe = oscillation_wavelength;
			//radiation_intensity = 1 + (psdsd * sin((3.14*(float) pl)/qwe));
			radiation_intensity = 1 + (psdsd * sin((3.14*(float) pl)/qwe)) + zl;
			if (i > 500)
				pl++;
			if (i > 1000)
				zl += 0.0003; // Global Warming
		} else if (radia == 2) {
			radiation_intensity = 1;
		} else if (radia == 3) {
			if (i > wthreeinc) {
				zl += 0.0003;
			}

			radiation_intensity = 1 + zl;
		}

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
