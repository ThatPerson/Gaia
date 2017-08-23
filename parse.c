extern int diploid;
extern int sexual;
extern int sim_length;
extern char output_folder[500];

extern float radiation_intensity;
extern float global_temperature;
extern int resources_for_reproducing;

extern int initial_mutation_rate;
extern float initial_colour;
extern float initial_t_opt;
extern int initial_dispersal;
extern int initial_progeny;
extern int initial_pop;

extern int verbose;
extern int verbose_s;
extern int cheat;
extern int goldschmidt;
extern int goldschmidt_freq;
extern int cheat_freq;
extern int runprint;
extern float dominance;
extern float sex_freq;

extern float oscillation_wavelength;
extern int edea;

float get_arg(char * c, int start) {

	int i;
	char str[256] = "";
	for (i = start; i < strlen(c); i++) {
		str[i-start] = c[i];
	}

	return atof(str);
}

int parse_settings(int argc, char*argv[]) {
	int i;
	for (i = 1; i < argc; i++) {

		switch (argv[i][0]) {
			case 'd': diploid = 1; break;
			case 'h': diploid = 0; break;
			case 's': sexual = 1; break;
			case 'e': edea = 1; break; // 
			case 'm': dominance = (float) get_arg(argv[i], 1); break; //
			case 'a': sexual = 0; break;
			case 'c': cheat = 1; cheat_freq = (int) get_arg(argv[i], 1); break;
			case 'l': sim_length = (int) get_arg(argv[i], 1); break;
			case 'v': verbose = 1; verbose_s = (int) get_arg(argv[i], 1); break; //
			case 'u': runprint = 0; break; //
			case 'o': 
				switch (argv[i][1]) {
					case 's': sex_freq = (float) get_arg(argv[i], 2); break;
					default:
						if (i < argc - 1) {
							strcpy(output_folder, argv[i+1]);
							i++;
						} else {
							printf("Please provide a folder name.");
						}
						break;
				}
				break;
			case 't': global_temperature = get_arg(argv[i], 1); break;
			case 'r': resources_for_reproducing = (int) get_arg(argv[i], 1); break;
			case 'k': radiation_intensity = get_arg(argv[i], 1); break;
			case 'g': goldschmidt = 1; goldschmidt_freq = (int) get_arg(argv[i], 1); break;
			case 'f': oscillation_wavelength = get_arg(argv[i], 1); break; 
			case 'i':
				switch (argv[i][1]) {
					case 'm': initial_mutation_rate = (int) (get_arg(argv[i], 2) * 1000); break;
					case 'c': initial_colour = get_arg(argv[i], 2); break;
					case 't': initial_t_opt = get_arg(argv[i], 2); break;
					case 'd': initial_dispersal = (int) get_arg(argv[i], 2); break;
					case 'p': initial_pop = (int) get_arg(argv[i], 2); break;
					case 'k': initial_progeny = (int) get_arg(argv[i], 2); break;
					default: break;
				}
				break;
			case '!': 
				printf("Welcome to Gaia, a Daisyworld simulation program.\ncNUMBER - mate individuals across the map to simulate human intervention. Number is times in 10,000. \ngNUMBER - introduce goldschmidt macromutations. Number is frequency in 10,000.\nd - create diploid daisies (DEFAULT)\nh - create haploid daisies\ne - activate EDEA\nmNUMBER - set dominance coefficient (DEFAULT = 0.5)\nvNUMBER - enable verbose mode, print detailed map every NUMBER times\nu - disable output\ns - create sexual daisies (DEFAULT)\na - create asexual daisies\nlNUMBER - set number of timesteps to run (DEFAULT = 6000)\noFOLDER - set output folder (DEFAULT = .)\ntNUMBER - set initial temperature\nrNUMBER - set resources needed to reproduce (DEFAULT = 30)\nkNUMBER - set radiation intensity (DEFAULT = 1.0)\nimNUMBER - set initial mutation rate (DEFAULT 0.05)\nicNUMBER - set initial colour (DEFAULT 0.5)\nitNUMBER - set initial optimum temperature\nidNUMBER - set initial dispersal (DEFAULT 3)\nikNUMBER - set initial number of progeny (DEFAULT 2)\nipNUMBER - set initial population (DEFAULT 10)\n");
				return -1;
				break;
			default: break;
		}
			
	}

	//printf("%d\n", 1+ atoi(argv[1]));
	if (sexual == 1 && diploid == 0) {
		printf("Can't have sexual haploids.\n");
		return -1;
	}
	return 1;
}


