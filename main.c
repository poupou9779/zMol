#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <float.h>

#define NUMBER_OF_ATOMS 112
#define LENGTH_ATOM_NAME_MAX 3

struct atom
{
    char name[LENGTH_ATOM_NAME_MAX+1];
    double molar_mass;
};

double get_molar_mass_atom(const struct atom list_of_atoms[], size_t length_of_list, char *atom);
int extract_formula(double *coefficient, char **formula, const char *chemical_formula, size_t length_of_chemical);
double get_molar_mass(const struct atom list_of_atoms[], size_t length_of_list, const char *buffer, size_t length_of_buffer);

/*
in :
    - list_of_atoms is a table of struct atom which contains the atoms you want to check
    - length_of_list contains the size of list-of_atoms (there are 133 atoms in the Periodic table, but you can choose to selct a few of those)
    - atom is a string which contains the chemical name of the atom
out :
    - -1 if atom was not found in list_of_atoms
    - the molar mass written in list_of_atoms of atom
*/
double get_molar_mass_atom(const struct atom list_of_atoms[], size_t length_of_list, char *atom)
{
    size_t i;
    for(i = 0; i < length_of_list; ++i)
        if(strcmp(atom, list_of_atoms[i].name) == 0)
            return list_of_atoms[i].molar_mass;
    return -1.;
}

/*
in :
    - coefficient is a pointer on the coefficient value (pointer because the function changes its value)
    - formula is a double-pointer on char because it is allocated in that function (then the value of the pointer changes)
    - chemical_formula is the buffer which contains the formula nd the coefficient
    - length_of_chemical is the size of the buffer
out :
    - -1 if allocation failed
    - 0 if allocation successed
*/
int extract_formula(double *coefficient, char **formula, const char *chemical_formula, size_t length_of_chemical)
{
    char *tmp;
    int i = 0;
    while(isalpha(chemical_formula[i]) == 0 && chemical_formula[i] != '(')
        ++i;
    if(i != 0)
    {
        *coefficient = strtod(chemical_formula, &tmp);
        if(*coefficient == 0)
            *coefficient = 1;
        while(isalpha(*tmp) == 0)
            ++tmp;
        *formula = malloc(strlen(tmp)+1);
        if(*formula == NULL)
            return -1;
        strcpy(*formula, tmp);
    }
    else
    {
        *formula = malloc(length_of_chemical+1);
        if(*formula == NULL)
            return -1;
        strncpy(*formula, chemical_formula, length_of_chemical+1);
        (*formula)[length_of_chemical] = '\0';
    }
    return 0;
}

/*
in :
    - list_of_atoms is a table of struct atom which contains the atoms you want to check
    - length_of_list contains the size of list-of_atoms (there are 133 atoms in the Periodic table, but you can choose to selct a few of those)
    - buffer is a string which contains the formula of the molecule you want get the molar mass of
    - length_of_buffer contains the size of the buffer (useful to treat the brackets '()')
out :
    - -1.0 if there has been an allocation problem,
    - -2.0 if an atom is unknown (or doesn't exist)
    - -1500.0 if there is a non-alphanum and non-brace char.
    - the molar mass (precision .01) of the molecule if no error has been found
*/
double get_molar_mass(const struct atom list_of_atoms[], size_t length_of_list, const char *buffer, size_t length_of_buffer)
{
    double ret = 0.;
    char *tmp = NULL;
    char *stat_buf = NULL;
    int i = 0,
        tmp_index = LENGTH_ATOM_NAME_MAX-1,
        j;
    int pow_10 = 1;
    double tmp_molar_mass;
    unsigned long tmp_coef = 0;
    double coefficient = 1.;

    if(extract_formula(&coefficient, &stat_buf, buffer, length_of_buffer) == -1)
        return -1.;

    tmp = malloc(LENGTH_ATOM_NAME_MAX+1);
    tmp[LENGTH_ATOM_NAME_MAX] = '\0';
    for(i = strlen(stat_buf)-1; i >= 0; --i)
    {
        if(isdigit(stat_buf[i]) != 0)
        {
            tmp_coef += (stat_buf[i] - '0') * pow_10;
            pow_10 *= 10;
        }
        else if(isalpha(stat_buf[i]) != 0)
        {
            if(tmp_index < 0)
            {
                free(tmp);
                free(stat_buf);
                return -2.0;
            }
            else
            {
                tmp[tmp_index--] = stat_buf[i];
                if(isupper(tmp[tmp_index+1])) {
                    tmp_molar_mass = get_molar_mass_atom(list_of_atoms, length_of_list, tmp+tmp_index+1);
                    if(tmp_molar_mass < DBL_EPSILON)
                    {
                        free(tmp);
                        free(stat_buf);
                        return -2.0;
                    }
                    if(tmp_coef == 0)
                        tmp_coef = 1;
                    ret += tmp_molar_mass*tmp_coef;
                    tmp_coef = 0;
                    pow_10 = 1;
                    tmp_index = 2;
                }
            }
        }
        else if(stat_buf[i] == ')')
        {
            j = i-1;
            while(stat_buf[j] != '(')
                --j;
            ret += tmp_coef*get_molar_mass(list_of_atoms, length_of_list, &stat_buf[j+1], i-j-1);
            i = j;
            tmp_coef = 0;
            pow_10 = 1;
            tmp_index = 2;
        }
        else
        {
            ret = -1500.;
            break;
        }
    }
    free(tmp);
    free(stat_buf);
    return coefficient*ret;
}

int main(void)
{
    struct atom list_of_atoms[NUMBER_OF_ATOMS] =
    {
        {"H", 	1.01}, 	 {"He", 4.00}, 	 {"Li", 6.94}, 		{"Be",  9.01},
        {"B", 	10.81},  {"C", 	12.01},  {"N", 	14.01}, 	{"O",   16.00},
        {"F", 	19.00},  {"Ne", 20.18},  {"Na", 22.99}, 	{"Mg",  24.31},
        {"Al", 	26.98},  {"Si", 28.09},  {"P",  30.97}, 	{"S",   32.07},
        {"Cl", 	35.45},  {"Ar", 39.95},  {"K",  39.10}, 	{"Ca",  40.08},
        {"Sc", 	44.96},  {"Ti", 47.87},  {"V",  50.94}, 	{"Cr",  52.00},
        {"Mn", 	54.94},  {"Fe", 55.85},  {"Co", 58.93}, 	{"Ni",  58.69},
        {"Cu", 	63.55},  {"Zn", 65.38},  {"Ga", 69.72}, 	{"Ge",  72.64},
        {"As", 	74.92},  {"Se", 78.96},  {"Br", 79.90}, 	{"Kr",  83.80},
        {"Rb", 	85.47},  {"Sr", 87.62},  {"Y",  88.91}, 	{"Zr",  91.22},
        {"Nb", 	92.91},  {"Mo", 95.94},	 {"Tc", 98.91}, 	{"Ru",  101.07},
        {"Rh", 	102.91}, {"Pd", 106.40}, {"Ag", 107.87}, 	{"Cd",  112.40},
        {"In", 	114.82}, {"Sn", 118.70}, {"Sb", 121.75}, 	{"Te",  127.60},
        {"I", 	126.90}, {"Xe", 131.30}, {"Cs", 132.91}, 	{"Ba",  137.34},
        {"La", 	138.91}, {"Ce", 140.12}, {"Pr", 140.91}, 	{"Nd",  144.24},
        {"Pm", 	146.92}, {"Sm", 150.40}, {"Eu", 151.96}, 	{"Gd",  157.25},
        {"Tb", 	158.93}, {"Dy", 162.50}, {"Ho", 164.93}, 	{"Er",  167.26},
        {"Tm", 	168.93}, {"Yb", 173.04}, {"Lu", 174.97}, 	{"Hf",  178.49},
        {"Ta", 	180.95}, {"W", 	183.85}, {"Re", 186.21}, 	{"Os",  190.20},
        {"Ir", 	192.22}, {"Pt", 195.10}, {"Au", 196.97}, 	{"Hg",  200.60},
        {"Tl", 	204.37}, {"Pb", 207.20}, {"Bi", 208.98}, 	{"Po",  209.00},
        {"At", 	210.00}, {"Rn", 222.00}, {"Fr", 223.00}, 	{"Ra",  226.03},
        {"Ac", 	227.00}, {"Th", 232.04}, {"Pa", 231.04}, 	{"U",   238.03},
        {"Np", 	237.05}, {"Pu", 239.05}, {"Am", 241.06}, 	{"Cm",  247.07},
        {"Bk", 	249.08}, {"Cf", 251.08}, {"Es", 254.09}, 	{"Fm",  257.10},
        {"Md", 	258.10}, {"No", 255.00}, {"Lr", 262.10}, 	{"Rf",  261.00},
        {"Db", 	262.00}, {"Sg", 263.00}, {"Bh", 264.00}, 	{"Hs",  265.00},
        {"Mt", 	266.00}, {"Ds", 281.00}, {"Uuu", 272.00},	{"Uub", 285.00}
    };

    int i;
    const char *buffer[6] = { "(C21H25ClN2O3)6", "CO", "C2O", "C6H6", "(CH3)2CO", "CO(CH2OH)2"};
    double tests[6] = { 2333.58, 28.01, 40.02, 78.11, 58.08, 90.08};
    for(i = 0; i < 6; ++i)
        printf("%s --> %.10g\t\t(verif : %.10g)\n", buffer[i], get_molar_mass(list_of_atoms, NUMBER_OF_ATOMS, buffer[i], strlen(buffer[i])), tests[i]);
    return EXIT_SUCCESS;
}
