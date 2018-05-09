//Author: Sarah P. Flanagan
//Purpose: Simulate SNPs and mating to evaluate how many SNPs are necessary and sufficient for parentage

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random_numbers.h"

using namespace std;

class individual
{
public:
	vector<int> maternal;
	vector<int> paternal;
	int mom, dad, num_mates;

	individual()
	{
		maternal = paternal = vector<int>();
		mom = dad = num_mates = int();
	}

	void initialize(vector<vector<int>>& alleles)
	{
		int j;

		for (j = 0; j < alleles.size(); j++)
		{
			maternal.push_back(alleles[j][randnum(alleles[j].size())]);
			paternal.push_back(alleles[j][randnum(alleles[j].size())]);
		}
		num_mates = 0;
	}
};

vector<int> errors(vector<int>& alleles, int alleleA, int alleleB, double allelic_dropout, double scoring)
{
	vector<int> out_allele;
	out_allele.push_back(alleleA);
	out_allele.push_back(alleleB);
	double rand_ad, rand_sc;
	if (allelic_dropout > 0) 
	{
		rand_ad = genrand();
		if (rand_ad <= allelic_dropout)
		{
			//choose one allele to be lost
			if (genrand() < 0.5)
			{
				out_allele[1] = out_allele[0]; //they both are alleleA
			}
			else
			{
				out_allele[0] = out_allele[1]; //they both are alleleB
			}
		}
	}
	if (scoring > 0)
	{
		rand_sc = genrand();
		if (rand_sc <= scoring)
		{
			//choose one allele to be affected
			if (genrand() < 0.5)
			{
				for (int j = 0; j < alleles.size(); j++) 
				{
					if (alleles[j] != out_allele[0])//they're biallelic so this works
						out_allele[0] = alleles[j];
				}
			}
			else
			{
				for (int j = 0; j < alleles.size(); j++)
				{
					if (alleles[j] != out_allele[0])//they're biallelic so this works
						out_allele[1] = alleles[j];
				}
			}
		}
	}
	return out_allele;
}

void help_message()
{
	cout << "\n\t\tHelp Menu for ParentageSim\n";
	cout << "\nThis model creates unlinked SNPs for parents and generates a CERVUS input file for their offspring\n";
	cout << "\nFemales mate with one male but males can mate multiply\n";
	cout << "-o:\tOutput file name (simulated_genotypes.txt)\n";
	cout << "-F:\tnumber of females (50)\n";
	cout << "-M:\tnumber of males (50)\n";
	cout << "-S:\tnumber of SNPs per locus (1)\n";
	cout << "-L:\tnumber of loci (10)\n";
	cout << "-x:\tmaximum number of mates per male (5)\n";
	cout << "-f:\tfecundity, aka number of offspring produced per mating (4)\n";
	cout << "-d:\tdirectory for output -- this is the full path to include in .crv file and should correspond to the relative directory\n";
	cout << "-r:\trelative directory for output -- this is the relative path to output files to, important for windows/linux compatibility (default of ../../results/)\n";
	cout << "-m:\tNumber of microsatellites to simulate (10)\n";
	cout << "-ad:\tSNP Allelic dropout rate (0)\n";
	cout << "-se:\tSNP scoring error rate (0)\n";
	cout << "-h:\tPrint this help message.\n";
}

int main(int argc, char*argv[])
{
	int i, ii, iii, num_snps,num_loci, num_males, num_females, mate, max_num_mates, fecundity, num_offspring;
	int num_microsats, microsat_nalleles;
	double prop_moms_sampled = 0.05;
	double allelic_dropout_prob, scoring_error_prob;
	vector<individual> males, females;
	string base_name = "simulated";
	string genotype_name, offspring_name, candmoms_name, canddads_name, crv_name,dir,rel_dir;
	ofstream genotypes,offspring,cand_males,cand_females,crv;

	//default settings
	microsat_nalleles = num_microsats = 10;
	num_snps = 1; //the number of snps per locus
	num_loci = 100;
	num_males = 50;
	num_females = 50;
	max_num_mates = 5;
	fecundity = 4;
	allelic_dropout_prob = scoring_error_prob = 0;
	dir = "../../results/";
	rel_dir = "../../results/";
	//read in parameters
	string tempstring1, tempstring2;
	if (argc == 1)
	{
		help_message();
		//return 0;
	}
	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h" || tempstring1 == "--help" || tempstring1 == "--H" || tempstring1 == "--Help" || tempstring1 == "--HELP")
		{
			help_message();
			return 0;
		}
		else
		{			
			for (i = 1; i < argc; i++)
			{
				tempstring1 = argv[i];
				tempstring2 = argv[i+ 1];
				i++;
				if (tempstring1 == "-o")
					base_name = tempstring2;
				if (tempstring1 == "-M")
					num_males = atoi(tempstring2.c_str());
				if (tempstring1 == "-F")
					num_females = atoi(tempstring2.c_str());
				if (tempstring1 == "-x")
					max_num_mates= atoi(tempstring2.c_str());
				if (tempstring1 == "-f")
					fecundity = atoi(tempstring2.c_str());
				if (tempstring1 == "-S")
					num_snps = atoi(tempstring2.c_str());
				if (tempstring1 == "-L")
					num_loci = atoi(tempstring2.c_str());
				if (tempstring1 == "-d")
					dir = tempstring2;
				if (tempstring1 == "-r")
					rel_dir = tempstring2;
				if (tempstring1 == "-m")
					num_microsats = atoi(tempstring2.c_str());
				if (tempstring1 == "-ad")
					allelic_dropout_prob = atof(tempstring2.c_str());
				if (tempstring1 == "-se")
					scoring_error_prob = atof(tempstring2.c_str());
			}
		}
	}
	
	//initialize loci
	vector<vector<int>> alleles;
	for (i = 0; i < num_loci; i++)
	{
		alleles.push_back(vector<int>());
		int nsnps = randnum(num_snps) + 1;
		while (nsnps < 1) nsnps = randnum(num_snps); //sanity check
		for (ii = 0; ii < pow(2,nsnps); ii++)
		{
			alleles[i].push_back(ii+1);
		}
	}
	for (i = 0; i < num_microsats; i++)
	{
		alleles.push_back(vector<int>());
		for (ii = 0; ii < microsat_nalleles; ii++)
			alleles[num_loci + i].push_back(ii + 1);
	}
	//initialize adults
	genotype_name = rel_dir + base_name + "_genotypes.txt";
	genotypes.open(genotype_name);
	genotypes << "ID\tMom\tDad";
	for (i = 0; i < num_loci; i++)
		genotypes << "\tA" << i << "\tB" << i;
	for (i = 0; i < num_microsats; i++)
		genotypes << "\tMSATa" << i << "\tMSATb" << i;
	candmoms_name = rel_dir + base_name + "_candidate_mothers.txt";
	canddads_name = rel_dir + base_name + "_candidate_fathers.txt";
	cand_females.open(candmoms_name);
	cand_males.open(canddads_name);
	vector<int> err_alleles;
	for (i = 0; i < num_males; i++)
	{
		males.push_back(individual());
		males[i].initialize(alleles);
		cand_males << "MAL" << i << '\n';
		genotypes << '\n' << "MAL" << i << "\tNA\tNA";
		for (ii = 0; ii < num_loci; ii++)
		{
			err_alleles = errors(alleles[ii], males[i].maternal[ii], males[i].paternal[ii], allelic_dropout_prob, scoring_error_prob);
			genotypes << '\t' << err_alleles[0] << '\t' << err_alleles[1];
		}
		for (ii = num_loci; ii < num_loci + num_microsats; ii++)
		{
			genotypes << '\t' << males[i].maternal[ii] << '\t' << males[i].paternal[ii];
		}
	}
	for (i = 0; i < num_females; i++)
	{
		females.push_back(individual());
		females[i].initialize(alleles);
		cand_females << "FEM" << i << '\n';
		genotypes << '\n' << "FEM" << i << "\tNA\tNA";
		for (ii = 0; ii < num_loci; ii++) 
		{
			err_alleles = errors(alleles[ii], females[i].maternal[ii], females[i].paternal[ii], allelic_dropout_prob, scoring_error_prob);
			genotypes << '\t' << err_alleles[0] << '\t' << err_alleles[1];
		}
		for (ii = num_loci; ii < num_loci+num_microsats; ii++)
		{
			genotypes << '\t' << females[i].maternal[ii] << '\t' << females[i].paternal[ii];
		}
	}
	

	//mating
	offspring_name = rel_dir + base_name + "_offspring.txt";
	offspring.open(offspring_name);
	offspring << "OffspringID\tMomID\tDadID";	
	num_offspring = 0;
	cout << "\nCommencing mating";
	int mat_allele, pat_allele;
	for (i = 0; i < num_females; i++) //mate singly
	{
		while (females[i].num_mates == 0)
		{
			mate = randnum(num_males);
			if (males[mate].num_mates < max_num_mates)
			{
				for (iii = 0; iii < fecundity; iii++)//females mate with one male
				{									//males can mate multiply
					genotypes << '\n' << "OFFSPRING" << num_offspring << "\tFEM" << i << "\tMAL" << mate;
					offspring << '\n' << "OFFSPRING" << num_offspring << "\tFEM" << i << "\tMAL" << mate;

					for (ii = 0; ii < num_loci; ii++)
					{ 
						//from the mom
						if (genrand() < 0.5)
							mat_allele = females[i].maternal[ii];
						else
							mat_allele = females[i].paternal[ii];
						//from the dad
						if (genrand() < 0.5)
							pat_allele = males[i].maternal[ii];
						else
							pat_allele = males[i].paternal[ii];
						err_alleles = errors(alleles[ii], mat_allele, pat_allele, allelic_dropout_prob, scoring_error_prob);
						genotypes << '\t' << err_alleles[0] << '\t' << err_alleles[1];
					}
					for (ii = num_loci; ii < num_loci+num_microsats; ii++)
					{
						//from the mom
						if (genrand() < 0.5)
							genotypes << '\t' << females[i].maternal[ii];
						else
							genotypes << '\t' << females[i].paternal[ii];
						//from the dad
						if (genrand() < 0.5)
							genotypes << '\t' << males[i].maternal[ii];
						else
							genotypes << '\t' << males[i].paternal[ii];
					}
					num_offspring++;
					if (num_offspring % 100 == 0)
						cout << '\n' << num_offspring << " offspring have been created" << std::flush;
				}
				males[mate].num_mates++;
				females[i].num_mates++;
			}
		}
		
	}
	genotypes.close();
	offspring.close();
	cand_females.close();
	cand_males.close();
	cout << "\nSimulation completed.\n" << std::flush;


	return 0;
}