/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./cancer_immune_3D.h"
#include <stdlib.h>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip> 
#include <cmath>

Cell_Definition* pImmuneCell_EGFR; 
Cell_Definition* pImmuneCell_IL13;
Cell_Definition* pImmuneCell_HER2;

void create_immune_cell_type_EGFR( void )
{
	pImmuneCell_EGFR = find_cell_definition( "immune cell EGFR" ); 
	
	static int oxygen_ID = microenvironment.find_density_index( "oxygen" ); 
	static int immuno_ID = microenvironment.find_density_index( "immunostimulatory factor" ); 
	
	// reduce o2 uptake 
	
	pImmuneCell_EGFR->phenotype.secretion.uptake_rates[oxygen_ID] *= 
		parameters.doubles("immune_o2_relative_uptake");  
	
	pImmuneCell_EGFR->phenotype.mechanics.cell_cell_adhesion_strength *= 
		parameters.doubles("immune_relative_adhesion"); 
	pImmuneCell_EGFR->phenotype.mechanics.cell_cell_repulsion_strength *= 
		parameters.doubles("immune_relative_repulsion"); 
		
	// figure out mechanics parameters 
	
	pImmuneCell_EGFR->phenotype.mechanics.relative_maximum_attachment_distance 
		= pImmuneCell_EGFR->custom_data["max_attachment_distance"] / pImmuneCell_EGFR->phenotype.geometry.radius ; 
		
	pImmuneCell_EGFR->phenotype.mechanics.attachment_elastic_constant 
		= pImmuneCell_EGFR->custom_data["elastic_coefficient"]; 		
	
	pImmuneCell_EGFR->phenotype.mechanics.relative_detachment_distance 
		= pImmuneCell_EGFR->custom_data["max_attachment_distance" ] / pImmuneCell_EGFR->phenotype.geometry.radius ; 		
	
	// set functions 
	
	pImmuneCell_EGFR->functions.update_phenotype = NULL; 
	pImmuneCell_EGFR->functions.custom_cell_rule = immune_cell_rule; 
	pImmuneCell_EGFR->functions.update_migration_bias = immune_cell_motility;
	pImmuneCell_EGFR->functions.contact_function = adhesion_contact_function; 
	
	// set custom data values 
	
	return; 
}

void create_immune_cell_type_IL13(void)
{
	pImmuneCell_IL13 = find_cell_definition("immune cell IL13");

	static int oxygen_ID = microenvironment.find_density_index("oxygen");
	static int immuno_ID = microenvironment.find_density_index("immunostimulatory factor");

	// reduce o2 uptake 

	pImmuneCell_IL13->phenotype.secretion.uptake_rates[oxygen_ID] *=
		parameters.doubles("immune_o2_relative_uptake");

	pImmuneCell_IL13->phenotype.mechanics.cell_cell_adhesion_strength *=
		parameters.doubles("immune_relative_adhesion");
	pImmuneCell_IL13->phenotype.mechanics.cell_cell_repulsion_strength *=
		parameters.doubles("immune_relative_repulsion");

	// figure out mechanics parameters 

	pImmuneCell_IL13->phenotype.mechanics.relative_maximum_attachment_distance
		= pImmuneCell_IL13->custom_data["max_attachment_distance"] / pImmuneCell_IL13->phenotype.geometry.radius;

	pImmuneCell_IL13->phenotype.mechanics.attachment_elastic_constant
		= pImmuneCell_IL13->custom_data["elastic_coefficient"];

	pImmuneCell_IL13->phenotype.mechanics.relative_detachment_distance
		= pImmuneCell_IL13->custom_data["max_attachment_distance"] / pImmuneCell_IL13->phenotype.geometry.radius;

	// set functions 

	pImmuneCell_IL13->functions.update_phenotype = NULL;
	pImmuneCell_IL13->functions.custom_cell_rule = immune_cell_rule;
	pImmuneCell_IL13->functions.update_migration_bias = immune_cell_motility;
	pImmuneCell_IL13->functions.contact_function = adhesion_contact_function;

	// set custom data values 

	return;
}

void create_immune_cell_type_HER2(void)
{
	pImmuneCell_HER2 = find_cell_definition("immune cell HER2");

	static int oxygen_ID = microenvironment.find_density_index("oxygen");
	static int immuno_ID = microenvironment.find_density_index("immunostimulatory factor");

	// reduce o2 uptake 

	pImmuneCell_HER2->phenotype.secretion.uptake_rates[oxygen_ID] *=
		parameters.doubles("immune_o2_relative_uptake");

	pImmuneCell_HER2->phenotype.mechanics.cell_cell_adhesion_strength *=
		parameters.doubles("immune_relative_adhesion");
	pImmuneCell_HER2->phenotype.mechanics.cell_cell_repulsion_strength *=
		parameters.doubles("immune_relative_repulsion");

	// figure out mechanics parameters 

	pImmuneCell_HER2->phenotype.mechanics.relative_maximum_attachment_distance
		= pImmuneCell_HER2->custom_data["max_attachment_distance"] / pImmuneCell_HER2->phenotype.geometry.radius;

	pImmuneCell_HER2->phenotype.mechanics.attachment_elastic_constant
		= pImmuneCell_HER2->custom_data["elastic_coefficient"];

	pImmuneCell_HER2->phenotype.mechanics.relative_detachment_distance
		= pImmuneCell_HER2->custom_data["max_attachment_distance"] / pImmuneCell_HER2->phenotype.geometry.radius;

	// set functions 

	pImmuneCell_HER2->functions.update_phenotype = NULL;
	pImmuneCell_HER2->functions.custom_cell_rule = immune_cell_rule;
	pImmuneCell_HER2->functions.update_migration_bias = immune_cell_motility;
	pImmuneCell_HER2->functions.contact_function = adhesion_contact_function;

	// set custom data values 

	return;
}
void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	SeedRandom( parameters.ints("random_seed") ); 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	
	
	cell_defaults.parameters.o2_proliferation_saturation = 38.0;  
	cell_defaults.parameters.o2_reference = 38.0; 

	static int oxygen_ID = microenvironment.find_density_index( "oxygen" ); // 0 
	static int immuno_ID = microenvironment.find_density_index( "immunostimulatory factor" ); // 1
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/

	initialize_cell_definitions_from_pugixml(); 
	
	// change the max cell-cell adhesion distance 
	cell_defaults.phenotype.mechanics.relative_maximum_attachment_distance = 
		cell_defaults.custom_data["max_attachment_distance"] / cell_defaults.phenotype.geometry.radius;
		
	cell_defaults.phenotype.mechanics.relative_detachment_distance 
		= cell_defaults.custom_data["max_attachment_distance"] / cell_defaults.phenotype.geometry.radius ; 
		
	cell_defaults.phenotype.mechanics.attachment_elastic_constant 
		= cell_defaults.custom_data[ "elastic_coefficient" ];	
		
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_and_immune_stimulation; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = adhesion_contact_function; 
	cell_defaults.functions.update_migration_bias = NULL; 

	// create the immune cell type 
	create_immune_cell_type_EGFR(); 
	create_immune_cell_type_IL13();
	create_immune_cell_type_HER2();

	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		std::cout << "Warning: overriding 2D setting to return to 3D" << std::endl; 
		default_microenvironment_options.simulate_2D = false; 
	}
	
	initialize_microenvironment(); 	

	return; 
}	

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc = 0, yc = 0, zc = 0;
	double x_spacing = cell_radius * sqrt(3);
	double y_spacing = cell_radius * 2;
	double z_spacing = cell_radius * sqrt(3);


	std::vector<double> tempPoint(3, 0.0);
	// std::vector<double> cylinder_center(3,0.0);

	for (double z = -sphere_radius; z < sphere_radius; z += z_spacing, zc++)
	{
		for (double x = -sphere_radius; x < sphere_radius; x += x_spacing, xc++)
		{
			for (double y = -sphere_radius; y < sphere_radius; y += y_spacing, yc++)
			{
				tempPoint[0] = x + (zc % 2) * 0.5 * cell_radius;
				tempPoint[1] = y + (xc % 2) * cell_radius;
				tempPoint[2] = z;

				if (sqrt(norm_squared(tempPoint)) < sphere_radius)
				{
					cells.push_back(tempPoint);
				}
			}

		}
	}
	return cells;

}

void introduce_il13_immune_cells( std::vector<double> il13_center )
{
	double tumor_radius = -9e9; // 250.0; 
	double temp_radius = 0.0; 
	
	// for the loop, deal with the (faster) norm squared 
	for( int i=0; i < (*all_cells).size() ; i++ )
	{
		temp_radius = norm_squared( (*all_cells)[i]->position ); 
		if( temp_radius > tumor_radius )
		{ tumor_radius = temp_radius; }
	}
	// now square root to get to radius 
	tumor_radius = sqrt( tumor_radius ); 
	
	// if this goes wackadoodle, choose 250 
	if( tumor_radius < 50.0 )
	{ tumor_radius = 50.0; }
	
//	std::cout << "current tumor radius: " << tumor_radius << std::endl; 
	
	// now seed immune cells 
	
	int number_of_immune_cells = 
		parameters.ints("number_of_il13_immune_cells"); // 7500; // 100; // 40; 
	double radius_inner = tumor_radius + 
		parameters.doubles("initial_min_immune_distance_from_tumor"); 30.0; // 75 // 50; 
	double radius_outer = radius_inner + 
		parameters.doubles("thickness_of_immune_seeding_region"); // 75.0; // 100; // 1000 - 50.0; 
	
	double mean_radius = 0.5*(radius_inner + radius_outer); 
	double std_radius = 0.33*( radius_outer-radius_inner)/2.0;

	//Cell* pCell = create_cell(*pImmuneCell_EGFR)
	
	//there are about 14 cells per unit of radius, so given the number of cells we want to introduce, we divide it by 14 to have an approximated radius of sphere
	int sphere_radius = 4.5;
	std::vector<std::vector<double>> sphere_positions = create_cell_sphere_positions(4.5, sphere_radius); //cell radius 8, sphere radius 100
	while (sphere_positions.size() < number_of_immune_cells)
	{
		sphere_radius += 1;
		sphere_positions = create_cell_sphere_positions(4.5,sphere_radius);
	}
	//std::cout << sphere_positions.size() << " IL13 immune cells are introduced" << std::endl;
	//int number_of_rows = sizeof sphere_positions / sizeof sphere_positions[0];
	//int number_of_columns = sizeof sphere_positions[0] / sizeof sphere_positions[0][0];
	//std::cout << "rows of sphere positions: " << number_of_rows << std::endl;
	//std::cout << "columns of sphere positions: " << number_of_columns<< std::endl;
	//std::cout << "x=" << sphere_positions[0][0] << ", y=" << sphere_positions[0][1] << ", z=" << sphere_positions[0][2] << std::endl;
  //std::cout << "sphere radius= "<<sphere_radius<<std::endl;
	for (int i = 0; i < number_of_immune_cells; i++)
	{
		Cell* pCell = create_cell(*pImmuneCell_IL13);
		
		double x = sphere_positions[i][0] - 250; //-250 for PBT018 treating IL13 region;+100 for PBT025; -200 for PBT030
		double y = sphere_positions[i][1] + 150; //+150 for PBT018 treating IL13 region;+100 for PBT025; +100 for PBT030
		double z = sphere_positions[i][2];
		pCell -> assign_position(x, y, z);
	}
	
	


	/*for (int i = 0; i < number_of_immune_cells; i++)
	{
		//double theta = UniformRandom() * 6.283185307179586476925286766559;
		//double phi = acos(2.0 * UniformRandom() - 1.0);

		//double radius = NormalRandom(mean_radius, std_radius);

		//Cell* pCell = create_cell(*pImmuneCell);
		
		//pCell->assign_position(radius * cos(theta) * sin(phi), radius * sin(theta) * sin(phi), radius * cos(phi));
		if (i % 3 == 0) {
			double x = -700;
			double y = 1500 * UniformRandom() - 750;
			double z = 0;
			Cell* pCell = create_cell(*pImmuneCell_EGFR);
		}
		else if (i % 3 == 1) {
			double x = -700;
			double y = 1500 * UniformRandom() - 750;
			double z = 0;
			Cell* pCell = create_cell(*pImmuneCell_HER2);

		}
		else {
			double x = -700;
			double y = 1500 * UniformRandom() - 750;
			double z = 0;
			Cell* pCell = create_cell(*pImmuneCell_IL13);
		}


		
	}*/

	
	// std::vector<double> cylinder_center(3,0.0);

	
	
	return; 
}

void introduce_il13_immune_cells_1( void )
{
	double tumor_radius = -9e9; // 250.0; 
	double temp_radius = 0.0; 
	
	// for the loop, deal with the (faster) norm squared 
	for( int i=0; i < (*all_cells).size() ; i++ )
	{
		temp_radius = norm_squared( (*all_cells)[i]->position ); 
		if( temp_radius > tumor_radius )
		{ tumor_radius = temp_radius; }
	}
	// now square root to get to radius 
	tumor_radius = sqrt( tumor_radius ); 
	
	// if this goes wackadoodle, choose 250 
	if( tumor_radius < 50.0 )
	{ tumor_radius = 50.0; }
	
//	std::cout << "current tumor radius: " << tumor_radius << std::endl; 
	
	// now seed immune cells 
	
	int number_of_immune_cells = 
		parameters.ints("number_of_il13_immune_cells_1"); // 7500; // 100; // 40; 
	double radius_inner = tumor_radius + 
		parameters.doubles("initial_min_immune_distance_from_tumor"); 30.0; // 75 // 50; 
	double radius_outer = radius_inner + 
		parameters.doubles("thickness_of_immune_seeding_region"); // 75.0; // 100; // 1000 - 50.0; 
	
	double mean_radius = 0.5*(radius_inner + radius_outer); 
	double std_radius = 0.33*( radius_outer-radius_inner)/2.0;

	//Cell* pCell = create_cell(*pImmuneCell_EGFR)
	
	//there are about 14 cells per unit of radius, so given the number of cells we want to introduce, we divide it by 14 to have an approximated radius of sphere
	int sphere_radius = 4.5;
	std::vector<std::vector<double>> sphere_positions = create_cell_sphere_positions(4.5, sphere_radius); //cell radius 8, sphere radius 100
	while (sphere_positions.size() < number_of_immune_cells)
	{
		sphere_radius += 1;
		sphere_positions = create_cell_sphere_positions(4.5,sphere_radius);
	}
	//std::cout << sphere_positions.size() << " IL13 immune cells are introduced" << std::endl;
	//int number_of_rows = sizeof sphere_positions / sizeof sphere_positions[0];
	//int number_of_columns = sizeof sphere_positions[0] / sizeof sphere_positions[0][0];
	//std::cout << "rows of sphere positions: " << number_of_rows << std::endl;
	//std::cout << "columns of sphere positions: " << number_of_columns<< std::endl;
	//std::cout << "x=" << sphere_positions[0][0] << ", y=" << sphere_positions[0][1] << ", z=" << sphere_positions[0][2] << std::endl;
  //std::cout << "sphere radius= "<<sphere_radius<<std::endl;
	for (int i = 0; i < number_of_immune_cells; i++)
	{
		Cell* pCell = create_cell(*pImmuneCell_IL13);
		double x = sphere_positions[i][0]; //+100 for PBT018 treating IL13 region;+100 for PBT025; -200 for PBT030
		double y = sphere_positions[i][1]; //+200 for PBT018 treating IL13 region;+100 for PBT025; +100 for PBT030
		double z = sphere_positions[i][2];
		pCell -> assign_position(x, y, z);
	}
	
	


	/*for (int i = 0; i < number_of_immune_cells; i++)
	{
		//double theta = UniformRandom() * 6.283185307179586476925286766559;
		//double phi = acos(2.0 * UniformRandom() - 1.0);

		//double radius = NormalRandom(mean_radius, std_radius);

		//Cell* pCell = create_cell(*pImmuneCell);
		
		//pCell->assign_position(radius * cos(theta) * sin(phi), radius * sin(theta) * sin(phi), radius * cos(phi));
		if (i % 3 == 0) {
			double x = -700;
			double y = 1500 * UniformRandom() - 750;
			double z = 0;
			Cell* pCell = create_cell(*pImmuneCell_EGFR);
		}
		else if (i % 3 == 1) {
			double x = -700;
			double y = 1500 * UniformRandom() - 750;
			double z = 0;
			Cell* pCell = create_cell(*pImmuneCell_HER2);

		}
		else {
			double x = -700;
			double y = 1500 * UniformRandom() - 750;
			double z = 0;
			Cell* pCell = create_cell(*pImmuneCell_IL13);
		}


		
	}*/

	
	// std::vector<double> cylinder_center(3,0.0);

	
	
	return; 
}



void introduce_her2_immune_cells(std::vector<double> her2_center)
{
	double tumor_radius = -9e9; // 250.0; 
	double temp_radius = 0.0;

	// for the loop, deal with the (faster) norm squared 
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		temp_radius = norm_squared((*all_cells)[i]->position);
		if (temp_radius > tumor_radius)
		{
			tumor_radius = temp_radius;
		}
	}
	// now square root to get to radius 
	tumor_radius = sqrt(tumor_radius);

	// if this goes wackadoodle, choose 250 
	if (tumor_radius < 50.0)
	{
		tumor_radius = 50.0;
	}

	//std::cout << "current tumor radius: " << tumor_radius << std::endl;

	// now seed immune cells 

	int number_of_immune_cells =
		parameters.ints("number_of_her2_immune_cells"); // 7500; // 100; // 40; 
	double radius_inner = tumor_radius +
		parameters.doubles("initial_min_immune_distance_from_tumor"); 30.0; // 75 // 50; 
	double radius_outer = radius_inner +
		parameters.doubles("thickness_of_immune_seeding_region"); // 75.0; // 100; // 1000 - 50.0; 

	double mean_radius = 0.5 * (radius_inner + radius_outer);
	double std_radius = 0.33 * (radius_outer - radius_inner) / 2.0;
	//Cell* pCell = create_cell(*pImmuneCell_EGFR)
	int sphere_radius = 4.5;
	std::vector<std::vector<double>> sphere_positions = create_cell_sphere_positions(4.5, sphere_radius); //cell radius 8, sphere radius 8
	while (sphere_positions.size() < number_of_immune_cells)
	{
		sphere_radius += 1;
		sphere_positions = create_cell_sphere_positions(4.5, sphere_radius);
	}
	//std::cout << sphere_radius << " HER2 immune cells are introduced" << std::endl;
	//std::cout <<"Number of spaces for HER2: " << sphere_positions.size() << std::endl;
		
	//int number_of_rows = sizeof sphere_positions / sizeof sphere_positions[0];
	//int number_of_columns = sizeof sphere_positions[0] / sizeof sphere_positions[0][0];
	//std::cout << "rows of sphere positions: " << number_of_rows << std::endl;
	//std::cout << "columns of sphere positions: " << number_of_columns<< std::endl;
	//std::cout << "x=" << sphere_positions[0][0] << ", y=" << sphere_positions[0][1] << ", z=" << sphere_positions[0][2] << std::endl;
	
	for (int i = 0; i < number_of_immune_cells; i++)
	{
		Cell* pCell = create_cell(*pImmuneCell_HER2);
		double x = sphere_positions[i][0] - 50; //-250 for PBT018 treating IL13 region; -400 for PBT025; -50 for treating HER2 region in PBT018
		double y = sphere_positions[i][1] - 50; //+150 for PBT018 treating IL13 region; -300 for PBT025; -50 for treating HER2 region in PBT018
		double z = sphere_positions[i][2];
		pCell->assign_position(x, y, z);
		//std::cout << "i=" << i << std::endl;
	}



	/*for (int i = 0; i < number_of_immune_cells; i++)
	{
		//double theta = UniformRandom() * 6.283185307179586476925286766559;
		//double phi = acos(2.0 * UniformRandom() - 1.0);

		//double radius = NormalRandom(mean_radius, std_radius);

		//Cell* pCell = create_cell(*pImmuneCell);

		//pCell->assign_position(radius * cos(theta) * sin(phi), radius * sin(theta) * sin(phi), radius * cos(phi));
		if (i % 3 == 0) {
			double x = -700;
			double y = 1500 * UniformRandom() - 750;
			double z = 0;
			Cell* pCell = create_cell(*pImmuneCell_EGFR);
		}
		else if (i % 3 == 1) {
			double x = -700;
			double y = 1500 * UniformRandom() - 750;
			double z = 0;
			Cell* pCell = create_cell(*pImmuneCell_HER2);

		}
		else {
			double x = -700;
			double y = 1500 * UniformRandom() - 750;
			double z = 0;
			Cell* pCell = create_cell(*pImmuneCell_IL13);
		}



	}*/


	// std::vector<double> cylinder_center(3,0.0);



	return;
}

void introduce_egfr_immune_cells(std::vector<double> egfr_center)
{
	double tumor_radius = -9e9; // 250.0; 
	double temp_radius = 0.0;

	// for the loop, deal with the (faster) norm squared 
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		temp_radius = norm_squared((*all_cells)[i]->position);
		if (temp_radius > tumor_radius)
		{
			tumor_radius = temp_radius;
		}
	}
	// now square root to get to radius 
	tumor_radius = sqrt(tumor_radius);

	// if this goes wackadoodle, choose 250 
	if (tumor_radius < 50.0)
	{
		tumor_radius = 50.0;
	}

	//std::cout << "current tumor radius: " << tumor_radius << std::endl;

	// now seed immune cells 

	int number_of_immune_cells =
		parameters.ints("number_of_egfr_immune_cells"); // 7500; // 100; // 40; 
	double radius_inner = tumor_radius +
		parameters.doubles("initial_min_immune_distance_from_tumor"); 30.0; // 75 // 50; 
	double radius_outer = radius_inner +
		parameters.doubles("thickness_of_immune_seeding_region"); // 75.0; // 100; // 1000 - 50.0; 

	double mean_radius = 0.5 * (radius_inner + radius_outer);
	double std_radius = 0.33 * (radius_outer - radius_inner) / 2.0;
	//Cell* pCell = create_cell(*pImmuneCell_EGFR)
	int sphere_radius = 4.5;
	std::vector<std::vector<double>> sphere_positions = create_cell_sphere_positions(4.5, sphere_radius); //cell radius 8, sphere radius 100
	while (sphere_positions.size() < number_of_immune_cells)
	{
		sphere_radius += 1;
		sphere_positions = create_cell_sphere_positions(4.5, sphere_radius);
	}
	//std::cout << sphere_radius << " HER2 immune cells are introduced" << std::endl;
	//std::cout <<"Number of spaces for HER2: " << sphere_positions.size() << std::endl;

	//int number_of_rows = sizeof sphere_positions / sizeof sphere_positions[0];
	//int number_of_columns = sizeof sphere_positions[0] / sizeof sphere_positions[0][0];
	//std::cout << "rows of sphere positions: " << number_of_rows << std::endl;
	//std::cout << "columns of sphere positions: " << number_of_columns<< std::endl;
	//std::cout << "x=" << sphere_positions[0][0] << ", y=" << sphere_positions[0][1] << ", z=" << sphere_positions[0][2] << std::endl;

	for (int i = 0; i < number_of_immune_cells; i++)
	{
		Cell* pCell = create_cell(*pImmuneCell_EGFR);
		double x = sphere_positions[i][0];
		double y = sphere_positions[i][1];
		double z = sphere_positions[i][2];
		pCell->assign_position(x, y, z);
		//std::cout << "i=" << i << std::endl;
	}


/*
	for (int i = 0; i < number_of_immune_cells; i++)
	{
		//double theta = UniformRandom() * 6.283185307179586476925286766559;
		//double phi = acos(2.0 * UniformRandom() - 1.0);

		//double radius = NormalRandom(mean_radius, std_radius);

		//Cell* pCell = create_cell(*pImmuneCell);

		//pCell->assign_position(radius * cos(theta) * sin(phi), radius * sin(theta) * sin(phi), radius * cos(phi));
		if (i % 3 == 0) {
			double x = -700;
			double y = 1500 * UniformRandom() - 750;
			double z = 0;
			Cell* pCell = create_cell(*pImmuneCell_EGFR);
		}
		else if (i % 3 == 1) {
			double x = -700;
			double y = 1500 * UniformRandom() - 750;
			double z = 0;
			Cell* pCell = create_cell(*pImmuneCell_HER2);

		}
		else {
			double x = -700;
			double y = 1500 * UniformRandom() - 750;
			double z = 0;
			Cell* pCell = create_cell(*pImmuneCell_IL13);
		}



	}
*/

	// std::vector<double> cylinder_center(3,0.0);



	return;
}

//write my own function to create cells positions 
std::vector<std::vector<double>> create_cell_position( void )
{
	std::fstream coordinates;
	std::vector<std::vector<double>> cells;
	//fstream antigens;
	coordinates.open("/bigdata/cholab/rli122/Physicell/PhysiCell/Glioma_Tissue_Michael_Barish_2022-10-30/PBT018/coor_anti.txt", std::ios::in);
	//antigens.open("D:\\MyCellClass\\MyCellClass\\antigens.txt", ios::in);

	if (coordinates.is_open()) {
		std::string line;
		while (getline(coordinates, line)) {

			std::string temp;

			std::vector<double> vec ;


			size_t i = 0, start = 0, end;

			do {
				end = line.find_first_of(' ', start);
				temp = line.substr(start, end);
				if (isdigit(temp[0]))
				{
					vec.push_back(atof(temp.c_str()));
					++i;
					
				}
				start = end + 1;
			} while (start);

			
			//vec[0] = vec[0] - 44.4244;
			//vec[1] = vec[1] - 80.6364;
			//vec[2] = vec[2] - 15.5;
			//std::cout << "x=" << vec[0] << ", y=" << vec[1] << ", z=" << vec[2] << std::endl;
			
			cells.push_back(vec);
			//std::cout << "x=" << cells[1][0] << ", y=" << cells[1][1] << ", z=" << cells[1][2] << std::endl;


		}

		coordinates.close();

	}
	return cells;
}

//write a function to place adjacent cells 


std::vector<std::vector<double>> setup_tissue( void )
{
	// place a cluster of tumor cells at the center 
	
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	//double cell_spacing = 200 * cell_radius;
	double tumor_radius = 
		parameters.doubles("tumor_radius"); // 250.0; 
	
	Cell* pCell = NULL; 
	double ratio = 1.7;

	
  int number_of_cancer_cells = 12500;
  int sphere_radius = 8;
  
	
  std::vector<std::vector<double>> sphere_positions = create_cell_sphere_positions(8, sphere_radius); //cell radius 8, sphere radius 100
	while (sphere_positions.size() < number_of_cancer_cells)
	{
		sphere_radius += 1;
		sphere_positions = create_cell_sphere_positions(8, sphere_radius);
	}
  
 
	std::vector<std::vector<double>> positions = create_cell_position();
	std::cout << "creating " << positions.size()-1 << " closely-packed tumor cells ... " << std::endl; 
 
	//std::cout << "showing cells: " << positions[0,0] << std::endl;

	static double imm_mean = parameters.doubles("tumor_mean_immunogenicity"); 
	static double imm_sd = parameters.doubles("tumor_immunogenicity_standard_deviation"); 
	cell_defaults.phenotype.volume.multiply_by_ratio(ratio);
  double il13_sum = 0.0;
  double her2_sum = 0.0;
  double egfr_sum = 0.0;
  double il13_x_sum = 0.0;
  double il13_y_sum = 0.0;
  double il13_z_sum = 0.0;
  double her2_x_sum = 0.0;
  double her2_y_sum = 0.0;
  double her2_z_sum = 0.0;
  double egfr_x_sum = 0.0;
  double egfr_y_sum = 0.0;
  double egfr_z_sum = 0.0;
  
		
	/* use the following loop for the given initial condition*/
 
  for( int i=0; i < positions.size()-1; i++ )
	{
		pCell = create_cell(); // tumor cell 
	
		int index_egfr = pCell->custom_data.find_variable_index("EGFR_antigen_level");
		int index_il13 = pCell->custom_data.find_variable_index("IL13_antigen_level");
		int index_her2 = pCell->custom_data.find_variable_index("HER2_antigen_level");
    double egfr_antigen_threshold = pCell->custom_data["egfr_antigen_threshold"];
    double her2_antigen_threshold = pCell->custom_data["her2_antigen_threshold"];
    double il13_antigen_threshold = pCell->custom_data["il13_antigen_threshold"];
	
		double x = 1.5 * cell_defaults.phenotype.geometry.radius * (positions[i][0]- 50.5);
	
		double y = 1.5 * cell_defaults.phenotype.geometry.radius * (positions[i][1]- 40.5);
	
		double z = 1.5 * cell_defaults.phenotype.geometry.radius * (positions[i][2]-5.5);

		pCell->assign_position( x, y, z);
    //std::cout<<"egfr_antigen_threshold="<<egfr_antigen_threshold<<std::endl;
    
		pCell->custom_data["oncoprotein"] = NormalRandom( imm_mean, imm_sd );
		/*
		if (positions[i][3] >= 30 || positions[i][4]>=30 || positions[i][5]>=30) {
			pCell -> custom_data[index_il13] = 65;
		
		pCell->custom_data[index_her2] = 0;
		pCell->custom_data[index_egfr] = 0;}*/
   
   if (positions[i][3]>=egfr_antigen_threshold){
   pCell->custom_data[index_egfr]=positions[i][3];
   egfr_x_sum += positions[i][3]*x;
   egfr_y_sum += positions[i][3]*y;
   egfr_z_sum += positions[i][3]*z;
   egfr_sum += positions[i][3];
   }
   if (positions[i][4]>=il13_antigen_threshold){
   pCell->custom_data[index_il13]=positions[i][4];
    il13_x_sum += positions[i][4]*x;
   il13_y_sum += positions[i][4]*y; 
   il13_z_sum += positions[i][4]*z;
   il13_sum += positions[i][4];
   }
   if (positions[i][5]>=her2_antigen_threshold){
   pCell->custom_data[index_her2]=positions[i][5];
    her2_x_sum += positions[i][5]*x;
   her2_y_sum += positions[i][5]*y;
   her2_z_sum += positions[i][5]*z;
   her2_sum += positions[i][5];
   }
   
	 
		if( pCell->custom_data["oncoprotein"] < 0.0 )
		{ pCell->custom_data["oncoprotein"] = 0.0; } 
   
	}
 
  double x_avg_egfr = round(egfr_x_sum/egfr_sum);
  double y_avg_egfr = round(egfr_y_sum/egfr_sum);
  double z_avg_egfr = round(egfr_z_sum/egfr_sum);
  double x_avg_il13 = round(il13_x_sum/il13_sum);
  double y_avg_il13 = round(il13_y_sum/il13_sum);
  double z_avg_il13 = round(il13_z_sum/il13_sum);
  double x_avg_her2 = round(her2_x_sum/her2_sum);
  double y_avg_her2 = round(her2_y_sum/her2_sum);
  double z_avg_her2 = round(her2_z_sum/her2_sum);
  std::vector<double> egfr_center = { x_avg_egfr,  y_avg_egfr,  z_avg_egfr};
  std::vector<double> il13_center = { x_avg_il13,  y_avg_il13,  z_avg_il13};
  std::vector<double> her2_center = { x_avg_her2,  y_avg_her2,  z_avg_her2};
  std::vector<std::vector<double>> three_centers={{egfr_center},{il13_center},{her2_center}};


 /* use the following loop for parameter calibration*/
 
 /*
  for( int i=0; i < number_of_cancer_cells; i++ )
	 {
	 	 pCell = create_cell(); // tumor cell 
		 pCell->assign_position( sphere_positions[i] );
     int index_egfr = pCell->custom_data.find_variable_index("EGFR_antigen_level");
		 int index_il13 = pCell->custom_data.find_variable_index("IL13_antigen_level");
		 int index_her2 = pCell->custom_data.find_variable_index("HER2_antigen_level");
	
		 //std::cout << "positions: " << positions[i] << std::endl;
		 pCell->custom_data["oncoprotein"] = NormalRandom( imm_mean, imm_sd );
     pCell -> custom_data[index_il13] = 35;
     pCell->custom_data[index_her2] = 0;
		 pCell->custom_data[index_egfr] = 0;
     
		 if( pCell->custom_data["oncoprotein"] < 0.0 )
		 { pCell->custom_data["oncoprotein"] = 0.0; } 
	 }
 */
	
	double sum = 0.0; 
	double min = 9e9; 
	double max = -9e9; 
	for( int i=0; i < all_cells->size() ; i++ )
	{
		double r = (*all_cells)[i]->custom_data["oncoprotein"]; 
		sum += r;
		if( r < min )
		{ min = r; } 
		if( r > max )
		{ max = r; }
		//std::cout << sum << std::endl;
	}
	double mean = sum / ( all_cells->size() + 1e-15 ); 
	// compute standard deviation 
	sum = 0.0; 
	for( int i=0; i < all_cells->size(); i++ )
	{
		sum +=  ( (*all_cells)[i]->custom_data["oncoprotein"] - mean )*( (*all_cells)[i]->custom_data["oncoprotein"] - mean );
		//std::cout << sum << std:endl;
	}
	double standard_deviation = sqrt( sum / ( all_cells->size() - 1.0 + 1e-15 ) ); 
	
	std::cout << std::endl << "Oncoprotein summary: " << std::endl
			  << "===================" << std::endl; 
	std::cout << "mean: " << mean << std::endl; 
	std::cout << "standard deviation: " << standard_deviation << std::endl; 
	std::cout << "[min max]: [" << min << " " << max << "]" << std::endl << std::endl; 
 
  
	return three_centers; 
}

// custom cell phenotype function to scale immunostimulatory factor with hypoxia 
void tumor_cell_phenotype_with_and_immune_stimulation( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" ); 
	
	// update secretion rates based on hypoxia 
	
	static int o2_index = microenvironment.find_density_index( "oxygen" ); 
	static int immune_factor_index = microenvironment.find_density_index( "immunostimulatory factor" ); 
	double o2 = pCell->nearest_density_vector()[o2_index];	

	phenotype.secretion.secretion_rates[immune_factor_index] = 10.0; 
	
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	
	// if cell is dead, don't bother with future phenotype changes. 
	// set it to secrete the immunostimulatory factor 
	if( phenotype.death.dead == true )
	{
		phenotype.secretion.secretion_rates[immune_factor_index] = 10; 
		pCell->functions.update_phenotype = NULL; 		
		return; 
	}

	// multiply proliferation rate by the oncoprotein 
	phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ; 
	
	return; 
}

std::vector<std::string> cancer_immune_coloring_function( Cell* pCell )
{
	static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" );
	static int index_egfr = pCell->custom_data.find_variable_index("EGFR_antigen_level");
	static int index_il13 = pCell->custom_data.find_variable_index("IL13_antigen_level");
	static int index_her2 = pCell->custom_data.find_variable_index("HER2_antigen_level");
 double egfr_antigen_threshold = pCell->custom_data["egfr_antigen_threshold"];
    double her2_antigen_threshold = pCell->custom_data["her2_antigen_threshold"];
    double il13_antigen_threshold = pCell->custom_data["il13_antigen_threshold"];
	//std::cout << "HER2: " << pCell->custom_data[index_her2] << std::endl;
	//std::cout << "IL13: "<< pCell->custom_data[index_il13] << std::endl;
	//std::cout << "EGFR: " << pCell-> custom_data[index_egfr] << std::endl;
	
	// immune are black
	std::vector< std::string > output( 4, "black" ); 
	
	if( pCell->type == 1 )
	{ 
		//type 1 gives EGFR
		//output[0] = "lime";
		//output[1] = "lime";
		//output[2] = "green"; 
		output[0] = "rgb(189, 183, 107)"; //dark khaki
		output[1] = "rgb(189, 183, 107)";
		output[2] = "rgb(255, 218, 185)"; //peachpuff
		return output;
	} 
	if (pCell->type == 2)
	{
		//type 2 gives IL13
		output[0] = "rgb(255, 255, 0)"; //yellow
		output[1] = "rgb(255, 255, 0)";
		output[2] = "rgb(255, 255, 224)"; //lightyellow
		return output;
	}
	if (pCell->type == 3)
	{
		//type 3 gives HER2
		output[0] = "rgb(139, 0, 139)"; //darkmagenta
		output[1] = "rgb(139, 0, 139)"; 
		output[2] = "rgb(255, 0, 255)"; //magenta
		return output;
	}

	// if I'm under attack, color me 
	if( pCell->state.attached_cells.size() > 0 )
	{
		output[0] = "darkcyan"; // orangered // "purple"; // 128,0,128
		output[1] = "black"; // "magenta"; 
		output[2] = "cyan"; // "magenta"; //255,0,255
		return output; 
	}
	
	// live cells are green, but shaded by oncoprotein value 
	/*if (pCell->phenotype.death.dead == false)
	{
		int oncoprotein = (int) round( 0.5 * pCell->custom_data[oncoprotein_i] * 255.0 ); 
		char szTempString [128];
		sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein );
		std::cout << szTempString << std::endl;
		output[0].assign( szTempString );
		output[1].assign( szTempString );

		sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/2.0) , (int)round(output[0][1]/2.0) , (int)round(output[0][2]/2.0) );
		output[2].assign( szTempString );
		
		return output; 
	}*/
	//IL13 red; HER2 green; EGFR blue
	if (pCell->phenotype.death.dead == false)
	{
		//int egfr = (int)round(0.5 * pCell->custom_data[antigen_egfr] * 255.0);
		//int il13 = (int)round(0.5 * pCell->custom_data[antigen_il13] * 255.0);
		//int her2 = (int)round(0.5 * pCell->custom_data[antigen_her2] * 255.0);
		
		 
		char szTempString[128];
		
		//sprintf(szTempString, "rgb(%u,%u,%u)", egfr, il13, her2);
		//output[1].assign(szTempString);
		/*std::cout << "EGFR: " << parameters.doubles("EGFR_antigen_level") << "; IL13: " << parameters.doubles("IL13_antigen_level") <<
			"; HER2: " << parameters.doubles("HER2_antigen_level") << "; Oncoprotein: " << pCell->custom_data[oncoprotein_i] << std::endl;*/
		if (pCell->custom_data[index_egfr] >= egfr_antigen_threshold && pCell->custom_data[index_il13] >= il13_antigen_threshold && pCell->custom_data[index_her2] >=her2_antigen_threshold )

		{
			int r = 114;  //these values come from the matlab code
			int g = 90;
			int b = 34;
			sprintf(szTempString, "rgb(%u,%u,%u)", r, g, b);
			//std::cout << szTempString << std::endl;
			output[0].assign(szTempString);
			output[1].assign(szTempString);
			output[2].assign(szTempString);
		}
		if (pCell->custom_data[index_egfr] >= egfr_antigen_threshold && pCell->custom_data[index_il13] >= il13_antigen_threshold && pCell-> custom_data[index_her2]<her2_antigen_threshold)
		{
			int r = 166;  //these values come from the matlab code
			int g = 47;
			int b = 46;
			sprintf(szTempString, "rgb(%u,%u,%u)", r, g, b);
			//std::cout << szTempString << std::endl;
			output[0].assign(szTempString);
			output[1].assign(szTempString);
			output[2].assign(szTempString);
		}
		if (pCell->custom_data[index_egfr] >= egfr_antigen_threshold && pCell->custom_data[index_her2] >= her2_antigen_threshold && pCell->custom_data[index_il13]<il13_antigen_threshold)
		{
			int r = 48;  //these values come from the matlab code
			int g = 131;
			int b = 47;
			sprintf(szTempString, "rgb(%u,%u,%u)", r, g, b); 
			//std::cout << szTempString << std::endl;
			output[0].assign(szTempString);
			output[1].assign(szTempString);
			output[2].assign(szTempString);
		}
		if (pCell->custom_data[index_il13] >= il13_antigen_threshold && pCell->custom_data[index_her2] >= her2_antigen_threshold && pCell->custom_data[index_egfr]<egfr_antigen_threshold)
		{
			int r = 129;  //these values come from the matlab code
			int g = 93;
			int b = 8;
			sprintf(szTempString, "rgb(%u,%u,%u)", r, g, b);
			//std::cout << szTempString << std::endl;
			output[0].assign(szTempString);
			output[1].assign(szTempString);
			output[2].assign(szTempString);
		}
		if (pCell->custom_data[index_egfr] >= egfr_antigen_threshold && pCell->custom_data[index_il13]<il13_antigen_threshold && pCell->custom_data[index_her2]<her2_antigen_threshold)
		{
			int r = 85;  //these values come from the matlab code
			int g = 85;
			int b = 85;
			sprintf(szTempString, "rgb(%u,%u,%u)", r, g, b);
			//std::cout << szTempString << std::endl;
			output[0].assign(szTempString);
			output[1].assign(szTempString);
			output[2].assign(szTempString);
		}
		if (pCell->custom_data[index_il13] >= il13_antigen_threshold && pCell-> custom_data[index_her2]<her2_antigen_threshold && pCell->custom_data[index_egfr]<egfr_antigen_threshold)
		{
			int r = 247;  //these values come from the matlab code
			int g = 9;
			int b = 7;
			sprintf(szTempString, "rgb(%u,%u,%u)", r, g, b);
			//std::cout << szTempString << std::endl;
			output[0].assign(szTempString);
			output[1].assign(szTempString);
			output[2].assign(szTempString);
		}
		if (pCell->custom_data[index_her2] >= her2_antigen_threshold && pCell->custom_data[index_il13]<il13_antigen_threshold && pCell->custom_data[index_egfr]<egfr_antigen_threshold)
		{
			int r = 11;  //these values come from the matlab code
			int g = 177;
			int b = 9;
			sprintf(szTempString, "rgb(%u,%u,%u)", r, g, b);
			//std::cout << szTempString << std::endl;
			output[0].assign(szTempString);
			output[1].assign(szTempString);
			output[2].assign(szTempString);
		}
		if (pCell->custom_data[index_her2]<her2_antigen_threshold && pCell->custom_data[index_il13]<il13_antigen_threshold && pCell->custom_data[index_egfr]<egfr_antigen_threshold){
			int r = 128;
			int g = 128;
			int b = 128;
			sprintf(szTempString, "rgb(%u,%u,%u)", r, g, b);
			//std::cout << szTempString << std::endl;
			output[0].assign(szTempString);
			output[1].assign(szTempString);
			output[2].assign(szTempString);

		}
		//sprintf(szTempString, "rgb(%u,%u,%u)", (int)round(output[0][0] / 2.0), (int)round(output[0][1] / 2.0), (int)round(output[0][2] / 2.0));
		//output[2].assign(szTempString);
		//std::cout << szTempString << std::endl;
		return output;
	}

	// if not, dead colors 
	
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )  // Apoptotic - Red
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
	}
	
	// Necrotic - Brown
	

	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "rgb(250,138,38)"; //orange 
		output[2] = "rgb(139,69,19)";
	}	
	
	return output; 
}

/*
void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant )
{
	std::vector<double> displacement = pAttachedTo->position - pActingOn->position; 
	axpy( &(pActingOn->velocity) , elastic_constant , displacement ); 
	
	return; 
}

void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	for( int i=0; i < pCell->state.attached_cells.size() ; i++ )
	{
		add_elastic_velocity( pCell, pCell->state.attached_cells[i], pCell->custom_data["elastic_coefficient"] ); 
	}

	return; 
}	

void attach_cells( Cell* pCell_1, Cell* pCell_2 )
{
	#pragma omp critical
	{
		
	bool already_attached = false; 
	for( int i=0 ; i < pCell_1->state.attached_cells.size() ; i++ )
	{
		if( pCell_1->state.attached_cells[i] == pCell_2 )
		{ already_attached = true; }
	}
	if( already_attached == false )
	{ pCell_1->state.attached_cells.push_back( pCell_2 ); }
	
	already_attached = false; 
	for( int i=0 ; i < pCell_2->state.attached_cells.size() ; i++ )
	{
		if( pCell_2->state.attached_cells[i] == pCell_1 )
		{ already_attached = true; }
	}
	if( already_attached == false )
	{ pCell_2->state.attached_cells.push_back( pCell_1 ); }

	}

	return; 
}

void dettach_cells( Cell* pCell_1 , Cell* pCell_2 )
{
	#pragma omp critical
	{
		bool found = false; 
		int i = 0; 
		while( !found && i < pCell_1->state.attached_cells.size() )
		{
			// if cell 2 is in cell 1's list, remove it
			if( pCell_1->state.attached_cells[i] == pCell_2 )
			{
				int n = pCell_1->state.attached_cells.size(); 
				// copy last entry to current position 
				pCell_1->state.attached_cells[i] = pCell_1->state.attached_cells[n-1]; 
				// shrink by one 
				pCell_1->state.attached_cells.pop_back(); 
				found = true; 
			}
			i++; 
		}
	
		found = false; 
		i = 0; 
		while( !found && i < pCell_2->state.attached_cells.size() )
		{
			// if cell 1 is in cell 2's list, remove it
			if( pCell_2->state.attached_cells[i] == pCell_1 )
			{
				int n = pCell_2->state.attached_cells.size(); 
				// copy last entry to current position 
				pCell_2->state.attached_cells[i] = pCell_2->state.attached_cells[n-1]; 
				// shrink by one 
				pCell_2->state.attached_cells.pop_back(); 
				found = true; 
			}
			i++; 
		}

	}
	
	return; 
}
*/

void immune_cell_motility( Cell* pCell, Phenotype& phenotype, double dt )
{
	// if attached, biased motility towards director chemoattractant 
	// otherwise, biased motility towards cargo chemoattractant 
	
	static int immune_factor_index = microenvironment.find_density_index( "immunostimulatory factor" ); 

	// if not docked, attempt biased chemotaxis 
	if( pCell->state.attached_cells.size() == 0 )
	{
		phenotype.motility.is_motile = true; 
		
		phenotype.motility.migration_bias_direction = pCell->nearest_gradient(immune_factor_index);	
		normalize( &( phenotype.motility.migration_bias_direction ) );			
	}
	else
	{
		phenotype.motility.is_motile = false; 
	}
	
	return; 
}

Cell* immune_cell_check_neighbors_for_attachment( Cell* pAttacker , double dt )
{
	std::vector<Cell*> nearby = pAttacker->cells_in_my_container(); 
	int i = 0; 
 //std::cout<< "number of nearby cells" <<nearby.size()<<std::endl;
	while( i < nearby.size() )
	{
		// don't try to kill yourself 
		if( nearby[i] != pAttacker )
		{
			if( immune_cell_attempt_attachment( pAttacker, nearby[i] , dt ) )
			{ return nearby[i]; }
		}
		i++; 
	}
	
	return NULL; 
}

bool immune_cell_attempt_attachment( Cell* pAttacker, Cell* pTarget , double dt )
{
	/*static int oncoprotein_i = pTarget->custom_data.find_variable_index("oncoprotein_type");*/
	static int attach_rate_i = pAttacker->custom_data.find_variable_index( "attachment_rate" ); 
	static int egfr_antigen_index = pTarget->custom_data.find_variable_index("EGFR_antigen_level");
	static int il13_antigen_index = pTarget->custom_data.find_variable_index("IL13_antigen_level");
	static int her2_antigen_index = pTarget->custom_data.find_variable_index("HER2_antigen_level");
	//int egfr_antigen_level = pTarget->custom_data[egfr_antigen_index];
	//int il13_antigen_level = pTarget->custom_data[il13_antigen_index];
	//int her2_antigen_level = pTarget->custom_data[her2_antigen_index];
	//std::cout << "HER2: "<<her2_antigen_level << std::endl;
	//std::cout << "IL13: "<<il13_antigen_level << std::endl;
	//std::cout << "EGFR: "<<egfr_antigen_level << std::endl;
	/*double oncoprotein_saturation =
		pAttacker->custom_data["oncoprotein_saturation"];  */
	double antigen_saturation = pAttacker->custom_data["antigen_saturation"];

	/*double oncoprotein_threshold =
		pAttacker->custom_data["oncoprotein_threshold"];   */
	double egfr_antigen_threshold = pAttacker->custom_data["egfr_antigen_threshold"];
 double her2_antigen_threshold = pAttacker->custom_data["her2_antigen_threshold"];
 double il13_antigen_threshold = pAttacker->custom_data["il13_antigen_threshold"];
	/*double oncoprotein_difference = oncoprotein_saturation - oncoprotein_threshold;*/
	double egfr_antigen_difference = antigen_saturation - egfr_antigen_threshold;
 double il13_antigen_difference = antigen_saturation - il13_antigen_threshold;
 double her2_antigen_difference = antigen_saturation - her2_antigen_threshold;

	
	double max_attachment_distance = 
		pAttacker->custom_data["max_attachment_distance"];   
	double min_attachment_distance = 
		pAttacker->custom_data["min_attachment_distance"];   
	double attachment_difference = max_attachment_distance - min_attachment_distance; 
	
	/*if (pTarget->custom_data[oncoprotein_i] > oncoprotein_threshold && pTarget->phenotype.death.dead == false)
	{
		std::vector<double> displacement = pTarget->position - pAttacker->position;
		double distance_scale = norm( displacement ); 
		if( distance_scale > max_attachment_distance )
		{ return false; } 
	
		double scale = pTarget->custom_data[oncoprotein_i];
		scale -= oncoprotein_threshold; 
		scale /= oncoprotein_difference;
		if( scale > 1.0 )
		{ scale = 1.0; } 
		
		distance_scale *= -1.0; 
		distance_scale += max_attachment_distance; 
		distance_scale /= attachment_difference; 
		if( distance_scale > 1.0 )
		{ distance_scale = 1.0; } 
		
		if( UniformRandom() < pAttacker->custom_data[attach_rate_i] * scale * dt * distance_scale )
		{
//			std::cout << "\t attach!" << " " << pTarget->custom_data[oncoprotein_i] << std::endl; 
			attach_cells( pAttacker, pTarget ); 
		}
		
		return true; 
	}*/
	if (pTarget->custom_data[egfr_antigen_index] > egfr_antigen_threshold && pTarget->phenotype.death.dead == false && pAttacker->type ==1)
	{
		//std::cout << "egfr antigen working" << std::endl;
		std::vector<double> displacement = pTarget->position - pAttacker->position;
		double distance_scale = norm(displacement);
		if (distance_scale > max_attachment_distance)
		{
			return false;
		}

		double scale = pTarget->custom_data[egfr_antigen_index];
		scale -= egfr_antigen_threshold;
		scale /= egfr_antigen_difference;
		if (scale > 1.0)
		{
			scale = 1.0;
		}

		distance_scale *= -1.0;
		distance_scale += max_attachment_distance;
		distance_scale /= attachment_difference;
		if (distance_scale > 1.0)
		{
			distance_scale = 1.0;
		}

		if (UniformRandom() < pAttacker->custom_data[attach_rate_i] * scale * dt * distance_scale)
		{
			//			std::cout << "\t attach!" << " " << pTarget->custom_data[oncoprotein_i] << std::endl; 
			attach_cells(pAttacker, pTarget);
		}

		return true;
	}

	if (pTarget->custom_data[il13_antigen_index] > il13_antigen_threshold && pTarget->phenotype.death.dead == false && pAttacker->type == 2)
	{
		//std::cout << "il13 antigen working" << std::endl;
		std::vector<double> displacement = pTarget->position - pAttacker->position;
		double distance_scale = norm(displacement);
		if (distance_scale > max_attachment_distance)
		{
			return false;
		}

		double scale = pTarget->custom_data[il13_antigen_index];
		scale -= il13_antigen_threshold;
		scale /= il13_antigen_difference;
		if (scale > 1.0)
		{
			scale = 1.0;
		}

		distance_scale *= -1.0;
		distance_scale += max_attachment_distance;
		distance_scale /= attachment_difference;
		if (distance_scale > 1.0)
		{
			distance_scale = 1.0;
		}
    

		if (UniformRandom() < pAttacker->custom_data[attach_rate_i] * scale * dt * distance_scale)
		{
			//			std::cout << "\t attach!" << " " << pTarget->custom_data[oncoprotein_i] << std::endl; 
			attach_cells(pAttacker, pTarget);
		}

		return true;
	}

	if (pTarget->custom_data[her2_antigen_index] > her2_antigen_threshold && pTarget->phenotype.death.dead == false && pAttacker->type == 3)
	{
		//std::cout << "her2 antigen working" << std::endl;
		std::vector<double> displacement = pTarget->position - pAttacker->position;
		double distance_scale = norm(displacement);
		if (distance_scale > max_attachment_distance)
		{
			return false;
		}

		double scale = pTarget->custom_data[her2_antigen_index];
		scale -= her2_antigen_threshold;
		scale /= her2_antigen_difference;
		if (scale > 1.0)
		{
			scale = 1.0;
		}

		distance_scale *= -1.0;
		distance_scale += max_attachment_distance;
		distance_scale /= attachment_difference;
		if (distance_scale > 1.0)
		{
			distance_scale = 1.0;
		}

		if (UniformRandom() < pAttacker->custom_data[attach_rate_i] * scale * dt * distance_scale)
		{
			//			std::cout << "\t attach!" << " " << pTarget->custom_data[oncoprotein_i] << std::endl; 
			attach_cells(pAttacker, pTarget);
		}

		return true;
	}
	
	return false; 
}

bool immune_cell_attempt_apoptosis( Cell* pAttacker, Cell* pTarget, double dt )
{
	static int oncoprotein_i = pTarget->custom_data.find_variable_index( "oncoprotein" ); 
	static int egfr_antigen_index = pTarget->custom_data.find_variable_index("EGFR_antigen_level");
	static int il13_antigen_index = pTarget->custom_data.find_variable_index("IL13_antigen_level");
	static int her2_antigen_index = pTarget->custom_data.find_variable_index("HER2_antigen_level");
	static int apoptosis_model_index = pTarget->phenotype.death.find_death_model_index( "apoptosis" );	
	static int kill_rate_index = pAttacker->custom_data.find_variable_index( "kill_rate" ); 
	
	double antigen_saturation = pAttacker->custom_data["antigen_saturation"];
	double egfr_antigen_threshold = pAttacker->custom_data["egfr_antigen_threshold"];
 double il13_antigen_threshold = pAttacker->custom_data["il13_antigen_threshold"];
 double her2_antigen_threshold = pAttacker->custom_data["her2_antigen_threshold"];
	double egfr_antigen_difference = antigen_saturation - egfr_antigen_threshold;
 double il13_antigen_difference = antigen_saturation - il13_antigen_threshold;
 double her2_antigen_difference = antigen_saturation - her2_antigen_threshold;
	double oncoprotein_saturation =
		pAttacker->custom_data["oncoprotein_saturation"]; // 2.0; 
	double oncoprotein_threshold =  
		pAttacker->custom_data["oncoprotein_threshold"]; // 0.5; // 0.1; 
	double oncoprotein_difference = oncoprotein_saturation - oncoprotein_threshold;
	

	// new 
	/*if (pTarget->custom_data[oncoprotein_i] < oncoprotein_threshold)
	{ return false; }
	double scale = pTarget->custom_data[oncoprotein_i];
	scale -= oncoprotein_threshold;
	scale /= oncoprotein_difference;
	if (scale > 1.0)
	{
		scale = 1.0;
	}*/
	
	if (pTarget->custom_data[egfr_antigen_index]<egfr_antigen_threshold && pTarget->custom_data[il13_antigen_index] < il13_antigen_threshold &&
		pTarget->custom_data[her2_antigen_index] < her2_antigen_threshold)
	{
		return false;
	}
	// new 
	double scale_egfr = pTarget->custom_data[egfr_antigen_index];
	scale_egfr -= egfr_antigen_threshold; 
	scale_egfr /= egfr_antigen_difference;

	if( scale_egfr > 1.0 )
	{ scale_egfr = 1.0; } 

	double scale_il13 = pTarget->custom_data[il13_antigen_index];
	scale_il13 -= il13_antigen_threshold;
	scale_il13 /= il13_antigen_difference;

	if (scale_il13 > 1.0)
	{
		scale_il13 = 1.0;
	}

	double scale_her2 = pTarget->custom_data[her2_antigen_index];
	scale_her2 -= her2_antigen_threshold;
	scale_her2 /= her2_antigen_difference;

	if (scale_her2 > 1.0)
	{
		scale_her2 = 1.0;
	}
	
	/*if (UniformRandom() < pAttacker->custom_data[kill_rate_index] * scale * dt)
	{ 
//		std::cout << "\t\t kill!" << " " << pTarget->custom_data[oncoprotein_i] << std::endl; 
		return true; 
	}*/
	//std::cout << "condition to kill (egfr): " << pAttacker->custom_data[kill_rate_index] * scale_egfr * dt << std::endl;
	//std::cout << "condition to kill (il13): " << pAttacker->custom_data[kill_rate_index] * scale_il13 * dt << std::endl;
	//std::cout << "condition to kill (her2): " << pAttacker->custom_data[kill_rate_index] * scale_her2 * dt << std::endl;
	//std::cout << "condition to kill (by oncoprotein): " << pAttacker->custom_data[kill_rate_index] * scale * dt << std::endl;
	if (UniformRandom() < pAttacker->custom_data[kill_rate_index] * scale_egfr * dt && pAttacker->type == 1)
	{
		//		std::cout << "\t\t kill!" << " " << pTarget->custom_data[oncoprotein_i] << std::endl; 
		//std::cout << "egfr killing" << std::endl;
		return true;
	}
	if (UniformRandom() < pAttacker->custom_data[kill_rate_index] * scale_il13 * dt && pAttacker->type==2)
	{
		//		std::cout << "\t\t kill!" << " " << pTarget->custom_data[oncoprotein_i] << std::endl; 
		//std::cout << "il13 killing" << std::endl;
		return true;
	}
	if (UniformRandom() < pAttacker->custom_data[kill_rate_index] * scale_her2 * dt && pAttacker->type ==3)
	{
		//		std::cout << "\t\t kill!" << " " << pTarget->custom_data[oncoprotein_i] << std::endl; 
		//std::cout << "her2 killing" << std::endl;

		return true;
	}
	return false; 
}

bool immune_cell_trigger_apoptosis( Cell* pAttacker, Cell* pTarget )
{
	static int apoptosis_model_index = pTarget->phenotype.death.find_death_model_index( "apoptosis" );	
	
	// if the Target cell is already dead, don't bother!
	if( pTarget->phenotype.death.dead == true )
	{ return false; }

	pTarget->start_death( apoptosis_model_index );
	return true; 
}

void immune_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int attach_lifetime_i = pCell->custom_data.find_variable_index( "attachment_lifetime" ); 
	
	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions, 
		// since those are part of mechanics. 
		
		// Let's just fully disable now. 
		pCell->functions.custom_cell_rule = NULL; 
		return; 
	}
	
	// if I'm docked
	if( pCell->state.number_of_attached_cells() > 0 )
	{
		// attempt to kill my attached cell
		
		bool detach_me = false; 
		
		if( immune_cell_attempt_apoptosis( pCell, pCell->state.attached_cells[0], dt ) )
		{
			immune_cell_trigger_apoptosis( pCell, pCell->state.attached_cells[0] ); 
			detach_me = true; 
		}
		
		// decide whether to detach 
		
		if( UniformRandom() < dt / ( pCell->custom_data[attach_lifetime_i] + 1e-15 ) )
		{ detach_me = true; }
		
		// if I dettach, resume motile behavior 
		
		if( detach_me )
		{
			detach_cells( pCell, pCell->state.attached_cells[0] ); 
			phenotype.motility.is_motile = true; 
		}
		return; 
	}
	
	// I'm not docked, look for cells nearby and try to docked
	
	// if this returns non-NULL, we're now attached to a cell 
	if( immune_cell_check_neighbors_for_attachment( pCell , dt) )
	{
		// set motility off 
		phenotype.motility.is_motile = false; 
		return; 
	}
	phenotype.motility.is_motile = true; 
	
	return; 
}

void adhesion_contact_function( Cell* pActingOn, Phenotype& pao, Cell* pAttachedTo, Phenotype& pat , double dt )
{
	std::vector<double> displacement = pAttachedTo->position - pActingOn->position; 
	
	static double max_elastic_displacement = pao.geometry.radius * pao.mechanics.relative_detachment_distance; 
	static double max_displacement_squared = max_elastic_displacement*max_elastic_displacement; 
	
	// detach cells if too far apart 
	
	if( norm_squared( displacement ) > max_displacement_squared )
	{
		detach_cells( pActingOn , pAttachedTo );
		return; 
	}
	
	axpy( &(pActingOn->velocity) , pao.mechanics.attachment_elastic_constant , displacement ); 
	
	return; 
}

