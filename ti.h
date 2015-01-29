
/** @file
 * @author James W. Barnett jbarnet4@tulane.edu
 * @date January 29, 2015
 * @brief Header file for ti.cpp
 */

#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <vector>

#include "boost/program_options/cmdline.hpp"
#include "boost/program_options/config.hpp"
#include "boost/program_options/environment_iterator.hpp"
#include "boost/program_options/eof_iterator.hpp"
#include "boost/program_options/errors.hpp"
#include "boost/program_options/option.hpp"
#include "boost/program_options/options_description.hpp"
#include "boost/program_options/parsers.hpp"
#include "boost/program_options/positional_options.hpp"
#include "boost/program_options/value_semantic.hpp"
#include "boost/program_options/variables_map.hpp"
#include "boost/program_options/version.hpp"

#include "gromacs/fileio/enxio.h"
#include "gromacs/utility/smalloc.h"

using namespace std;

/**
 * @brief Class containing dV/dl and lambda values.
 *
 * @details Contains all relevant info for thermodynamic integration to occur.
 */

class dVdl_energy_info {
	private:
		/**
		 * @brief The number of dV/dl groups in the energy files from the simulation.
		 */
		int n_groups;

		/**
		 * @brief The number of lambda values.
		 */
		int n_lambdas;

		/**
		 * @brief Routine for saving dVdl values from energy files.
		 * @param filename Array of names of energy files to read.
		 * @parap skip Number of saved frames to skip between data points.
		 */
		int save_dVdl(char** filename, int skip);

		/**
		 * @brief The units (kJ/mol or kcal/mol) of dV/dl
		 */
		string units;


		/**
		 * @brief Numeric locations of dV/dl energy groups in energy files.
		 */
        vector <int> locations;


		/**
		 * @brief Names of dV/dl groups, from energy files.
		 */
        vector <string> group_names;

		/**
		 * @brief Lambda values stored from the toplogy file.
		 * @details The first dimension is the group (i.e. corresponds with vdw-lambdas, etc.) 
		 * The second is the location within that group.
		 */
		vector < vector <double> > lambda;

		/**
		 * @brief All dV/dl values from all energy files.
		 * @details The first dimension is the group (i.e. corresponds with vdw-lambdas, etc.),
		 * the second is the lambda number, and the third is the data as a time series. So, dV/dl
		 * values are stored by lambda and energy group.
		 */
		vector < vector < vector <double> > > dVdl;

		int get_group_info(ener_file *fp);

    public:

		dVdl_energy_info();

		/** 
		 * @brief Construction which reads in all info from energy and topology files.
		 * @details If no tpr files are specified all relevant info is attempted to be read from the energy files.
		 * However, if dhdl-separate-file is set to 'yes' (the default), it's necessary to gather info from
		 * the tpr files.
		 * @param edrfiles Vector of energy files
		 */
		dVdl_energy_info(vector <string> edrfiles);

		/** @brief Gets the dV/dl value for a specific group, lambda, and frame.
		 * @param group Group number.
		 * @param lambda Lambda number.
		 * @param i Location in vector.
		 * @return dV/dl
		 */
		double get_dVdl(int group, int lambda, int i);

		/** @brief Gets the lambda for a group and specific simulation.
		 * @param group Group number, corresponding with vdw-lambdas, etc.
		 * @param i Number corresponding with simulation.
		 * @return lambda
		 */
		double get_lambda(int group,int i);

		/** @brief Gets the number of dV/dl values for a specific lambda and group.
		 * @param group Group number, corresponding with vdw-lambdas, etc.
		 * @param lambda Lambda number.
		 * @return Number of dV/dl values.
		 */
		int get_dVdl_size(int group, int lambda);

		/** @brief Gets the energy file location for a saved group.
		 * @details Since only a couple of groups are saved (the dV/dl groups) from the energy file
		 * we can get the location in the energy file by giving our saved group number. That is,
		 * if dVvdw/dl is saved in the 2nd spot in all of our vectors, we do get_group_location(1), which will
		 * tell us the numeric value of where it is located in the energy file (maybe 12, for example).
		 * @param i Group number.
		 * @return Energy file group number.
		 */
        int get_group_location(int i);

		/** @brief Get the number of dV/dl groups that were saved.
		 * @return Number of saved dV/dl groups.
		 */
		int get_n_groups();
		
		/** @brief Get the number of lambdas that were saved.
		 * @return Number of lambdas.
		 */
		int get_n_lambdas();

		/** @brief Sorts edr files by lambda stored in tpr files
		 * @param edrfiles Array of energy files
		 * @param tprfiles Array of tpr files, in same order as energy files
		 * @param nedrfile Number of energy files
		 * @return 0 on success
		 */
		int sort_files(char** edrfiles, char** tprfiles, int nedrfile);

		/** @brief Get the name of a saved energy group.
		 * @param i Group number.
		 * @return Name of saved energy group (i.e. dVvdw/dl).
		 */
        string get_group_name(int i);

		/** @brief Gets the units for the free energy.
		 */
		string get_units();

		/** @brief Get all the lambdas for a specific group.
		 * @details This will return all of the lambdas for a specific group (i.e., vdw-lambdas).
		 * @return Vector of lambdas for specified group.
		 */
		vector <double> get_lambdas(int group);

};

/**
 * @brief Performs integration based on predetermined weights.
 * @param data The data to be integrated.
 * @param weights The weights for each interval.
 * @return The summation result.
 */
double integrate(vector <double> data, vector <double> weights);

/**
 * @brief Gets the weights for Gaussian-Legendre integration.
 * @param lambda The set of lambdas for this integration.
 * @return The weights, used in integrate function.
 */
vector <double> get_gaussian_quadurature_weights(vector <double> lambda);

/**
 * @brief Prints the points of Gaussian-Legendre integration.
 * @param npts Number of points to be used in integration.
 * @return 0 on sucess.
 */
int get_gaussian_quad_points(int npts);

/**
 * @brief The interior portion of Gaussian-Legendre weight/point calculation.
 * @param i Iteration number.
 * @param n Number of points.
 * @param t Returned.
 * @param pp Returned.
 */
int do_legendre_calc(int i, int n, double &t, double &pp);

/**
 * @brief Gets integration weights for Simpson's rule.
 * @param lambda The set of lambdas for this integration.
 * @return The weights, used in integrate function.
 */
vector <double> get_simpsons_weights(vector <double> lambda);

/**
 * @brief Gets integration weights for the trapezoid rule.
 * @param lambda The set of lambdas for this integration.
 * @return The weights, used in integrate function.
 */
vector <double> get_trapezoid_weights(vector <double> lambda);

vector <double> get_weights(vector <double> lambdas, string integration_type);

/**
 * @brief Gets lambdas for Gaussian quadrature
 * @details Calculates points (the lambdas) for Gaussian-Legendre quadrature for the range 0.0 to
 * 1.0 and prints to screen.
 * @param Number of points.
 * @return 0 on success.
 */
int get_free_energy_change(dVdl_energy_info &dei, double &free_energy_change);

#endif
