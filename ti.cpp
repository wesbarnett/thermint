
/**
 * @file
 *
 * @author James W. Barnett jbarnet4@tulane.edu
 * @date January 10, 2015
 * @brief Analysis module for performing thermodynamic integration
 * @see ti.h
 *
 */

#include "ti.h"

const double bootstrap_tol = 1.0e-1;
const vector <string> dV_group_names = 
{
    "dVvdw/dl",
    "dVcoul/dl",
    "dVbonded/dl",
    "dVrestraint/dl"
};
const vector <string> lambda_group_names = 
{
    "vdw-lambdas",
    "coul-lambdas",
    "bonded-lambdas",
    "restraint-lambdas"
};
const vector <string> group_short_names =
{
    "vdw",
    "coul",
    "bonded",
    "restraint"
};

int main(int argc, char* argv[]) {

    double bin;
    double bin_size;
    double bootstrap_avg = 0.0;
    double bootstrap_var = 0.0;
    double bootstrap_min;
    double bootstrap_max;
    double bootstrap_result;
    double mean;
    double result = 0.0;
    double sum = 0.0;
    double sum2 = 0.0;
    double uncertainty = 0.0;
    double var;
    double x;
    int block;
    int block_count;
    int block_n = 5;
    int bootstrap_n = 1000;
    int bins_n = 200;
    int i;
    int j;
    int group_count = 0;
    int lambda_count;
    int N;
    int pos;
    ofstream oFS;
	string integration_type;
    string bootstrap_file;
    vector <double> bootstrap_values;
    vector <double> bootstrap_hist;
    vector <double> dVdlAvg;
    vector <double> dVdlVar;
    vector <double> weights;
	vector <string> energy_file_names;

	namespace po = boost::program_options;

	string description = "See README.md for more details. Usage";

	po::options_description desc(description);
	po::positional_options_description p;
	p.add("file", -1);
	desc.add_options()
		("help,h","Show usage.")
		("file,f",po::value< vector <string> >(&energy_file_names),"Energy files to be read in.")
		("method,m",po::value<string>(&integration_type)->default_value("trapezoid"),"Method: trapezoid, simpsons, or gaussian.")
		("output,o",po::value<string>(&bootstrap_file)->default_value("bootstrap.dat"),"Name of bootstrap histogram file to be output.")
		("nboot",po::value<int>(&bootstrap_n)->default_value(200),"Number of iterations to be performed in bootstrap calculation.")
		("nblocks",po::value<int>(&block_n)->default_value(5),"Number of blocks for bootstrap calculation.")
		("nbins",po::value<int>(&bins_n)->default_value(200),"Number of bins for bootstrap histograms.")
	;

	po::variables_map vm;
	po::store(po::command_line_parser(argc,argv).options(desc).positional(p).run(),vm);
	po::notify(vm);

	if (vm.count("h"))
	{
		cout << desc << endl;
		return -1;
	}

	if (!vm.count("file"))
	{
		cout << "ERROR: You must specifiy energy files." << endl;
		cout << desc << endl;
		return -1;
	}

    bootstrap_hist.resize(bins_n+1,0.0);

    if (energy_file_names.size() == 1)
    {
        cout << "ERROR: Only one energy file specified! You need to specify all energy files from all simulations for this free energy change calculation." << endl;
		return -1;
    }


    switch (integration_type.at(0)) 
    {
        case 'g':
			cout << "Using Gaussian-Legendre quadrature." << endl;
	        cout << "NOTE: You selected Gaussian quadrature, so the lambdas you used for your simulation should match Gaussian-Legendre points." << endl;
            break;
        case 's':
            cout << "Using Simpons's rule." << endl;
            break;
        case 't':
			cout << "Using the trapezoid rule." << endl;
            break;
        default:
			cout << "Error: not a valid integration type." << endl;
			return -1;
            break;
    }

	cout << endl;

    srand(time(0));

    try
    {
    dVdl_energy_info dei(energy_file_names);

	cout << "Results in " << dei.get_units() << endl << endl;

    for (group_count = 0; group_count < dei.get_n_groups(); group_count++) 
    {

        pos = find(dV_group_names.begin(), dV_group_names.end(), dei.get_group_name(group_count)) - dV_group_names.begin();
		cout << group_short_names.at(pos) << ":" << endl;
        try
        {
            weights = get_weights(dei.get_lambdas(group_count),integration_type);
        }
        catch(exception const& e)
        {
            cout << "ERROR: " << e.what() << "\n";
            return -1;
        }

		cout << right;
		cout << setw(3) << "Sim";
        cout << setw(10) << "lambda";
		cout << "  <" << setw(8) << dei.get_group_name(group_count) << "> ± " << setw(8) << " std.dev.";
        cout << setw(12) << "weight";
        cout << setw(5) << "w*<" << dei.get_group_name(group_count) << ">";
        cout << setw(9) << "tot" << endl;

        for (lambda_count = 0; lambda_count < dei.get_n_lambdas(); lambda_count++) 
        {
            sum = 0.0;
            N = dei.get_dVdl_size(group_count,lambda_count);
            for (j = 0; j < N; j++) 
            {
                sum += dei.get_dVdl(group_count,lambda_count,j);
            }
            mean = sum / N;

            var = 0.0;
            sum2 = 0.0;
            for (j = 0; j < N; j++) 
            {
                x = dei.get_dVdl(group_count,lambda_count,j);
                sum2 += (x - mean) * (x - mean);
            }
            var = sum2 / (N-1);

            result += mean*weights.at(lambda_count);

			cout << fixed << setprecision(0);
			cout << setw(3) << lambda_count+1;
			cout << fixed << setprecision(3);
			cout << setw(10) << dei.get_lambda(group_count,lambda_count);
			cout << setw(12) << mean;
			cout << " ± " << setw(8) << sqrt(var);
			cout << setw(12) << weights.at(lambda_count);
			cout << setw(12) << mean*weights.at(lambda_count);
			cout << setw(12) << result << endl;
        }
		cout << endl;
    }

    cerr << "Doing bootstrap calculation for uncertainty (could take a few minutes)..." << endl << endl;

    /* Bootstrap section starts here.
     * Block bootstrapping is used to help get an accurate uncertainty, since
     * there could be (probably is) correlated dV/dl's in each simulation.
     */

    bootstrap_avg = 0.0;
    bootstrap_var = 0.0;
    bootstrap_values.resize(bootstrap_n);
    #pragma omp parallel for private(i,bootstrap_result,group_count,weights,lambda_count,sum,N,block_count,block,j,mean)
    for (i = 0; i < bootstrap_n; i++) 
    {

        bootstrap_result = 0.0;
        for (group_count = 0; group_count < dei.get_n_groups(); group_count++) 
        {

            weights = get_weights(dei.get_lambdas(group_count),integration_type);

            for (lambda_count = 0; lambda_count < dei.get_n_lambdas(); lambda_count++) 
            {

                sum = 0.0;
                N = dei.get_dVdl_size(group_count,lambda_count) / block_n;
                for (block_count = 0; block_count < block_n; block_count++)
                {
                    block = rand() % block_n;
                    for (j = 0+N*block; j < N*(block+1); j++) 
                    {
                        sum += dei.get_dVdl(group_count,lambda_count,j);
                    }

                }
                mean = sum / dei.get_dVdl_size(group_count,lambda_count);

                bootstrap_result += mean * weights.at(lambda_count);

            }

        }

        bootstrap_values.at(i) = bootstrap_result;
        bootstrap_avg += bootstrap_result;
    }

    bootstrap_avg /= (bootstrap_n-1.0);

    if (abs(bootstrap_avg - result) > bootstrap_tol) {
        cout << endl << "WARNING: Bootstrap result differs greatly from original result (>" << bootstrap_tol << "), so the calculated uncertainty may not ";
        cout << "be accurate. You should rerun with more iterations using -nboot flag." << endl;
    }

    for (i = 0; i < bootstrap_n; i++) 
    {
        bootstrap_var += pow(bootstrap_values.at(i)-bootstrap_avg,2);
    }

    bootstrap_var /= (bootstrap_n-1.0);
    uncertainty = sqrt(bootstrap_var);

    bootstrap_min = *min_element(bootstrap_values.begin(), bootstrap_values.end());
    bootstrap_max = *max_element(bootstrap_values.begin(), bootstrap_values.end());
        
    bin_size = abs(bootstrap_max - bootstrap_min) / bins_n;
    for (i = 0; i < bootstrap_n; i++)
    {
        bin = ceil((bootstrap_values.at(i) - bootstrap_min) / bin_size);
        if (bin < bootstrap_hist.size()) {
            bootstrap_hist.at(bin) += 1.0;
        }
    }

    for (i = 0; i < bins_n; i++)
    {
        bootstrap_hist.at(i) /= bootstrap_n;
    }

	cout << "total: " << result << " ± " << uncertainty << " " << dei.get_units() << endl << endl;

    oFS.open(bootstrap_file.c_str());
    oFS << setprecision(6);
    oFS << "# Bootstrapped Free Energy Difference Distribution" << endl;
    oFS << "# X: Free energy difference" << endl;
    oFS << "# Y: Probability" << endl;
    for (i = 0; i < bins_n; i++)
    {
        oFS << bootstrap_min+i*bin_size << "  " << bootstrap_hist.at(i) << endl;
    }
    oFS.close();

    }

    catch(exception const& e)
    {
        cout << "ERROR: " << e.what() << "\n";
        return -1;
    }

    return 0;
}

vector <double> get_trapezoid_weights(vector <double> lambda) {

    vector <double> w(lambda.size());
    int i;

    for (i = 0; i < lambda.size(); i++)   
    {
        if (i == 0) 
        {
            w.at(i) = (lambda.at(i+1) - lambda.at(i)) / 2.0;
        } 
        else if (i == lambda.size()-1)
        {
            w.at(i) = (lambda.at(i) - lambda.at(i-1)) / 2.0;
        }
        else
        {
            w.at(i) = (lambda.at(i+1) - lambda.at(i-1)) / 2.0;
        }
    }

    return w;
}

vector <double> get_simpsons_weights(vector <double> lambda) {

    vector <double> w(lambda.size());
    int i;

    if (lambda.size() % 2 == 0) 
    {
        throw runtime_error("Simpson's rule can only be used with an odd number of points! Choose trapezoid method and try again.");
    }

    for (i = 0; i < lambda.size(); i++)   
    {
        if (i == 0) 
        {
            w.at(i) = (lambda.at(i+1) - lambda.at(i)) / 3.0;
        } 
        else if (i == lambda.size() - 1)
        {
            w.at(i) = (lambda.at(i) - lambda.at(i-1)) / 3.0;
        }
        else if (i % 2 == 1) 
        {
            w.at(i) = (lambda.at(i+1) - lambda.at(i-1)) * 2.0 / 3.0;
        } 
        else 
        {
            w.at(i) = (lambda.at(i+1) - lambda.at(i-1)) * 1.0 / 3.0;
        }
    }

    return w;
}

vector <double> get_gaussian_quadurature_weights(vector <double> lambda) 
{

    int m;
    int i;
    double t;
    double pp;

    int n = lambda.size();

    vector <double> w(n);

    m = (n + 1) / 2;

    for (i = 1; i<= m; i++) 
    {  
        do_legendre_calc(i,n,t,pp);
        w.at(i-1) = 2.0/((1.0-t*t)*pp*pp);
        w.at(n-i) = w.at(i-1.0);
    }

    for (i = 0; i < n ; i++)  
    {
        w.at(i) = w.at(i) * 0.5;
    }

	return w;
}

int get_gaussian_quad_points(int n) 
{

    int i;
    int m;
    double t;
    double pp;

    vector <double> x(n);

    m = (n + 1) / 2;

    for (i = 1; i <= m; i++) 
    {  
        do_legendre_calc(i,n,t,pp);
        x.at(i-1) = -t;
        x.at(n-i) = t;
    }
    for (i = 0; i < n; i++) 
    {
        x.at(i) = (x.at(i) + 1.0) * 0.5;
    }

    cout << "Lambdas for " << n << " point Gaussian-Legendre quadrature:" << endl << endl;
    for (i = 0; i < n; i++) 
    {
		cout << x.at(i) << endl;
    }
	cout << endl << endl;
    cout << "Now place the values in the lambdas parameter in your .mdp files (e.g., vdw-lambdas, etc.)." << endl;

	return 0;
}

int do_legendre_calc(int i, int n, double &t, double &pp)
{
    const double eps= 3.0e-14;
	double t1;
	double p1;
	double p2;
	double p3;
    int j;

    t = cos( M_PI * (i - 0.25) / (n + 0.5) );
    t1 = 1.0;

    while(abs(t-t1) > eps) 
    { 

        p1 = 1.0;
        p2 = 0.0;

        for (j = 1; j <= n; j++) 
        {

            p3 = p2;
            p2 = p1;
            p1 = ((2.0*j-1.0)*t*p2-(j-1.0)*p3)/j;

        }

        pp = n*(t*p1-p2)/(t*t-1.0);
        t1 = t;
        t = t1 - p1/pp;

    }   

    return 0;

}

vector <double> get_weights(vector <double> lambdas, string integration_type)
{
    if (integration_type.at(0) == 'g') return get_gaussian_quadurature_weights(lambdas);
    if (integration_type.at(0) == 's') return get_simpsons_weights(lambdas);
    return get_trapezoid_weights(lambdas);
}

/*
 * dVdl_energy_info class methods start here. See gmx_ti.h for more info.
 */
dVdl_energy_info::dVdl_energy_info(vector <string> edrfiles) 
{

    ener_file_t fp;
    t_enxframe *fr;
    double this_lambda;
    int fep_state = 0;
    int count;
    int i;
    int j;
    int k;
    int pos;
    gmx_enxnm_t *nm = NULL;
    int nre;

    /* This is the preferred way to read in lambda values, directly from the
     * energy files. However, by default lambdas are not saved there (see
     * below), so we may have to get them from the tpr files).
     */
	n_lambdas = edrfiles.size();

    for (j = 0; j < edrfiles.size(); j++)
    {

		cerr << endl;
        fp = open_enx(edrfiles.at(j).c_str(),"r");

        // Get energy groups with dV/dl in first energy file
		if (j == 0)
        {
			get_group_info(fp);
            dVdl.resize(n_groups);
            lambda.resize(n_groups);
        } 
        else
        {
			do_enxnms(fp,&nre,&nm);
        }

        snew(fr,1);
        do_enx(fp,fr);
        if (fr->nblock == 0)
        {
            throw runtime_error("Energy file does not have needed values. Most likely separate-dhdl-file was set to yes in the .mdp file.");
        }
        // Get lambda from subblock
        for (i = 0; i < fr->nblock; i++)
        {
			if (fr->block[i].id == enxDHCOLL)
            {
				cerr << "lambdas:" << endl;
                fep_state = fr->block[i].sub[1].ival[0];
                for (k = 0; k < n_groups; k++)
                {
					lambda.at(k).resize(edrfiles.size());
                    this_lambda = fr->block[i].sub[0].dval[5+k];
                    pos = find(dV_group_names.begin(), dV_group_names.end(), group_names.at(k)) - dV_group_names.begin();
					cerr <<" " << group_short_names.at(pos) << ": " << this_lambda << endl;
                    lambda.at(k).at(fep_state) = this_lambda;
                }
			}
        }

        // Save data from first frame since we already read it
        for (i = 0; i < n_groups; i++)
        {
			dVdl.at(i).resize(edrfiles.size());
            dVdl.at(i).at(fep_state).push_back(fr->ener[locations.at(i)].e);
        }

        // Now go through all of the frames
        count = 0;
        while (do_enx(fp,fr)) 
        {
			for (i = 0; i < n_groups; i++)
            {
				dVdl.at(i).at(fep_state).push_back(fr->ener[locations.at(i)].e);
            }
        count++;
        }
		cerr << endl;

        free_enxframe(fr);
        sfree(fr);

        close_enx(fp);
    }

}

int dVdl_energy_info::get_group_info(ener_file *fp)
{
    gmx_enxnm_t *nm = NULL;
    int nre;

    do_enxnms(fp,&nre,&nm);

    // Filter the groups until we get the dV/dl groups
    for (int i = 0; i < nre; i++) 
    {
        string mystring(nm[i].name);
        if (find(dV_group_names.begin(), dV_group_names.end(),mystring) != dV_group_names.end())
        { 
            group_names.push_back(mystring);
            locations.push_back(i);
            units = nm[i].unit;
        }
    }

    n_groups = group_names.size();

    if (n_groups > 0) 
    {
        cerr << "the following dv/dl energy groups were found: " << endl;
        for (int i =0; i < n_groups; i++) 
        {
			cerr <<"  " << group_names.at(i) << endl;
        }
		cerr << endl;
        return 0;
    } 
    else 
    {
        throw runtime_error("No groups with dV/dl present in the energy file.");
    }

}


string dVdl_energy_info::get_group_name(int i) 
{
    return group_names.at(i);
}

int dVdl_energy_info::get_group_location(int i) 
{
    return locations.at(i);
}

int dVdl_energy_info::get_n_groups() 
{
    return n_groups;
}

int dVdl_energy_info::get_n_lambdas() 
{
    return n_lambdas;
}

double dVdl_energy_info::get_dVdl(int group, int lambda, int i) 
{
    return dVdl.at(group).at(lambda).at(i);

}

int dVdl_energy_info::get_dVdl_size(int group, int lambda) 
{
    return dVdl.at(group).at(lambda).size();

}

double dVdl_energy_info::get_lambda(int group, int i) 
{
    return lambda.at(group).at(i);
}

vector <double> dVdl_energy_info::get_lambdas(int group)
{
    return lambda.at(group);
}

string dVdl_energy_info::get_units()
{
    return units;
}
